# HOPG（Highly Ordered Pyrolytic Graphite）X線分光器データ解析プログラム
# 事前キャリブレーション済みパラメータを使用してスペクトラム解析を実行

import pprint
import numpy as np
import re
import os
import tkinter as tk
from tkinter import filedialog, messagebox
from time import sleep
from scipy.optimize import differential_evolution
from scipy import integrate
from scipy.interpolate import interp1d
from numpy.linalg import solve
import matplotlib.pyplot as plt
import csv
import requests
from bs4 import BeautifulSoup
from datetime import datetime, timedelta
import urllib.parse
import numpy as np
import re
import os
import tkinter as tk
from tkinter import filedialog, messagebox
from time import sleep
from scipy.optimize import differential_evolution
from scipy import integrate
from scipy.interpolate import interp1d
from numpy.linalg import solve
import matplotlib.pyplot as plt
import csv

# 時間計算関数をenergy_conversion_HOPG_calibration.pyからインポート
try:
    from energy_conversion_HOPG_calibration import (
        extract_reading_time_from_filename,
        fetch_shot_time_from_web,
        calculate_time_delay_auto
    )
    AUTO_TIME_AVAILABLE = True
except ImportError:
    AUTO_TIME_AVAILABLE = False

# HOPGの格子間隔 [Å]
D_HOPG = 3.357 

def ip_time_correction(t):
    """時間に依存するIPの感度補正関数
    
    Args:
        t (float): 測定時間遅延 [時間]
    
    Returns:
        float: 時間補正係数
    """
    return 0.297 * np.exp(-t/57.6) + 0.723

def ip_energy_correction(energy):
    """エネルギーに依存するIPの感度補正関数
    
    Args:
        energy (float or array): X線エネルギー [keV]
    
    Returns:
        float or array: エネルギー補正係数
    """
    m = 0.17   # 傾き
    c = 2.27   # 切片
    return m * energy + c

def load_filter_data():
    """フィルターデータの読み込みと補間関数作成
    
    Returns:
        tuple: (filter1, filter2) - ベリリウムとポリエチレンフィルターの補間関数
    """
    # Be（ベリリウム）フィルター 500μm厚の透過率データ
    filter1_data = np.loadtxt('Be_500um_transmittance.dat', skiprows=2, unpack=True)
    filter1 = interp1d(filter1_data[0,:]*1.0E-03, filter1_data[1,:], kind='cubic')
    
    # CH2（ポリエチレン）フィルター 50μm厚の透過率データ  
    filter2_data = np.loadtxt('CH2_50um_transmittance.dat', skiprows=2, unpack=True)
    filter2 = interp1d(filter2_data[0,:]*1.0E-03, filter2_data[1,:], kind='cubic')
    
    return filter1, filter2

def get_user_input(file_path=None, shot_num=None, laser_type=None):
    """ユーザーからの入力取得
    
    Args:
        file_path (str, optional): データファイルパス（自動計算用）
        shot_num (str, optional): ショット番号（自動計算用）
        laser_type (str, optional): レーザータイプ（自動計算用）
    
    Returns:
        float: 測定時間遅延 [時間]
    """
    if AUTO_TIME_AVAILABLE and file_path and shot_num and laser_type:
        print("測定時間遅延の入力方法を選択してください:")
        print("1. 自動計算（ファイル名とWebページから）")
        print("2. 手動入力")
        
        while True:
            choice = input("選択 (1 または 2): ").strip()
            if choice == "1":
                print("   - 自動時間遅延計算を開始します...")
                time_delay = calculate_time_delay_auto(file_path, shot_num, laser_type)
                
                if time_delay is None:
                    print("   - 自動計算に失敗しました。手動入力に切り替えます。")
                    t = input("   手動で時間遅延を入力してください（時間）: ")
                    return float(t)
                else:
                    print("   - 自動計算完了: {:.3f} 時間".format(time_delay))
                    return time_delay
            elif choice == "2":
                t = input("Input the scanning time delay (hours): ")
                return float(t)
            else:
                print("1 または 2 を入力してください。")
    else:
        # 自動計算が利用できない場合は手動入力のみ
        if not AUTO_TIME_AVAILABLE:
            print("注意: 自動時間計算機能が利用できません（キャリブレーションスクリプトが必要）")
        t = input("Input the scanning time delay (hours): ")
        return float(t)

def load_experimental_data(shot_num, file_path):
    """実験データの読み込み
    
    Args:
        shot_num (str): ショット番号
        file_path (str): データファイルパス
    
    Returns:
        numpy.ndarray: IPスキャンデータ
    """
    return np.genfromtxt('{}.csv'.format(file_path), delimiter=',', 
                        skip_header=1, skip_footer=0, autostrip=True)

def calibrate_parameters():
    """キャリブレーションパラメータの設定と微調整
    
    Returns:
        tuple: (s, E_def, theta_def, place) - キャリブレーションパラメータ
    """
    # 基準線のエネルギー [keV]
    E_def = 8.048
    # 基準線のブラッグ角 [rad]
    theta_def = np.arcsin(12.3984/(2.0 * D_HOPG * E_def))
    # 基準線のIP上での位置 [cm]
    place = 1.005  # 手動で設定された基準線位置
    
    # 事前に決定されたキャリブレーションパラメータ
    s = [-0.01988027, -0.79709931]
    
    # キャリブレーションパラメータの微調整
    diff = 13 + s[0] - place - (50 - s[1]) * np.tan(np.arcsin(12.3984/(2*D_HOPG*E_def)))
    s[0] = s[0] - diff
    
    return s, E_def, theta_def, place

def convert_position_to_energy(data, s):
    """IP上の位置座標からエネルギーへの変換
    
    Args:
        data (numpy.ndarray): IPスキャンデータ
        s (list): キャリブレーションパラメータ
    
    Returns:
        tuple: (E, theta_rad, theta_deg) - エネルギー、ブラッグ角
    """
    # IP上の位置座標からブラッグ角への変換
    theta_rad = np.arctan((13 + s[0] - data[:,0])/(50 - s[1]))
    theta_deg = theta_rad * 180 / np.pi  # 度に変換
    
    # ブラッグの式を用いてエネルギーを計算 [keV]
    E = 12.3984 / (2 * D_HOPG * np.sin(theta_rad))
    
    return E, theta_rad, theta_deg

def calculate_energy_resolution(s):
    """エネルギー分解能計算関数の生成
    
    Args:
        s (list): キャリブレーションパラメータ
    
    Returns:
        function: dE/dx計算関数
    """
    def dEdx(x):
        a = 12.3984  # hc [keV·Å]
        b = 13 + s[0]  # 変換パラメータ
        c = 50 - s[1]  # 変換パラメータ
        theta = np.arctan((b-x)/c)
        
        # dθ/dx の計算
        def dtdx(x):
            return -c/(c**2 + (b - x)**2)
        
        # dE/dx = (∂E/∂θ) × (∂θ/∂x)
        return (a/(2*D_HOPG)) * (2 * np.cos(theta))/(-1 + np.cos(2*theta)) * dtdx(x)
    
    return dEdx

def convert_intensity_to_photon(data, E, t, dEdx_func, filter1, filter2):
    """IPの生データをフォトン数密度に変換
    
    Args:
        data (numpy.ndarray): IPスキャンデータ
        E (numpy.ndarray): エネルギー配列
        t (float): 時間遅延
        dEdx_func (function): エネルギー分解能関数
        filter1, filter2: フィルター透過率関数
    
    Returns:
        numpy.ndarray: フォトン数密度
    """
    # 各種補正を適用
    correction_factor = (ip_time_correction(t) * ip_energy_correction(E) * 
                        0.0025 * dEdx_func(data[:,0]) * filter1(E) * filter2(E))
    
    return data[:,1] * 1000 / correction_factor

def save_data(E, Photon, Energy, shot_num):
    """データの保存
    
    Args:
        E (numpy.ndarray): エネルギー配列
        Photon (numpy.ndarray): フォトン数密度
        Energy (float): レーザーエネルギー
        shot_num (str): ショット番号
    """
    # 生スペクトラムデータをCSVファイルに保存
    out_raw = np.array([E, Photon]).T
    np.savetxt('data/HOPG_{}.csv'.format(shot_num), out_raw, delimiter=',')
    
    # エネルギー規格化スペクトラムデータをCSVファイルに保存
    out_normalized = np.array([E, Photon/Energy]).T
    np.savetxt('data_normalize/HOPG_energy_normalize_{}.csv'.format(shot_num), 
               out_normalized, delimiter=',')

def plot_raw_spectrum(E, Photon, shot_num, Energy):
    """生スペクトラムのグラフ作成と保存
    
    Args:
        E (numpy.ndarray): エネルギー配列
        Photon (numpy.ndarray): フォトン数密度
        shot_num (str): ショット番号
        Energy (float): レーザーエネルギー
    """
    plt.rcParams.update({'font.size': 22})
    fig = plt.figure(figsize=(19.20, 10.80))
    ax = fig.add_subplot(111)
    
    # ショット番号とエネルギー情報をグラフに表示
    ax.text(np.min(E), 4.5E+5, "{}\n{} kJ".format(shot_num, Energy), 
            size=35, color='black')
    
    # 軸ラベルとフォーマット設定
    ax.set_xlabel("Energy keV", size=30)
    ax.set_ylabel("Photon/keV", size=30)
    ax.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.ylim([0, 6E+5])
    
    # スペクトラムをプロット
    ax.plot(E, Photon)
    
    # 生データのグラフを保存
    plt.savefig("data/HOPG_{}.png".format(shot_num), bbox_inches='tight', pad_inches=0)

def plot_normalized_spectrum(E, Photon, Energy, shot_num):
    """エネルギー規格化スペクトラムのグラフ作成と保存
    
    Args:
        E (numpy.ndarray): エネルギー配列
        Photon (numpy.ndarray): フォトン数密度
        Energy (float): レーザーエネルギー
        shot_num (str): ショット番号
    """
    plt.figure()
    plt.rcParams.update({'font.size': 22})
    fig = plt.figure(figsize=(19.20, 10.80))
    ax = fig.add_subplot(111)
    
    # ショット番号をグラフに表示
    ax.text(np.min(E), 4.5E+5, "{}".format(shot_num), size=35, color='black')
    
    # 軸ラベルとフォーマット設定
    ax.set_xlabel("Energy keV", size=30)
    ax.set_ylabel("Photon/keV/kJ", size=30)
    ax.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.ylim([0, 6E+5])
    
    # エネルギー規格化スペクトラムをプロット
    ax.plot(E, Photon/Energy)
    
    # エネルギー規格化グラフを保存
    plt.savefig("data_normalize/HOPG_energy_normalize_{}.png".format(shot_num), 
                bbox_inches='tight', pad_inches=0)

def extract_shot_number_from_filename(file_path):
    """ファイル名からショット番号を抽出する
    
    Args:
        file_path (str): データファイルのパス
    
    Returns:
        str: 抽出されたショット番号 (例: G43798, L12345)
        
    Raises:
        ValueError: ショット番号が見つからない場合
    """
    # ファイル名のみを取得（拡張子を除く）
    filename = os.path.splitext(os.path.basename(file_path))[0]
    
    # GXIIレーザー（Gから始まる）またはLFEXレーザー（Lから始まる）のショット番号を検索
    # パターン: _G##### または _L##### （#は数字）
    pattern = r'_([GL]\d+)_'
    match = re.search(pattern, filename)
    
    if match:
        return match.group(1)
    else:
        # パターンが見つからない場合の詳細エラー
        raise ValueError(
            "ファイル名からショット番号を抽出できませんでした。\n"
            "ファイル名: {}\n"
            "期待される形式: XXX_GYYYY_ または XXX_LYYYY_\n"
            "G: GXIIレーザー, L: LFEXレーザー".format(filename)
        )

def get_laser_type_from_shot_number(shot_num):
    """ショット番号からレーザータイプを判定する
    
    Args:
        shot_num (str): ショット番号
    
    Returns:
        str: レーザータイプ ("GXII" または "LFEX")
    """
    if shot_num.startswith('G'):
        return "GXII"
    elif shot_num.startswith('L'):
        return "LFEX"
    else:
        return "Unknown"

def get_file_path_gui():
    """GUIでファイルパスの選択を行う
    
    Returns:
        str: 選択されたデータファイルのパス（拡張子なし）、キャンセルの場合はNone
    """
    # Tkinterのルートウィンドウを作成（非表示）
    root = tk.Tk()
    root.withdraw()  # メインウィンドウを非表示にする
    
    # ファイル選択ダイアログを表示
    filetypes = [
        ('CSV files', '*.csv'),
        ('All files', '*.*')
    ]
    
    file_path = filedialog.askopenfilename(
        title='HOPG データファイルを選択してください',
        filetypes=filetypes,
        initialdir=os.getcwd()  # 現在のディレクトリから開始
    )
    
    # ルートウィンドウを破棄
    root.destroy()
    
    if file_path:
        # 拡張子を除いたパスを返す
        return os.path.splitext(file_path)[0]
    else:
        return None

def get_file_path_input():
    """ユーザーからファイルパスの入力を取得（CLI版）
    
    Returns:
        str: データファイルのパス（拡張子なし）
    """
    file_path = input("Input the data file path (without .csv extension): ")
    return file_path.strip()

def get_file_path_with_gui_fallback():
    """GUIでファイル選択、失敗時はCLI入力にフォールバック
    
    Returns:
        str: データファイルのパス（拡張子なし）
    """
    try:
        print("ファイル選択ダイアログを開いています...")
        file_path = get_file_path_gui()
        
        if file_path is None:
            print("ファイル選択がキャンセルされました。")
            print("手動でファイルパスを入力してください。")
            return get_file_path_input()
        else:
            print("選択されたファイル: {}.csv".format(file_path))
            return file_path
            
    except Exception as e:
        print("GUI ファイル選択でエラーが発生しました: {}".format(e))
        print("手動でファイルパスを入力してください。")
        return get_file_path_input()

def validate_data_file(file_path):
    """データファイルの存在確認
    
    Args:
        file_path (str): データファイルのパス（拡張子なし）
    
    Returns:
        bool: ファイルが存在する場合True
    """
    csv_file_path = file_path + '.csv'
    return os.path.exists(csv_file_path)

def main():
    """HOPG X線分光器データ解析のメイン処理"""
    
    print("=" * 60)
    print("HOPG X線分光器データ解析プログラム")
    print("=" * 60)
    
    try:
        # 1. ファイルパスの入力とショット番号の自動抽出
        print("1. データファイルの指定とショット番号の抽出")
        file_path = get_file_path_with_gui_fallback()
        
        # ファイルの存在確認
        if not validate_data_file(file_path):
            raise FileNotFoundError("データファイルが見つかりません: {}.csv".format(file_path))
        
        # ショット番号の自動抽出
        shot_num = extract_shot_number_from_filename(file_path)
        laser_type = get_laser_type_from_shot_number(shot_num)
        
        print("   - データファイル: {}.csv".format(file_path))
        print("   - ショット番号: {}".format(shot_num))
        print("   - レーザータイプ: {}".format(laser_type))
        
        # レーザーエネルギーの設定（必要に応じて調整）
        laser_energy = 1.406  # デフォルト値 [kJ]
        print("   - レーザーエネルギー: {} kJ (デフォルト値)".format(laser_energy))
        
        # 2. フィルターデータの読み込み
        print("\n2. フィルターデータの読み込み中...")
        filter1, filter2 = load_filter_data()
        print("   - ベリリウムフィルター読み込み完了")
        print("   - ポリエチレンフィルター読み込み完了")
        
        # 3. ユーザー入力の取得
        print("\n3. 測定パラメータの入力")
        time_delay = get_user_input(file_path, shot_num, laser_type)
        print("   - 使用する時間遅延: {:.3f} 時間".format(time_delay))
        
        # 4. 実験データの読み込み
        print("\n4. 実験データの読み込み中...")
        data = load_experimental_data(shot_num, file_path)
        print("   - データ点数: {} 点".format(len(data)))
        
        # 5. キャリブレーションパラメータの設定
        print("\n5. キャリブレーションパラメータの設定...")
        s, E_def, theta_def, place = calibrate_parameters()
        print("   - 基準エネルギー: {:.3f} keV".format(E_def))
        print("   - 基準位置: {:.3f} cm".format(place))
        print("   - キャリブレーション係数: s = [{:.6f}, {:.6f}]".format(s[0], s[1]))
        
        # 6. エネルギー変換
        print("\n6. エネルギー変換処理中...")
        E, theta_rad, theta_deg = convert_position_to_energy(data, s)
        print("   - エネルギー範囲: {:.2f} - {:.2f} keV".format(np.min(E), np.max(E)))
        
        # 7. エネルギー分解能関数の生成
        print("\n7. エネルギー分解能関数の生成...")
        dEdx_func = calculate_energy_resolution(s)
        
        # 8. フォトン数密度への変換
        print("\n8. フォトン数密度への変換処理中...")
        Photon = convert_intensity_to_photon(data, E, time_delay, dEdx_func, filter1, filter2)
        print("   - 最大フォトン数密度: {:.2e} photons/keV".format(np.max(Photon)))
        
        # 9. データの保存
        print("\n9. データの保存中...")
        save_data(E, Photon, laser_energy, shot_num)
        print("   - 生データ保存: data/HOPG_{}.csv".format(shot_num))
        print("   - 規格化データ保存: data_normalize/HOPG_energy_normalize_{}.csv".format(shot_num))
        
        # 10. グラフの作成と保存
        print("\n10. グラフの作成と保存中...")
        plot_raw_spectrum(E, Photon, shot_num, laser_energy)
        print("    - 生スペクトラム: data/HOPG_{}.png".format(shot_num))
        
        plot_normalized_spectrum(E, Photon, laser_energy, shot_num)
        print("    - 規格化スペクトラム: data_normalize/HOPG_energy_normalize_{}.png".format(shot_num))
        
        print("\n" + "=" * 60)
        print("解析完了！")
        print("ショット番号: {} ({} レーザー)".format(shot_num, laser_type))
        print("=" * 60)
        
    except FileNotFoundError as e:
        print("エラー: ファイルが見つかりません - {}".format(e))
        print("必要なファイル:")
        print("  - Be_500um_transmittance.dat")
        print("  - CH2_50um_transmittance.dat")
        print("  - 指定されたデータファイル (.csv)")
    except ValueError as e:
        print("エラー: {}".format(e))
        print("ファイル名の例:")
        print("  - profile/200911_G43801_HOPG_1212  (GXIIレーザー)")
        print("  - profile/210315_L12345_HOPG_0900  (LFEXレーザー)")
    except Exception as e:
        print("エラー: 解析中に問題が発生しました - {}".format(e))

if __name__ == "__main__":
    main()
