# HOPG（Highly Ordered Pyrolytic Graphite）X線分光器のキャリブレーションと解析プログラム
# HOPGの回折を利用したX線エネルギー測定のデータ処理とキャリブレーション

import pprint
import numpy as np
import re
import os
import sys
import tkinter as tk
from tkinter import filedialog, messagebox
from time import sleep
from scipy.optimize import differential_evolution
from scipy import integrate
from numpy.linalg import solve
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import csv
import requests
from bs4 import BeautifulSoup
from datetime import datetime, timedelta
import urllib.parse

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

def get_user_input():
    """ユーザーからの入力取得
    
    Returns:
        float: 測定時間遅延 [時間]
    """
    print("測定時間遅延の入力方法を選択してください:")
    print("1. 自動計算（ファイル名とWebページから）")
    print("2. 手動入力")
    
    while True:
        choice = input("選択 (1 または 2): ").strip()
        if choice == "1":
            return "auto"
        elif choice == "2":
            t = input("Input the scanning time delay (hours): ")
            return float(t)
        else:
            print("1 または 2 を入力してください。")

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

def setup_calibration_references():
    """キャリブレーション用基準線の設定
    
    Returns:
        tuple: (E_def, theta_def, place) - 基準エネルギー、ブラッグ角、位置
    """
    # キャリブレーション用の基準線のエネルギー [keV]
    E_def = np.array([8.048, 8.391012])
    
    # 基準線のブラッグ角を計算 [rad]
    # ブラッグの式: nλ = 2d sinθ より θ = arcsin(12.3984/(2d*E))
    theta_def = np.arcsin(12.3984/(2.0 * D_HOPG * E_def))
    
    # 基準線のIP上での位置 [cm]
    place = np.array([1.005, 1.52])
    
    return E_def, theta_def, place

def calculate_calibration_parameters(E_def, theta_def, place):
    """キャリブレーションパラメータの計算
    
    Args:
        E_def (numpy.ndarray): 基準エネルギー配列
        theta_def (numpy.ndarray): 基準ブラッグ角配列
        place (numpy.ndarray): 基準位置配列
    
    Returns:
        numpy.ndarray: キャリブレーションパラメータ s
    """
    # 基準線の位置から座標変換パラメータを決定
    # 連立一次方程式を解いて、IP座標からブラッグ角への変換係数を求める
    left = [[1, np.tan(theta_def[0])],
            [1, np.tan(theta_def[1])]]
    
    right = [50 * np.tan(theta_def[0]) + place[0] - 13, 
             50 * np.tan(theta_def[1]) + place[1] - 13]
    
    # 連立方程式の解を求める
    s = solve(left, right)
    
    return s

def convert_position_to_energy(data, s):
    """IP上の位置座標からエネルギーへの変換
    
    Args:
        data (numpy.ndarray): IPスキャンデータ
        s (numpy.ndarray): キャリブレーションパラメータ
    
    Returns:
        tuple: (E, theta_rad, theta_deg) - エネルギー、ブラッグ角
    """
    # IP上の位置座標からブラッグ角への変換
    theta_rad = np.arctan((13 + s[0] - data[:,0])/(50 - s[1]))
    theta_deg = theta_rad * 180 / np.pi  # 度に変換
    
    # ブラッグの式を用いてエネルギーを計算 [keV]
    # E = hc/λ = 12.3984 / (2d sinθ)
    E = 12.3984 / (2 * D_HOPG * np.sin(theta_rad))
    
    return E, theta_rad, theta_deg

def calculate_energy_resolution(s):
    """エネルギー分解能計算関数の生成
    
    Args:
        s (numpy.ndarray): キャリブレーションパラメータ
    
    Returns:
        function: dE/dx計算関数
    """
    def dEdx(x):
        a = 12.3984  # hc [keV·Å]
        b = 13 + s[0]  # 変換パラメータ
        c = 50 - s[1]   # 変換パラメータ
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

def update_analysis_script_calibration(s, shot_num, file_path):
    """energy_conversion_HOPG.pyのキャリブレーションパラメータを直接更新
    
    Args:
        s (numpy.ndarray): キャリブレーションパラメータ
        shot_num (str): ショット番号
        file_path (str): データファイルパス
    """
    import datetime
    import re
    
    # 現在の日時を取得
    now = datetime.datetime.now()
    timestamp = now.strftime("%Y年%m月%d日 %H時%M分%S秒")
    
    script_path = 'energy_conversion_HOPG.py'
    
    try:
        # energy_conversion_HOPG.pyファイルを読み込み
        if not os.path.exists(script_path):
            print("    - energy_conversion_HOPG.py が見つかりません")
            return False
        
        with open(script_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # バックアップを作成
        backup_path = 'energy_conversion_HOPG_backup_{}.py'.format(
            now.strftime("%Y%m%d_%H%M%S")
        )
        with open(backup_path, 'w', encoding='utf-8') as f:
            f.write(content)
        
        # キャリブレーションパラメータの行を探して置換
        # パターン1: calibrate_parameters関数内のsの定義を探す
        pattern1 = r'(\s+# 事前に決定されたキャリブレーションパラメータ.*?\n)(\s+s = \[.*?\])(.*?\n)'
        
        new_s_line = '    s = [{:.8f}, {:.8f}]  # 自動更新: {} (ショット: {})'.format(
            s[0], s[1], timestamp, shot_num
        )
        
        def replace_calibration_params(match):
            indent = '    '  # calibrate_parameters関数内のインデント
            comment_line = match.group(1)
            old_line = match.group(2)
            rest = match.group(3)
            return comment_line + indent + new_s_line + rest
        
        new_content = re.sub(pattern1, replace_calibration_params, content, flags=re.DOTALL)
        
        # パターン1で見つからない場合、より一般的なパターンで探す
        if new_content == content:
            pattern2 = r'(\s+s = \[.*?\].*?\n)'
            replacement = '    {}  # 自動更新: {} (ショット: {})\n'.format(
                new_s_line.strip(), timestamp, shot_num
            )
            new_content = re.sub(pattern2, replacement, content)
        
        # パターン2でも見つからない場合、手動で特定の行を探す
        if new_content == content:
            lines = content.split('\n')
            for i, line in enumerate(lines):
                if 's = [' in line and 'キャリブレーション' in lines[max(0, i-3):i+1]:
                    # 見つかった行を置換
                    indent = '    '
                    lines[i] = indent + new_s_line
                    new_content = '\n'.join(lines)
                    break
        
        # ファイルの更新
        if new_content != content:
            with open(script_path, 'w', encoding='utf-8') as f:
                f.write(new_content)
            
            print("    - energy_conversion_HOPG.py を自動更新しました")
            print("    - バックアップ作成: {}".format(backup_path))
            print("    - 新しいパラメータ: s = [{:.8f}, {:.8f}]".format(s[0], s[1]))
            print("    - 更新日時: {}".format(timestamp))
            print("    - キャリブレーション元: ショット {}".format(shot_num))
            return True
        else:
            print("    - キャリブレーションパラメータの更新対象が見つかりませんでした")
            # バックアップファイルを削除（変更がないため）
            if os.path.exists(backup_path):
                os.remove(backup_path)
            return False
            
    except Exception as e:
        print("    - 自動更新中にエラーが発生: {}".format(e))
        return False

def extract_reading_time_from_filename(file_path):
    """ファイル名からIP読み取り時刻を抽出する
    
    Args:
        file_path (str): データファイルのパス
    
    Returns:
        tuple: (hours, minutes) - 読み取り時刻の時間と分
        
    Raises:
        ValueError: 読み取り時刻が見つからない場合
    """
    # ファイル名のみを取得
    filename = os.path.basename(file_path)
    
    # _HOPG_XXXX パターンを検索（XXXXは4桁の数字）
    pattern = r'_HOPG_(\d{4})'
    match = re.search(pattern, filename)
    
    if match:
        time_str = match.group(1)
        hours = int(time_str[:2])    # 最初の2桁が時間
        minutes = int(time_str[2:])  # 残りの2桁が分
        
        # 時間と分の妥当性チェック
        if 0 <= hours <= 23 and 0 <= minutes <= 59:
            return hours, minutes
        else:
            raise ValueError("無効な時刻: {}時{}分".format(hours, minutes))
    else:
        raise ValueError("ファイル名にIP読み取り時刻（_HOPG_XXXX）が見つかりません")

def fetch_shot_time_from_web(shot_num, laser_type):
    """Webページからショット時刻を取得する
    
    Args:
        shot_num (str): ショット番号（GまたはLプレフィックス付き）
        laser_type (str): レーザータイプ（"GXII" または "LFEX"）
    
    Returns:
        tuple: (hours, minutes) - ショット時刻の時間と分、取得失敗時は (None, None)
    """
    try:
        # ショット番号からGまたはLプレフィックスを除去
        shot_id = shot_num[1:] if shot_num.startswith(('G', 'L')) else shot_num
        
        # URLの構築
        if laser_type == "GXII":
            url = "http://god.ile.osaka-u.ac.jp/gxii/ShotDataViewer/ShotDataViewer.php?Id={}".format(shot_id)
        elif laser_type == "LFEX":
            url = "http://god.ile.osaka-u.ac.jp/lfex-fe/ShotDataViewer/ShotDataViewer?Id={}".format(shot_id)
        else:
            print("    - ✗ エラー: 不明なレーザータイプ: {}".format(laser_type))
            print("    - 対応しているレーザータイプ: GXII, LFEX")
            print("    - プログラムを終了します")
            sys.exit(1)
        
        print("    - ショット時刻取得中: {}".format(url))
        
        # Webページを取得
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
        }
        response = requests.get(url, headers=headers, timeout=10)
        response.raise_for_status()
        
        # HTMLを解析
        soup = BeautifulSoup(response.content, 'html.parser')
        
        # HTMLから直接大文字のMAINを検索
        html_content = response.text
        main_pos = html_content.find('MAIN')
        
        if main_pos != -1:
            print("    - 大文字MAIN発見（HTML直接検索）")
            
            # MAINより後の部分で最初の時刻を探す
            after_main = html_content[main_pos:]
            time_pattern = r'(\d{1,2}):(\d{2})'
            time_matches = re.findall(time_pattern, after_main)
            
            if time_matches:
                print("    - MAIN後の時刻候補: {}".format(time_matches))
                
                # 最初の有効な時刻を使用（MAINの直後）
                for time_match in time_matches:
                    try:
                        hours = int(time_match[0])
                        minutes = int(time_match[1])
                        
                        if 0 <= hours <= 23 and 0 <= minutes <= 59:
                            print("    - ✓ ショット時刻（MAIN直後）: {}時{}分".format(hours, minutes))
                            return hours, minutes
                    except ValueError:
                        continue
            
            # MAIN後に時刻が見つからない場合
            print("    - ✗ エラー: MAIN後に有効な時刻が見つかりませんでした")
            print("    - 手動入力が必要です")
            return "TIME_NOT_FOUND", None
        
        # MAINが見つからない場合
        print("    - ✗ エラー: MAINキーワードが見つかりませんでした")
        print("    - 手動入力が必要です")
        return "MAIN_NOT_FOUND", None
            
    except requests.exceptions.RequestException as e:
        print("    - ✗ Webページ取得エラー: {}".format(e))
        print("    - 手動入力が必要です")
        return "WEB_ERROR", None
    except Exception as e:
        print("    - ✗ ショット時刻取得エラー: {}".format(e))
        print("    - 手動入力が必要です")
        return "UNKNOWN_ERROR", None

def calculate_time_delay_auto(file_path, shot_num, laser_type):
    """ファイル名とWebページから自動的に時間遅延を計算する
    
    Args:
        file_path (str): データファイルのパス
        shot_num (str): ショット番号
        laser_type (str): レーザータイプ
    
    Returns:
        float: 時間遅延 [時間]、計算失敗時は None
    """
    try:
        # 1. ファイル名からIP読み取り時刻を抽出
        print("    - ファイル名からIP読み取り時刻を抽出中...")
        read_hours, read_minutes = extract_reading_time_from_filename(file_path)
        print("    - IP読み取り時刻: {}時{}分".format(read_hours, read_minutes))
        
        # 2. Webページからショット時刻を取得
        print("    - Webページからショット時刻を取得中...")
        shot_hours, shot_minutes = fetch_shot_time_from_web(shot_num, laser_type)
        
        # エラー時の手動入力フォールバック
        if isinstance(shot_hours, str):  # エラーコードが返された場合
            print("    - ✗ 自動取得に失敗しました")
            print("    - 手動でショット時刻を入力してください")
            
            while True:
                try:
                    shot_time_input = input("    ショット時刻を入力 (HH:MM形式): ").strip()
                    if re.match(r'^\d{1,2}:\d{2}$', shot_time_input):
                        time_parts = shot_time_input.split(':')
                        shot_hours = int(time_parts[0])
                        shot_minutes = int(time_parts[1])
                        
                        if 0 <= shot_hours <= 23 and 0 <= shot_minutes <= 59:
                            print("    - 手動入力されたショット時刻: {}時{}分".format(shot_hours, shot_minutes))
                            break
                        else:
                            print("    - 無効な時刻です。0-23時、0-59分の範囲で入力してください")
                    else:
                        print("    - HH:MM形式で入力してください（例：14:20）")
                except (ValueError, KeyboardInterrupt):
                    print("    - 入力が無効です。再度入力してください")
        
        if shot_hours is None or shot_minutes is None:
            print("    - ✗ ショット時刻が取得できませんでした")
            return None
        
        # 3. 時間差分を計算
        # 読み取り時刻とショット時刻を分単位で計算
        read_total_minutes = read_hours * 60 + read_minutes
        shot_total_minutes = shot_hours * 60 + shot_minutes
        
        # 時間差分（読み取り時刻 - ショット時刻）
        time_diff_minutes = read_total_minutes - shot_total_minutes
        
        # 日をまたいだ場合の処理（負の値の場合は翌日とみなす）
        if time_diff_minutes < 0:
            time_diff_minutes += 24 * 60  # 24時間を加算
        
        # 時間単位に変換
        time_delay_hours = time_diff_minutes / 60.0
        
        print("    - 時間差分計算結果:")
        print("      ショット時刻: {}時{}分".format(shot_hours, shot_minutes))
        print("      IP読み取り時刻: {}時{}分".format(read_hours, read_minutes))
        print("      時間遅延: {:.2f} 時間".format(time_delay_hours))
        
        return time_delay_hours
        
    except ValueError as e:
        print("    - 時間抽出エラー: {}".format(e))
        return None
    except Exception as e:
        print("    - 時間遅延計算エラー: {}".format(e))
        return None

def main():
    """HOPG X線分光器キャリブレーションと解析のメイン処理"""
    
    print("=" * 70)
    print("HOPG X線分光器キャリブレーションと解析プログラム")
    print("=" * 70)
    
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
        laser_energy = 1.325  # デフォルト値 [kJ]
        print("   - レーザーエネルギー: {} kJ (デフォルト値)".format(laser_energy))
        
        # 2. フィルターデータの読み込み
        print("\n2. フィルターデータの読み込み中...")
        filter1, filter2 = load_filter_data()
        print("   - ベリリウムフィルター読み込み完了")
        print("   - ポリエチレンフィルター読み込み完了")
        
        # 3. 測定パラメータの入力
        print("\n3. 測定パラメータの入力")
        time_input = get_user_input()
        
        if time_input == "auto":
            # 自動計算モード
            print("   - 自動時間遅延計算を開始します...")
            time_delay = calculate_time_delay_auto(file_path, shot_num, laser_type)
            
            if time_delay is None:
                print("   - 自動計算に失敗しました。手動入力に切り替えます。")
                time_delay = float(input("   手動で時間遅延を入力してください（時間）: "))
        else:
            # 手動入力モード
            time_delay = time_input
        
        print("   - 使用する時間遅延: {:.3f} 時間".format(time_delay))
        
        # 4. 実験データの読み込み
        print("\n4. 実験データの読み込み中...")
        data = load_experimental_data(shot_num, file_path)
        print("   - データ点数: {} 点".format(len(data)))
        
        # 5. キャリブレーション基準線の設定
        print("\n5. キャリブレーション基準線の設定...")
        E_def, theta_def, place = setup_calibration_references()
        print("   - 基準エネルギー 1: {:.3f} keV".format(E_def[0]))
        print("   - 基準エネルギー 2: {:.3f} keV".format(E_def[1]))
        print("   - 基準位置 1: {:.3f} cm".format(place[0]))
        print("   - 基準位置 2: {:.3f} cm".format(place[1]))
        
        # 6. キャリブレーションパラメータの計算
        print("\n6. キャリブレーションパラメータの計算...")
        s = calculate_calibration_parameters(E_def, theta_def, place)
        print("   - キャリブレーション係数: s = [{:.6f}, {:.6f}]".format(s[0], s[1]))
        
        # 7. エネルギー変換
        print("\n7. エネルギー変換処理中...")
        E, theta_rad, theta_deg = convert_position_to_energy(data, s)
        print("   - エネルギー範囲: {:.2f} - {:.2f} keV".format(np.min(E), np.max(E)))
        print("   - ブラッグ角範囲: {:.2f} - {:.2f} 度".format(np.min(theta_deg), np.max(theta_deg)))
        
        # 8. エネルギー分解能関数の生成
        print("\n8. エネルギー分解能関数の生成...")
        dEdx_func = calculate_energy_resolution(s)
        
        # 9. フォトン数密度への変換
        print("\n9. フォトン数密度への変換処理中...")
        Photon = convert_intensity_to_photon(data, E, time_delay, dEdx_func, filter1, filter2)
        print("   - 最大フォトン数密度: {:.2e} photons/keV".format(np.max(Photon)))
        print("   - 積分フォトン数: {:.2e} photons".format(np.trapz(Photon, E)))
        
        # 10. データの保存
        print("\n10. データの保存中...")
        save_data(E, Photon, laser_energy, shot_num)
        print("    - 生データ保存: data/HOPG_{}.csv".format(shot_num))
        print("    - 規格化データ保存: data_normalize/HOPG_energy_normalize_{}.csv".format(shot_num))
        
        # 11. グラフの作成と保存
        print("\n11. グラフの作成と保存中...")
        plot_raw_spectrum(E, Photon, shot_num, laser_energy)
        print("    - 生スペクトラム: data/HOPG_{}.png".format(shot_num))
        
        plot_normalized_spectrum(E, Photon, laser_energy, shot_num)
        print("    - 規格化スペクトラム: data_normalize/HOPG_energy_normalize_{}.png".format(shot_num))
        
        # 12. キャリブレーションパラメータの表示
        print("\n12. キャリブレーション結果")
        print("    - 計算されたキャリブレーションパラメータ:")
        pprint.pprint(s)
        print("    - これらの値を energy_conversion_HOPG.py で使用してください")
        
        # 13. energy_conversion_HOPG.py のキャリブレーションパラメータを自動更新
        print("\n13. energy_conversion_HOPG.py のキャリブレーションパラメータを自動更新中...")
        if update_analysis_script_calibration(s, shot_num, file_path):
            print("    - energy_conversion_HOPG.py のキャリブレーションパラメータを更新しました")
        else:
            print("    - energy_conversion_HOPG.py のキャリブレーションパラメータの更新に失敗しました")
        
        print("\n" + "=" * 70)
        print("キャリブレーションと解析完了！")
        print("ショット番号: {} ({} レーザー)".format(shot_num, laser_type))
        print("=" * 70)
        
    except FileNotFoundError as e:
        print("エラー: ファイルが見つかりません - {}".format(e))
        print("必要なファイル:")
        print("  - Be_500um_transmittance.dat")
        print("  - CH2_50um_transmittance.dat")
        print("  - 指定されたデータファイル (.csv)")
    except ValueError as e:
        print("エラー: {}".format(e))
        print("ファイル名の例:")
        print("  - profile/200910_G43798_HOPG_1505  (GXIIレーザー)")
        print("  - profile/210315_L12345_HOPG_0900  (LFEXレーザー)")
    except Exception as e:
        print("エラー: 解析中に問題が発生しました - {}".format(e))

if __name__ == "__main__":
    main()