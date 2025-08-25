# HOPG（Highly Ordered Pyrolytic Graphite）X線分光器データ解析プログラム
# 事前キャリブレーション済みパラメータを使用してスペクトラム解析を実行

import numpy as np
import re
import os
import sys
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

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
        t (float): 測定時間遅延 [分]
    
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
    filter1_data = np.loadtxt('filter/Be_500um_transmittance.dat', skiprows=2, unpack=True)
    filter1 = interp1d(filter1_data[0,:]*1.0E-03, filter1_data[1,:], kind='cubic')
    
    # CH2（ポリエチレン）フィルター 50μm厚の透過率データ  
    filter2_data = np.loadtxt('filter/CH2_50um_transmittance.dat', skiprows=2, unpack=True)
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
                    t = input("   手動で時間遅延を入力してください（分）: ")
                    return float(t)
                else:
                    print("   - 自動計算完了: {:.1f} 分 ({:.3f} 時間)".format(time_delay, time_delay/60.0))
                    return time_delay
            elif choice == "2":
                t = input("Input the scanning time delay (minutes): ")
                return float(t)
            else:
                print("1 または 2 を入力してください。")
    else:
        # 自動計算が利用できない場合は手動入力のみ
        if not AUTO_TIME_AVAILABLE:
            print("注意: 自動時間計算機能が利用できません（キャリブレーションスクリプトが必要）")
        t = input("Input the scanning time delay (minutes): ")
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

def interactive_single_calibration_selection(data, E_def):
    """対話的な単一キャリブレーション基準線位置選択
    
    Args:
        data (numpy.ndarray): IPスキャンデータ
        E_def (float): 基準エネルギー
    
    Returns:
        float: 選択された基準位置
    """
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Button
    
    print("   - 対話的な基準線位置選択を開始します...")
    print("   - グラフ上で1つの基準線位置をクリックして選択してください")
    print("   - 選択後、縦線をドラッグして微調整できます")
    print("   - 基準エネルギー: {:.3f} keV".format(E_def))
    
    # グラフの設定
    plt.rcParams.update({'font.size': 12})
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # 生データをプロット
    ax.plot(data[:,0], data[:,1], 'b-', linewidth=2, label='IP Scan Data')
    ax.set_xlabel('Position [cm]', fontsize=14)
    ax.set_ylabel('Intensity [PSL]', fontsize=14)
    ax.set_title('Select Single Calibration Reference Position\nClick to place, then drag to fine-tune position', fontsize=16)
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # 選択された位置を保存する変数
    selected_position = [None]  # リストで包んで参照渡しにする
    vertical_line = [None]
    text_annotation = [None]
    
    # ドラッグ可能な縦線クラス（単一版）
    class DraggableSingleLine:
        def __init__(self, line, text, ax, selected_position):
            self.line = line
            self.text = text
            self.ax = ax
            self.selected_position = selected_position
            self.press = None
            self.is_dragging = False
            
        def connect(self):
            """イベントを接続"""
            self.cidpress = self.line.figure.canvas.mpl_connect('button_press_event', self.on_press)
            self.cidrelease = self.line.figure.canvas.mpl_connect('button_release_event', self.on_release)
            self.cidmotion = self.line.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
            
        def on_press(self, event):
            """マウス押下時の処理"""
            if event.inaxes != self.ax:
                return
            
            # 縦線の位置を取得
            line_x = self.line.get_xdata()[0]
            
            # クリック位置が縦線に近いかチェック（許容範囲: 0.05cm）
            tolerance = 0.05
            if abs(event.xdata - line_x) <= tolerance:
                self.press = (event.xdata, event.ydata)
                self.is_dragging = True
                # カーソルを変更
                self.ax.figure.canvas.set_cursor(1)  # 手のカーソル
                
        def on_motion(self, event):
            """マウス移動時の処理"""
            if not self.is_dragging or self.press is None:
                return
            
            if event.inaxes != self.ax:
                return
                
            # 新しいx位置を取得
            new_x = event.xdata
            if new_x is None:
                return
                
            # 線の位置を更新
            self.line.set_xdata([new_x, new_x])
            
            # テキストの位置を更新
            y_max = np.max(data[:,1])
            self.text.set_position((new_x, y_max * 0.9))
            self.text.set_text(f'{E_def:.3f} keV\n{new_x:.3f} cm')
            
            # 選択位置を更新
            self.selected_position[0] = new_x
            
            # 画面を更新
            self.ax.figure.canvas.draw_idle()
            
        def on_release(self, event):
            """マウス離上時の処理"""
            if self.is_dragging:
                self.is_dragging = False
                self.press = None
                # カーソルを元に戻す
                self.ax.figure.canvas.set_cursor(0)  # 通常のカーソル
                
                if event.inaxes == self.ax and event.xdata is not None:
                    print(f"   - 基準線を {event.xdata:.3f} cm に調整しました ({E_def:.3f} keV)")
    
    # ドラッグ可能オブジェクト
    draggable_object = [None]
    
    def on_click(event):
        """クリックイベントハンドラ"""
        if event.inaxes != ax:
            return
        
        # 既存のドラッグ処理中は新しい線を追加しない
        if draggable_object[0] is not None and draggable_object[0].is_dragging:
            return
            
        if selected_position[0] is None:
            x_pos = event.xdata
            selected_position[0] = x_pos
            
            # 既存の線があれば削除
            if vertical_line[0] is not None:
                vertical_line[0].remove()
            if text_annotation[0] is not None:
                text_annotation[0].remove()
            
            # 縦線を描画（適度な太さでドラッグしやすく）
            line = ax.axvline(x=x_pos, color='red', linestyle='--', linewidth=2, 
                             alpha=0.8, picker=True, pickradius=15,
                             label=f'Reference: {E_def:.3f} keV')
            vertical_line[0] = line
            
            # テキスト表示
            y_max = np.max(data[:,1])
            text = ax.text(x_pos, y_max * 0.9, f'{E_def:.3f} keV\n{x_pos:.3f} cm', 
                          ha='center', va='top', fontsize=11, 
                          bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.8))
            text_annotation[0] = text
            
            # ドラッグ可能オブジェクトを作成
            draggable = DraggableSingleLine(line, text, ax, selected_position)
            draggable.connect()
            draggable_object[0] = draggable
            
            print(f"   - 基準線: {x_pos:.3f} cm ({E_def:.3f} keV)")
            print(f"   - 赤い線をドラッグして位置を微調整できます")
            print("   - 位置が決まったら 'Confirm Selection' をクリックしてください")
            
            confirm_button.label.set_text('Confirm Selection')
            plt.draw()
    
    def confirm_selection(event):
        """選択確定ボタンのイベントハンドラ"""
        if selected_position[0] is not None:
            plt.close(fig)
        else:
            print("   - 1つの位置を選択してから確定してください")
    
    # ボタンの追加（Confirmボタンのみ）
    ax_confirm = plt.axes([0.8, 0.02, 0.15, 0.04])
    confirm_button = Button(ax_confirm, 'Select Position')
    confirm_button.on_clicked(confirm_selection)
    
    # クリックイベントを接続
    fig.canvas.mpl_connect('button_press_event', on_click)
    
    # 使用方法の表示
    instruction_text = ("Instructions:\n"
                       "1. Click on the graph to place reference line\n"
                       "2. Drag the red line to fine-tune position\n"
                       "3. Click 'Confirm Selection' when done")
    ax.text(0.02, 0.98, instruction_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=dict(boxstyle='round,pad=0.5', facecolor='lightblue', alpha=0.8))
    
    # グラフを表示
    plt.tight_layout()
    plt.show()
    
    if selected_position[0] is not None:
        place = selected_position[0]
        print("   - 選択された基準位置: {:.3f} cm ({:.3f} keV)".format(place, E_def))
        return place
    else:
        print("   - 選択がキャンセルされました。デフォルト値を使用します")
        return 1.005  # デフォルト値

def calibrate_parameters(data=None, interactive_mode=False):
    """キャリブレーションパラメータの設定と微調整
    
    Args:
        data (numpy.ndarray, optional): IPスキャンデータ（インタラクティブモード用）
        interactive_mode (bool): インタラクティブキャリブレーションを使用するかどうか
    
    Returns:
        tuple: (s, E_def, theta_def, place) - キャリブレーションパラメータ
    """
    # 基準線のエネルギー [keV]
    E_def = 8.048
    # 基準線のブラッグ角 [rad]
    theta_def = np.arcsin(12.3984/(2.0 * D_HOPG * E_def))
    
    # 基準線のIP上での位置 [cm]
    if interactive_mode and data is not None:
        try:
            print("   - インタラクティブモードで基準線位置を選択します")
            place = interactive_single_calibration_selection(data, E_def)
        except Exception as e:
            print(f"   - インタラクティブ選択でエラーが発生しました: {e}")
            print("   - 最大値の位置を自動選択します")
            max_index = np.argmax(data[:,1])
            place = data[max_index, 0]
            print(f"   - 最大値位置: {place:.3f} cm (強度: {data[max_index, 1]:.2f})")
            place = data[max_index, 0]
            print(f"   - 最大値位置: {place:.3f} cm (強度: {data[max_index, 1]:.2f})")
    else:
        place = 1.005  # 手動で設定された基準線位置
    
    # 事前に決定されたキャリブレーションパラメータ
    s = [-0.00905962, -0.84425207]  # 自動更新: 2025年08月25日 17時38分55秒 (ショット: G43798)
    
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
    """IPの生データを絶対フォトン数密度に変換
    
    Args:
        data (numpy.ndarray): IPスキャンデータ
        E (numpy.ndarray): エネルギー配列
        t (float): 時間遅延 [分]
        dEdx_func (function): エネルギー分解能関数
        filter1, filter2: フィルター透過率関数
    
    Returns:
        numpy.ndarray: 絶対フォトン数密度
    """
    # 結晶反射率
    crystal_reflectivity = 5E-08
    
    # 各種補正を適用（結晶反射率を含む）
    correction_factor = (ip_time_correction(t) * ip_energy_correction(E) * 
                        0.0025 * dEdx_func(data[:,0]) * filter1(E) * filter2(E) * 
                        crystal_reflectivity)
    
    return data[:,1] * 1000 / correction_factor

def ensure_directory_exists(directory):
    """ディレクトリが存在しない場合は作成する
    
    Args:
        directory (str): 作成するディレクトリのパス
    """
    if not os.path.exists(directory):
        try:
            os.makedirs(directory)
            print(f"   - ディレクトリを作成しました: {directory}")
        except OSError as e:
            print(f"   - ディレクトリ作成エラー: {e}")
            raise

def save_data(E, Photon, shot_num):
    """データの保存
    
    Args:
        E (numpy.ndarray): エネルギー配列
        Photon (numpy.ndarray): 絶対フォトン数密度
        shot_num (str): ショット番号
    """
    # データフォルダの存在確認・作成
    ensure_directory_exists('data')
    
    # 絶対フォトン数スペクトラムデータをCSVファイルに保存
    out_raw = np.array([E, Photon]).T
    np.savetxt('data/HOPG_{}.csv'.format(shot_num), out_raw, delimiter=',')
    
    print("   - 絶対フォトン数スペクトラムを保存しました: data/HOPG_{}.csv".format(shot_num))
def plot_spectrum(E, Photon, shot_num):
    """絶対フォトン数スペクトラムのグラフ作成と保存
    
    Args:
        E (numpy.ndarray): エネルギー配列
        Photon (numpy.ndarray): 絶対フォトン数密度
        shot_num (str): ショット番号
    """
    # データフォルダの存在確認・作成
    ensure_directory_exists('data')
    
    # フォントとグラフの設定
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
        'font.size': 14,
        'axes.linewidth': 1.5,
        'axes.labelsize': 16,
        'axes.titlesize': 18,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'legend.fontsize': 12,
        'figure.dpi': 300
    })
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # スペクトラムをプロット（見やすいスタイル）
    ax.plot(E, Photon, linewidth=2, color='blue', alpha=0.8)
    
    # 軸ラベルとフォーマット設定
    ax.set_xlabel("Energy [keV]", fontweight='bold')
    ax.set_ylabel("Absolute Photon Number [Photon/str/keV]", fontweight='bold')
    ax.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    
    # グリッドを追加して見やすくする
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    
    # 軸の範囲を自動調整
    ax.set_xlim(np.min(E), np.max(E))
    ax.set_ylim(0, np.max(Photon) * 1.1)
    
    # レイアウトを調整
    plt.tight_layout()
    
    # PDFとして保存
    plt.savefig("data/HOPG_{}.pdf".format(shot_num), bbox_inches='tight', pad_inches=0.1)
    print("   - スペクトラムグラフを保存しました: data/HOPG_{}.pdf".format(shot_num))
    
    # メモリを解放
    plt.close(fig)

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
    """HOPG X線分光器データ解析のメイン処理
    
    コマンドライン使用例:
    python energy_conversion_HOPG.py                          # 通常モード（自動時間計算）
    python energy_conversion_HOPG.py -c                       # インタラクティブキャリブレーションモード
    python energy_conversion_HOPG.py --calibration-mode       # インタラクティブキャリブレーションモード
    python energy_conversion_HOPG.py --time-delay 120         # 手動時間指定（120分）
    python energy_conversion_HOPG.py -c --time-delay 90       # インタラクティブ + 手動時間指定
    """
    
    print("=" * 60)
    print("HOPG X線分光器データ解析プログラム")
    print("=" * 60)
    
    # コマンドライン引数の処理
    import argparse
    parser = argparse.ArgumentParser(description='HOPG X線分光器データ解析プログラム')
    parser.add_argument('-t', '--time-delay', type=float,
                       help='時間遅延を手動指定 [分] (指定しない場合は自動計算)')
    parser.add_argument('-c', '--calibration-mode', action='store_true',
                       help='インタラクティブキャリブレーションモードを有効にする')
    
    args = parser.parse_args()

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
        
        # 2. フィルターデータの読み込み
        print("\n2. フィルターデータの読み込み中...")
        filter1, filter2 = load_filter_data()
        print("   - ベリリウムフィルター読み込み完了")
        print("   - ポリエチレンフィルター読み込み完了")
        
        # 3. ユーザー入力の取得
        print("\n3. 測定パラメータの入力")
        
        # 時間遅延の取得（簡素化されたロジック）
        if args.time_delay is not None:
            # 手動指定された場合
            time_delay = args.time_delay
            print(f"   - 指定された時間遅延: {time_delay:.1f} 分 ({time_delay/60.0:.3f} 時間)")
        else:
            # 自動計算を試行
            print("   - 自動時間遅延計算を開始します...")
            time_delay = calculate_time_delay_auto(file_path, shot_num, laser_type)
            
            if time_delay is None:
                print("   - 自動計算に失敗しました。手動入力に切り替えます。")
                time_delay = float(input("   手動で時間遅延を入力してください（分）: "))
        
        print("   - 使用する時間遅延: {:.1f} 分 ({:.3f} 時間)".format(time_delay, time_delay/60.0))
        
        # 4. 実験データの読み込み
        print("\n4. 実験データの読み込み中...")
        data = load_experimental_data(shot_num, file_path)
        print("   - データ点数: {} 点".format(len(data)))
        
        # 5. キャリブレーションパラメータの設定
        print("\n5. キャリブレーションパラメータの設定...")
        if args.calibration_mode:
            print("   - インタラクティブキャリブレーションモードが有効です")
            s, E_def, theta_def, place = calibrate_parameters(data, interactive_mode=True)
        else:
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
        
        # 8. 絶対フォトン数密度への変換
        print("\n8. 絶対フォトン数密度への変換処理中...")
        Photon = convert_intensity_to_photon(data, E, time_delay, dEdx_func, filter1, filter2)
        print("   - 結晶反射率: 5E-08 適用済み")
        print("   - 最大絶対フォトン数密度: {:.2e} photons/keV".format(np.max(Photon)))
        
        # 9. データの保存とグラフ作成
        print("\n9. データの保存とグラフ作成中...")
        save_data(E, Photon, shot_num)
        plot_spectrum(E, Photon, shot_num)
        
        print("\n" + "=" * 60)
        print("解析完了！")
        print("ショット番号: {} ({} レーザー)".format(shot_num, laser_type))
        if args.calibration_mode:
            print("使用モード: インタラクティブキャリブレーション")
        else:
            print("使用モード: 自動キャリブレーション")
        print("出力ファイル:")
        print("  - 絶対フォトン数スペクトラム: data/HOPG_{}.csv".format(shot_num))
        print("  - スペクトラムグラフ: data/HOPG_{}.pdf".format(shot_num))
        if args.calibration_mode:
            print("ヒント: 次回は同じ位置を使う場合、通常モードでも解析できます")
        print("=" * 60)
        
    except FileNotFoundError as e:
        print("エラー: ファイルが見つかりません - {}".format(e))
        print("必要なファイル:")
        print("  - filter/Be_500um_transmittance.dat")
        print("  - filter/CH2_50um_transmittance.dat")
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
