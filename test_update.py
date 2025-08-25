#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
キャリブレーション自動更新機能の実行テスト
"""

import numpy as np
import os
import datetime
import re

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
        
        new_s_line = 's = [{:.8f}, {:.8f}]  # 自動更新: {} (ショット: {})'.format(
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
                new_s_line, timestamp, shot_num
            )
            new_content = re.sub(pattern2, replacement, content)
        
        # パターン2でも見つからない場合、手動で特定の行を探す
        if new_content == content:
            lines = content.split('\n')
            for i, line in enumerate(lines):
                if 's = [' in line and ('キャリブレーション' in '\n'.join(lines[max(0, i-3):i+1]) or i == 96):  # 97行目
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

# テスト実行
if __name__ == "__main__":
    print("=== キャリブレーション自動更新機能テスト ===")
    
    # テスト用のキャリブレーションパラメータ
    test_s = np.array([-0.02123456, -0.81234567])
    test_shot_num = "G12345"
    test_file_path = "test/G12345_data"
    
    print("テスト用パラメータ:")
    print("  s = [{:.8f}, {:.8f}]".format(test_s[0], test_s[1]))
    print("  ショット番号: {}".format(test_shot_num))
    print()
    
    # 現在のパラメータを確認
    print("更新前のパラメータ確認:")
    os.system("grep -n 's = \\[' energy_conversion_HOPG.py | head -1")
    
    # 自動更新を実行
    print("\n自動更新実行中...")
    success = update_analysis_script_calibration(test_s, test_shot_num, test_file_path)
    
    if success:
        print("\n✓ 自動更新が成功しました")
        print("\n更新後のパラメータ確認:")
        os.system("grep -n 's = \\[' energy_conversion_HOPG.py | head -1")
    else:
        print("\n✗ 自動更新に失敗しました")
