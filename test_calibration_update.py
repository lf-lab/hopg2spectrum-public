#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
テスト用キャリブレーション自動更新機能の検証プログラム
"""

import numpy as np
import os
import datetime

# energy_conversion_HOPG_calibration.pyから自動更新関数をインポート
import sys
sys.path.append('.')

def test_calibration_update():
    """キャリブレーション自動更新機能のテスト"""
    
    # テスト用のキャリブレーションパラメータ
    test_s = np.array([-0.02123456, -0.81234567])
    test_shot_num = "G12345"
    test_file_path = "test/G12345_data"
    
    print("=== キャリブレーション自動更新機能テスト ===")
    print("テスト用パラメータ:")
    print(f"  s = [{test_s[0]:.8f}, {test_s[1]:.8f}]")
    print(f"  ショット番号: {test_shot_num}")
    print()
    
    # energy_conversion_HOPG_calibration.pyから関数をインポート
    try:
        from energy_conversion_HOPG_calibration import update_analysis_script_calibration
        
        # 自動更新を実行
        print("自動更新実行中...")
        success = update_analysis_script_calibration(test_s, test_shot_num, test_file_path)
        
        if success:
            print("✓ 自動更新が成功しました")
            
            # 更新後のファイルを確認
            print("\n更新後のenergy_conversion_HOPG.pyの内容確認:")
            with open('energy_conversion_HOPG.py', 'r', encoding='utf-8') as f:
                lines = f.readlines()
            
            # sの定義行を探して表示
            for i, line in enumerate(lines):
                if 's = [' in line and '自動更新' in line:
                    print(f"行 {i+1}: {line.strip()}")
                    break
        else:
            print("✗ 自動更新に失敗しました")
            
    except ImportError as e:
        print(f"エラー: インポートに失敗しました - {e}")
    except Exception as e:
        print(f"エラー: {e}")

if __name__ == "__main__":
    test_calibration_update()
