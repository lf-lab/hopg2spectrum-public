# HOPG X-ray Spectrometer Data Analysis System (Public Version)

A collection of Python programs for analyzing and calibrating X-ray spectroscopy measurement data using HOPG (Highly Ordered Pyrolytic Graphite).

## 📝 Citation and Acknowledgment

**If you use this software in your research, please include the following citation and acknowledgment:**
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16946066.svg)](https://doi.org/10.5281/zenodo.16946066)

### Recommended Citation
```
HOPG X-ray Spectrometer Data Analysis System (Public Version)
Ryunosuke Takizawa
Institute of Laser Engineering, Osaka University
Available at: https://github.com/lf-lab/hopg2spectrum-public
```

### Acknowledgment
```
This research utilized the HOPG X-ray spectrometer data analysis tools 
developed by Ryunosuke Takizawa at the Institute of Laser Engineering, Osaka University.
```

**For academic publications, please also consider citing relevant papers on HOPG-based X-ray spectroscopy methods and laser fusion diagnostics from the Institute of Laser Engineering, Osaka University.**

## ⚠️ Important Notice (Public Version)

This is the public release version. The following features require access to specific institutional networks:

- **Automatic time delay calculation**: Requires connection to internal web database
- **Automatic shot time retrieval**: Only available in specific network environments

General users should manually specify the time delay using the `-t` option.

## Overview

This project provides automated tools for X-ray energy measurement data processing in laser fusion experiments, utilizing HOPG diffraction-based X-ray spectrometers. It converts position-intensity data acquired from Imaging Plates (IP) into energy-absolute photon number density spectra.

## Key Features

- **Automatic Energy Calibration**: Automated calibration using known reference lines
- **Interactive Calibration**: GUI-based manual selection and fine-tuning of reference line positions
- **Automatic Time Delay Calculation**: Automatic time delay calculation from filenames and web database (requires institutional network access)
- **Filter Correction**: Transmittance correction for beryllium and polyethylene filters
- **Graph Output**: Automatic generation of high-quality PDF graphs
- **Data Export**: CSV format spectrum data output

## File Structure

```
hopg2spectrum/
├── energy_conversion_HOPG.py              # Main analysis program
├── energy_conversion_HOPG_calibration.py  # Calibration-specific program
├── filter/                                # Filter transmittance data
│   ├── Be_500um_transmittance.dat         # Beryllium filter (500μm)
│   └── CH2_50um_transmittance.dat         # Polyethylene filter (50μm)
├── profile/                               # IP scan data directory
│   └── README.md                          # Data format instructions
├── data/                                  # Analysis results output directory
│   └── README.md                          # Output format description
└── README.md                              # This file
```

## Required Libraries

```bash
pip install numpy scipy matplotlib tkinter requests beautifulsoup4
```

## Usage

### 1. Standard Analysis (energy_conversion_HOPG.py)

Execute rapid analysis using existing calibration parameters:

```bash
# Basic usage (automatic time delay calculation)
python energy_conversion_HOPG.py

# Manual time delay specification (recommended for general users)
python energy_conversion_HOPG.py -t 23

# Interactive calibration mode
python energy_conversion_HOPG.py -c
```

**Operation Steps:**
1. Select IP scan data file (CSV file) via GUI
2. Extract shot number and calculate time delay (or manual input)
   - **General users**: Always specify time delay using `-t` option
3. Execute energy conversion and filter correction
4. Output results as CSV file and PDF graph

### 2. Calibration (energy_conversion_HOPG_calibration.py)

Set new calibration parameters:

```bash
# Execute calibration (requires institutional network access)
python energy_conversion_HOPG_calibration.py

# Manual time delay specification
python energy_conversion_HOPG_calibration.py -t 25
```

**Operation Steps:**
1. Select IP scan data file
2. Interactively select two reference line positions on the graph
3. Fine-tune positions by dragging
4. Calculate new calibration parameters
5. Automatically update main program parameters

## Data File Formats

### Input Data (IP Scan Data)
- Filename format: `YYMMDD_GXXXXX_HOPG_HHMM.csv`
  - `YYMMDD`: Measurement date
  - `GXXXXX or LXXXX`: Shot number
  - `HHMM`: IP reading time
- CSV format: `Position[cm], Intensity[PSL]`

### Output Data (Absolute Photon Number Spectrum)
- Filename: `HOPG_GXXXXX.csv`
- CSV format: `Energy[keV], Absolute_Photon_Number_Density[photons/keV]`

## Technical Specifications

### Calibration Method
- **Reference Energy Lines**: Cu-Kα line (8.048 keV) and He-α line (8.391 keV)
- **Lattice Spacing**: d = 3.357 Å
- **Crystal Reflectivity**: 5±1E-08 str (Calibrated by R.Takizawa 2021/08/24)

### Correction Processes
- **Time Decay Correction**: Accounts for time-dependent IP sensitivity
- **Energy Response Correction**: Corrects for IP energy response characteristics
- **Filter Transmittance Correction**: Removes effects of beryllium and polyethylene filters
- **Crystal Reflectivity Correction**: Accounts for HOPG crystal reflectivity

### Automation Features
- **Automatic Time Delay Calculation**: Extracts reading time from filename and retrieves shot time from web database
- **Automatic Laser Type Detection**: Automatically determines GXII/LFEX laser from shot number
- **Automatic Directory Creation**: Creates output directories if they don't exist

## Interactive Features

### Calibration Position Selection
- Click-based reference line position selection on graphs
- Drag-based position fine-tuning functionality
- Real-time position display and parameter updates

### File Selection
- Initial directory set to `profile/` folder
- Automatic CSV file filtering
- Intuitive GUI-based file selection

## Error Handling

- **Interactive Mode Failure**: Automatically falls back to maximum value position selection
- **Web Connection Errors**: Switches to manual time delay input
- **File Reading Errors**: Displays error messages and continues processing appropriately

## Development History

- August 2024: Initial version development
- August 2024: Automatic directory creation feature added
- August 2024: Interactive calibration feature added
- August 2024: Code refactoring and feature separation

## License

Intended for research use. Please contact us before commercial use.

## Authors

Institute of Laser Engineering, Osaka University

## Important Notes

- Automatic shot time retrieval requires connection to the Laser Science Institute network
- Please save IP scan data in the specified format
- Always test with sample data after calibration

## 概要

このプロジェクトは、レーザー核融合実験におけるX線エネルギー測定において、HOPGの回折を利用したX線分光器のデータ処理を自動化するためのツール群です。IP（Imaging Plate）で取得した位置-強度データをエネルギー-絶対フォトン数密度スペクトラムに変換します。

## 主要機能

- **自動エネルギー校正**: 既知の基準線を用いた自動キャリブレーション
- **対話的キャリブレーション**: GUIによる基準線位置の手動選択・微調整
- **時間遅延自動計算**: ファイル名とWebデータベースから自動的に時間遅延を算出 (レーザー研内部にて動作可能)
- **フィルター補正**: ベリリウムとポリエチレンフィルターの透過率補正
- **グラフ出力**: 高品質なPDFグラフの自動生成
- **データ保存**: CSV形式でのスペクトラムデータ出力

## ファイル構成

```
hopg2spectrum/
├── energy_conversion_HOPG.py              # メインの解析プログラム
├── energy_conversion_HOPG_calibration.py  # キャリブレーション専用プログラム
├── filter/                                # フィルター透過率データ
│   ├── Be_500um_transmittance.dat         # ベリリウムフィルター (500μm)
│   └── CH2_50um_transmittance.dat         # ポリエチレンフィルター (50μm)
├── profile/                               # IPスキャンデータ置き場
│   └── *.csv                              # 実験データファイル
├── data/                                  # 解析結果出力先
│   ├── *.csv                              # 絶対フォトン数スペクトラムデータ
│   └── *.pdf                              # スペクトラムグラフ
└── README.md                              # このファイル
```

## 必要なライブラリ

```bash
pip install numpy scipy matplotlib tkinter requests beautifulsoup4
```

## 使用方法

### 1. 通常の解析（energy_conversion_HOPG.py）

既存のキャリブレーションパラメータを使用して迅速に解析を実行：

```bash
# 基本的な使用法（自動時間遅延計算）
python energy_conversion_HOPG.py

# 時間遅延を手動指定（推奨 - 一般ユーザー向け）
python energy_conversion_HOPG.py -t 23

# インタラクティブキャリブレーションモード
python energy_conversion_HOPG.py -c
```

**動作手順:**
1. GUIでIPスキャンデータファイル（CSVファイル）を選択
2. ショット番号と時間遅延を自動計算（または手動入力）
   - **一般ユーザー**: `-t`オプションで時間遅延を必ず指定してください
3. エネルギー変換とフィルター補正を実行
4. 結果をCSVファイルとPDFグラフで出力

### 2. キャリブレーション（energy_conversion_HOPG_calibration.py）

新しいキャリブレーションパラメータを設定：

```bash
# キャリブレーション実行　(レーザー研内部にて動作可能)
python energy_conversion_HOPG_calibration.py

# 時間遅延を手動指定
python energy_conversion_HOPG_calibration.py -t 25
```

**動作手順:**
1. IPスキャンデータファイルを選択
2. グラフ上で2つの基準線位置を対話的に選択
3. ドラッグによる位置の微調整が可能
4. 新しいキャリブレーションパラメータを計算
5. メインプログラムのパラメータを自動更新

## データファイル形式

### 入力データ（IPスキャンデータ）
- ファイル名形式: `YYMMDD_GXXXXX_HOPG_HHMM.csv`
  - `YYMMDD`: 測定日
  - `GXXXXX or LXXXX`: ショット番号
  - `HHMM`: IP読み取り時刻
- CSV形式: `位置[cm], 強度[PSL]`

### 出力データ（絶対フォトン数スペクトラム）
- ファイル名: `HOPG_GXXXXX.csv`
- CSV形式: `エネルギー[keV], 絶対フォトン数密度[photons/keV]`

## 技術仕様

### キャリブレーション方式
- **基準エネルギー**: Cu-Kα線 (8.048 keV) と He-α線 (8.391 keV)
- **格子間隔**: d = 3.357 Å
- **結晶反射率**: 5±1E-08 str (Calibrated by R.Takizawa 2021/08/24)

### 補正処理
- **時間減衰補正**: IP感度の時間依存性を考慮
- **エネルギー応答補正**: IPのエネルギー応答特性を補正
- **フィルター透過率補正**: ベリリウムとポリエチレンフィルターの影響を除去
- **結晶反射率補正**: HOPG結晶の反射率を考慮

### 自動化機能
- **時間遅延自動計算**: ファイル名から読み取り時刻を抽出し、Webデータベースからショット時刻を取得
- **レーザータイプ自動判定**: ショット番号からGXII/LFEXレーザーを自動判別
- **ディレクトリ自動作成**: 出力ディレクトリが存在しない場合は自動作成

## 対話的機能

### キャリブレーション位置選択
- グラフ上でのクリックによる基準線位置選択
- ドラッグによる位置の微調整機能
- リアルタイムでの位置表示とパラメータ更新

### ファイル選択
- 初期ディレクトリが`profile/`フォルダに設定
- CSVファイルの自動フィルタリング
- GUIによる直感的なファイル選択

## エラーハンドリング

- **インタラクティブモード失敗時**: 自動的に最大値位置をフォールバック選択
- **Web接続エラー**: 手動での時間遅延入力に切り替え
- **ファイル読み込みエラー**: エラーメッセージの表示と適切な処理継続

## 開発履歴

- 2024年8月: 初期バージョン開発
- 2024年8月: 自動ディレクトリ作成機能追加
- 2024年8月: 対話的キャリブレーション機能追加
- 2024年8月: コードリファクタリングと機能分離

## ライセンス

研究用途での利用を想定しています。商用利用については事前にご相談ください。

## 作者

大阪大学レーザー科学研究所
瀧澤龍之介

## 注意事項

- ショット時刻の自動取得にはレーザー科学研究所のネットワークへの接続が必要です
- IPスキャンデータは指定された形式で保存してください
- キャリブレーションは必ず基準線が見えるショットをしてください
