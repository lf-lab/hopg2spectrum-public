# サンプルデータ

このディレクトリには実験データ（CSVファイル）を配置してください。

## ファイル名形式
```
YYMMDD_GXXXXX_HOPG_HHMM.csv
```

- `YYMMDD`: 測定日（例: 240826）
- `GXXXXX`: ショット番号（例: G12345、L12345）
- `HHMM`: IP読み取り時刻（例: 1430）

## データ形式
CSVファイルは以下の形式である必要があります：
```
Position[cm], Intensity[PSL]
0.500, 100.25
0.501, 102.43
...
```

## 注意事項
- ヘッダー行は1行目に必要
- 位置は cm 単位
- 強度は PSL (Photo-Stimulated Luminescence) 単位
