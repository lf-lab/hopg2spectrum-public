# Sample Data Directory

Please place experimental data (CSV files) in this directory.

## Filename Format
```
YYMMDD_GXXXXX_HOPG_HHMM.csv
```

- `YYMMDD`: Measurement date (e.g., 240826)
- `GXXXXX`: Shot number (e.g., G12345, L12345)
- `HHMM`: IP reading time (e.g., 1430)

## Data Format
CSV files must be in the following format:
```
Position[cm], Intensity[PSL]
0.500, 100.25
0.501, 102.43
...
```

## Notes
- Header row is required on the first line
- Position in cm units
- Intensity in PSL (Photo-Stimulated Luminescence) units
