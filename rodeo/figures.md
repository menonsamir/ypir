Figures to generate:
- `table-bit-retrieval`: Table comparing costs of single-bit retrieval
  - `python3 rodeo/plot.py table-bit-retrieval rodeo/data/*.json rodeo/latest_data/v2*`
  - `python3 rodeo/plot.py table-bit-retrieval --star-variants rodeo/data/*.json rodeo/latest_data/v2*`
- `bit-retrieval`: Graph of throughput for single-bit retrieval from various database sizes
  - `python3 rodeo/plot.py bit-retrieval --output-type tex rodeo/data/*.json rodeo/latest_data/v2*`
- `ypir-breakdown`: Table breaking down YPIR computation
- `ccb`: Graph of effective throughput with CCB
- `large-items`: Table comparing costs for large item retrieval
  - `python3 rodeo/plot.py large-items rodeo/data/*.json rodeo/latest_data/v2* rodeo/latest_data/extra/v1-large-items.json rodeo/data/extra/hintless-pir-real-large-items.json`
- `sct`: Table comparing costs for SCT auditing
- `comm-comp-tradeoff`: Graph showing communication-computation tradeoff
  - `python3 rodeo/plot.py comm-comp-tradeoff --output-type tex rodeo/data/*.json rodeo/latest_data/v2*`
