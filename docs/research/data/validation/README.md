# Validation data extracts

Small extracts and summaries of open datasets, made for validating the smartphone battery model.
The full record and the reference keys are in `docs/research/validation_data.md`. The scripts are in
`docs/research/scripts/validation_*.py`. Each script docstring gives its usage and the source URL.

| File | Derived from | Licence of the source | Script |
|---|---|---|---|
| `nasa_pcoe_summary.json` | NASA Ames PCoE "Battery Data Set" (Saha & Goebel 2007), 34 cells, 4/24/43 °C | No licence stated; NASA asks users to acknowledge the repository and donors (see note) | `validation_nasa_summary.py` |
| `nasa_B00xx_*.csv`, `nasa_curves_index.json` | Same source. Six discharge curves (t, V, I, T_cell), decimated to 10 s; the square-wave curve is kept at native resolution | as above | `validation_nasa_extract_curves.py` |
| `lg_hg2_raw/*.csv` | Kollmeyer et al. 2020, LG HG2 3 Ah cell, Mendeley Data doi:10.17632/cp3473x7xv.3. **Unmodified copies** of 18 files (C/20, 0.5C and 2C discharges at 6 temperatures) | CC BY 4.0, attribution required | (copied as downloaded) |
| `lg_hg2_summary.json` | Same source. Q, E, mean V and T rise at -20/-10/0/10/25/40 °C, plus an apparent-resistance Arrhenius fit | CC BY 4.0 | `validation_lg_hg2_summary.py` |
| `oxford_ExampleDC_C1_artemis_discharge.csv`, `oxford_example_summary.json` | Birkl & Howey 2017, Oxford Battery Degradation Dataset 1, doi:10.5287/bodleian:KO2kdmYGg. Kokam 740 mAh pouch cell, Artemis drive-cycle discharge at 40 °C, 1 s logging | ODbL v1.0 / DbCL v1.0. These derived files are also offered under ODbL v1.0 | `validation_oxford_example.py` |
| `bonn_iphone_summary.json`, `bonn_iphone_summary_soc_curves.csv` | Schweizer et al. 2026, iPhone 13 Pro / 14 Pro battery-level traces, bonndata doi:10.60507/FK2/LSJGWW | CC BY 4.0 | `validation_bonn_iphone_summary.py` |
| `device_benchmarks.json` | Our arithmetic on published runtimes (Apple, Google, Notebookcheck, Laptop Mag) and the Bonn summary: P_avg = E/TTE, plus an additivity check | Facts with citation (not open data) | `validation_device_benchmarks.py` |

Notes

* **Sign convention.** In all exported CSVs, I > 0 means discharge.
* **NASA licence.** data.nasa.gov lists the NASA set as "License not specified" with public access. The NASA
  PCoE repository page asks publications to acknowledge the repository and the data donors. We treat it as a
  US-Government work. That status is background knowledge and was **not verified** in this session.
* **Raw files not kept.** The raw downloads (NASA .mat, Bonn JSON, Oxford .mat) are not stored in the repo. Re-download
  them with the URLs in the script docstrings. Then run each script with `python3 -I` from a directory outside
  the download folder.
* **Bonn device models are inferred.** The Bonn device-to-model mapping (13 Pro vs 14 Pro) is not stored in the
  data files. We infer it from idle runtime: the three shortest are taken as 13 Pro, matching the paper's 3 × 13 Pro
  + 6 × 14 Pro.
