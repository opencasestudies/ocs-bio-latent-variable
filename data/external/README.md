# Reproduction Downloads

The main learner path uses files already included in this repository. For additional reproduction files, download only what you need from [Zenodo release 21208760](https://zenodo.org/records/21208760). The files are provided individually, not as a ZIP.

| Task | Additional download | Destination folder |
| --- | --- | --- |
| Rebuild the prepared data | `kang_counts_25k.h5ad` | `data/source/` |
| Recreate R tables from the saved model | `cogaps_K6_seed2_iter2000.rds` | `data/processed/selected_model_k6/r/` |
| Fit a Python model | `cogaps_input_genesxcells_hvg3000_float64.h5ad` | `data/processed/input/` |
| Inspect the saved Python model | `cogaps_K6_seed2_iter2000_python.h5ad` | `data/processed/selected_model_k6/python/` |

Destination folders are relative to the project root. Keep downloaded originals unchanged and copy the needed files into these folders. Rename the Python model's working copy to `cogaps_K6_seed2_iter2000.h5ad`.

Fitting an R model uses the included prepared data; no saved-model download is needed. The Python input does not need downloading if you generated it during preprocessing.

See the [Full Reproduction Guide](../../_full_reproduction_guide.qmd) for setup, file checks, and overwrite precautions. Recorded file sizes and checksums are in [reproduction_files.json](reproduction_files.json).
