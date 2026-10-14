# External or Local-Only Artifacts

The case study is designed to render from GitHub-safe files. Large full-local files should be supplied separately only when a reviewer, instructor, or learner wants to regenerate selected-model derivatives or inspect the complete model object.

## Availability

The reproduction files are available in [version-specific Zenodo record 21208760](https://zenodo.org/records/21208760), DOI [10.5281/zenodo.21208760](https://doi.org/10.5281/zenodo.21208760). Complete anonymous downloads of all ten files were verified against the expected byte sizes, MD5 digests, and SHA-256 digests on October 5, 2026. The included files support the main learner path without these additional downloads.

This confirms public retrieval and file identity, not the full upstream preparation history or every reproduction route. Owner confirmation of release metadata and redistribution/licensing remains a separate release check. Use the fixed record above rather than a draft sharing link or a concept DOI that may resolve to another version.

## Choose Only the Inputs You Need

| Route | Required inputs | Generated outputs, not prerequisite downloads |
| --- | --- | --- |
| Main learner path | Included prepared data and compact summaries | No archive files required |
| Source preprocessing | `kang_counts_25k.h5ad` | Prepared data, preprocessing record, dense input |
| R fitting | Prepared data, already included | Full R result and checkpoint |
| R export/full-object inspection | Prepared data and full R result | Compact tables/diagnostics |
| Optional Python fitting | Prepared data and dense input | Full Python result |
| Python full-object inspection | Prepared data and full Python result | None |

A full Python model is not an R-route prerequisite. The checkpoint is optional provenance; resuming from it is not a validated route in this revision.

## Layout and Safe Placement

The current deposit contains individual files, not a directory-preserving ZIP. If a future release uses an archive container, extract it first. Keep the obtained files in a separate reference folder; creating that folder does not download anything. Mount it read-only at `data/external/reproduction_archive/` inside Docker, and use a separate clone or ZIP extraction as the writable reproduction project.

Paths below are relative to that working project's root (`index.qmd`, `scripts/`, `data/`). Copy only what your route requires. Do not replace original reference files or assume a subfolder called `reproduction_work/` redirects the model scripts.

| Archive filename | Working-project destination |
| --- | --- |
| `kang_counts_25k.h5ad` | `data/source/kang_counts_25k.h5ad` |
| `preprocessed_cells_hvg3000.h5ad` | `data/processed/input/preprocessed_cells_hvg3000.h5ad` |
| `preprocess_config_hvg3000.json` | `data/processed/input/preprocess_config_hvg3000.json` |
| `cogaps_input_genesxcells_hvg3000_float64.h5ad` | `data/processed/input/cogaps_input_genesxcells_hvg3000_float64.h5ad` |
| `cogaps_K6_seed2_iter2000.rds` | `data/processed/selected_model_k6/r/cogaps_K6_seed2_iter2000.rds` |
| `cogaps_K6_seed2_iter2000.checkpoint.out` | `data/processed/selected_model_k6/r/cogaps_K6_seed2_iter2000.checkpoint.out` |
| `cogaps_K6_seed2_iter2000_python.h5ad` | `data/processed/selected_model_k6/python/cogaps_K6_seed2_iter2000.h5ad` |

Notice the Python result's filename change: remove `_python` when placing it in the `python/` directory. Do not download a full result before fitting merely to satisfy a generic file list; the launchers refuse to overwrite existing full-model files.

## Verify Before Running a Stage

`reproduction_files.json` records the exact byte sizes and SHA-256 digests of the seven audited reference files. It is an integrity index, not a download mechanism or a replacement for `data/artifact_manifest.json`. Its `public_access_verified` field and dated verification details record the complete anonymous downloads verified on October 5, 2026. This confirms public retrieval and file identity, not upstream provenance, redistribution approval, or end-to-end reproduction; the indexed file identities are unchanged.

In the **RStudio Terminal inside the container**, from `/home/rstudio/project`, run:

```bash
/opt/cogaps-venv/bin/python scripts/reproduction/check_reproduction_inputs.py \
  --route preprocess
```

Use `fit-r`, `export-r`, `fit-python`, or `inspect-python` for the other stages; `learner-r` and `learner-python` check their starting files. The same command supports both languages. It checks readable, nonempty working inputs and the size/hash of original archive inputs. When a required working file is absent or invalid, it checks for a valid flat-file archive copy and reports where to copy it. It never copies, writes, installs, downloads, or fits anything.

Missing/invalid required inputs and existing full-model output collisions cause a nonzero exit. An empty archive is explicitly identified, but does not block a route whose inputs are already available. Files not used by the selected stage are marked `not needed for this route`; files the stage produces are marked `output to be generated`.

After successfully generating and validating an input yourself, name it explicitly, for example `--route fit-r --generated-input prepared`, or `--route export-r --generated-input prepared --generated-input r_model`. Readability and nonzero size are still required, but those declared outputs are not compared with the original model/data checksum. The source file cannot be exempted this way. This option is not a substitute for validating generated object content or a workaround for a corrupted archive download.

The full HPC sweep output archive, older large K7 model objects, and developer-facing runtime dashboards are not required for the learner path or the selected-model full reproduction guide.
