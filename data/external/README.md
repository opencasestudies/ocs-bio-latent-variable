# External or Local-Only Artifacts

The case study is designed to render from GitHub-safe files. Large full-local files should be supplied separately only when a reviewer, instructor, or learner wants to regenerate selected-model derivatives or inspect the complete model object.

## Availability

The reproduction files are available in [version-specific Zenodo record 21208760](https://zenodo.org/records/21208760), DOI [10.5281/zenodo.21208760](https://doi.org/10.5281/zenodo.21208760). Complete anonymous downloads of all ten files were verified against the expected byte sizes, MD5 digests, and SHA-256 digests on October 5, 2026. The included files support the main learner path without these additional downloads.

This confirms public retrieval and file identity, not the full upstream preparation history or every reproduction route. Owner confirmation of release metadata and redistribution/licensing remains a separate release check. Use the fixed record above rather than a draft sharing link or a concept DOI that may resolve to another version.

## Source Object and Provenance

The source study is [Kang et al. (2018)](https://doi.org/10.1038/nbt.4042), with public study data in [GEO GSE96583](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583). The local reproduction scripts do not reconstruct an annotated object from those GEO files. They begin with the supplied `kang_counts_25k.h5ad`.

The older project acquisition instructions identify Figshare file `34464122`, in [record 19397624](https://figshare.com/articles/dataset/Kang_HM_Subramaniam_M_Targ_S_Nguyen_M_et_al_2017/19397624), version 2 dated March 22, 2022. The record credits **Antonia Schumacher and Lukas Heumos** and links the Kang article and GEO accession. These are dataset-deposit credits, not proof of who performed every annotation or conversion step. [Pertpy 1.0.2's loader](https://pertpy.scverse.org/en/1.0.2/_modules/pertpy/data/_datasets.html#kang_2018) uses that Figshare file. The local source is also identical to a separately downloaded copy from the current Pertpy/scverse distribution, documented in the [current loader](https://pertpy.readthedocs.io/en/stable/_modules/pertpy/data/_datasets.html#kang_2018). No successful fresh direct Figshare download is claimed by that comparison.

This case study's Zenodo record preserves a copy of that distributed source:

- Filename: `kang_counts_25k.h5ad`.
- Size: 38,356,412 bytes.
- SHA-256: `e6a5adac64dcdeb36eaba27db49b63e0c64bb0ed4a64c6705971506b41c39830`.
- Dimensions: 24,673 cells by 15,706 genes.

The [local preprocessing script](../../scripts/reproduction/prepare_preprocessed_hvg3000.py) filters genes, normalizes and log-transforms expression, and selects 3,000 highly variable genes using the retained counts. It retains the supplied cells, donor `replicate` labels, `cell_type` labels, and `X_pca`/`X_umap` coordinates. It copies the supplied `label` values into `condition` rather than inferring treatment assignments. The prepared dataset therefore contains the same 24,673 cells, with 3,000 modeled genes.

The complete GEO-to-H5AD workflow, exact earlier cell/gene filtering, donor-label assignment history, cell-type annotation method and responsible annotator(s), embedding inputs/settings, and conversion software/versions are **not yet known**. In particular, the source contains 368 more cells than the paper's reported 24,305 IFN-experiment singlets (12,138 control plus 12,167 stimulated). The reason is **not yet known**; local gene filtering does not account for it. Neither public distribution nor matching checksums resolves these upstream questions, and the supplied object should not be described as the paper's exact reported singlet set without further evidence.

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

Use Docker image `othomas2/pycogaps-runtime-guide:0.3.1` and the updated repository scripts. Source preprocessing and native R import have been tested in that runtime. A separate saved-R-result walkthrough also verified diagnostic/compact exports and figure generation without fitting a model. These checks do not certify a new full model fit, every optional Python/HPC route, or every platform. Successful downloads and local reproduction do not independently validate the inherited annotations or reconstruct their upstream history.

The full HPC sweep output archive, older large K7 model objects, and developer-facing runtime dashboards are not required for the learner path or the selected-model full reproduction guide.
