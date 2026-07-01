# External or Local-Only Artifacts

The case study is designed to render from GitHub-safe files. Large full-local files should be supplied separately only when a reviewer, instructor, or learner wants to regenerate selected-model derivatives or inspect the complete model object.

Recommended optional paths:

- `data/processed/input/cogaps_input_genesxcells_hvg3000_float64.h5ad`
- `data/processed/selected_model_k6/r/cogaps_K6_seed2_iter2000.rds`
- `data/processed/selected_model_k6/r/cogaps_K6_seed2_iter2000.checkpoint.out`
- `data/processed/selected_model_k6/python/cogaps_K6_seed2_iter2000.h5ad`

These files may be hosted on Zenodo or another artifact host, or kept local for private full reproduction. They should not be committed to ordinary GitHub history.

The full HPC sweep output archive, older large K7 model objects, and developer-facing runtime dashboards are not required for the learner path or the selected-model full reproduction guide.
