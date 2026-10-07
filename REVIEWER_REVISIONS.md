# Reviewer Revision Record

Updated: 2026-10-07

## Scope

This record covers the approved revisions through R11, including the Python
display repair and AI attribution. The identifiers refer to the project's
revision plan, not issue numbers assigned by reviewers. It summarizes the
changes and their disposition without reproducing private review correspondence.

This is an interim review checkpoint, not a declaration that the case study is
ready for scientific release. R12 and later revisions are not included. In
particular, the downstream expression-normalization correction in R13/F1
remains outstanding.

The selected model is unchanged: R CoGAPS, six patterns, seed 2, 2,000 iterations
per phase, and sparse optimization. The selected R IFN-associated label is
Pattern4; the parallel Python output uses Pattern5. Pattern numbers are
run-specific. No CoGAPS model was refitted for these revisions.

## Changes and Reviewer Responses

### R01: Starting Instructions

**Responds to:** Sarah and Hunter; setup findings in the full review.

Clarified basic R/Python prerequisites and reading-only access; added Docker
readiness checks, GitHub ZIP instructions, and project-root checks. Distinguished
Host Terminal, RStudio Terminal, RStudio Console, and Python session. Explained
their separate working directories, recommended a saved R script, identified
the included output folders, and collapsed only the optional full-reproduction
setup. Main setup remains visible.

**Files:** `_ocs_frontmatter.qmd`, `_packages.qmd`.

### Python Display Repair

**Responds to:** execution testing associated with R01/F7.

The shared display helper keeps formatted HTML during page rendering and prints
readable tables in an ordinary Python session. Plotting chunks save learner PNGs
under `learner_outputs/python/`, with explicit opening instructions. That folder
is ignored by Git; learner plots do not replace the scientific reference assets.

**Files:** `scripts/helpers/case_study_display.py`, Python chunks across the
learner modules, `index.qmd`, `.gitignore`.

### R02: Reproduction Archive

**Responds to:** Sarah's missing-source and empty-archive reports.

Linked the published, version-specific Zenodo record and documented individual
downloads, exact filenames, destinations, sizes, and SHA-256 checks. Explained
that creating a directory does not obtain the files. Added read-only,
route-aware input checking: optional Python models are not R-route prerequisites,
and fitting does not require a pre-existing result. Archived reference inputs
remain distinct from regenerated working outputs. Public file retrieval was
verified on October 5; this does not establish upstream annotation provenance.

**Files:** `_packages.qmd`, `_full_reproduction_guide.qmd`,
`data/external/README.md`, `data/external/reproduction_files.json`,
`scripts/reproduction/check_reproduction_inputs.py` and its tests.

### R03: Full Reproduction Execution

**Responds to:** Sarah's execution-context and environment failures.

Updated the documented runtime to `othomas2/pycogaps-runtime-guide:0.3.1`, which
includes the missing preprocessing dependency, `scikit-misc`. The preparation
script validates treatment labels and creates `condition` from `label` only when
needed. Added native R H5AD checks, explicit command locations, working-copy
warnings, and Docker availability/daemon guards. The R wrapper retains its
explicit overwrite protection. Explained the temporary fitting container and
return to RStudio for exports. Documented the figure builder's exact scope and
direct diagnostic-PNG workflow; learners do not rebuild the case-study page.

**Files:** `_full_reproduction_guide.qmd`, `_packages.qmd`, `README.md`,
`config_automation.yml`, reproduction scripts and tests.

### R04: Source and Annotation Boundary

**Responds to:** Ryan and earlier provenance requests from Carrie.

Identified the supplied annotated object as the Figshare/Pertpy distribution
linked to the Kang study. Corrected the claim that local preprocessing starts
directly from GEO. Added exact file identity and a provenance table, distinguished
inherited annotations/embeddings from local operations, and disclosed the
24,673-versus-24,305 cell-count discrepancy. Unknown upstream details are marked
"not yet known". Corrected donor wording: PBMC aliquots, not people, received
IFN-beta. Updated the preprocessing diagram and limitations accordingly.

**Files:** `_data_description.qmd`, `_data_wrangling.qmd`, `_limitations.qmd`,
`data/external/README.md`.

### R05: Filtering and the Modeled Gene Set

**Responds to:** Michael, Ryan, and Sarah.

Distinguished the minimum-three-cell gene filter, inherited cell-level QC, and
the practical choice of 3,000 HVGs. Explained that CoGAPS does not require this
number or HVG-only input and that excluded genes limit interpretation. Did not
assert an untested Poisson model, fit a 5,000-gene comparator, or change the model
input to address the separate downstream CPM issue. Gene-selection robustness
remains untested and is disclosed as a limitation.

**Files:** `_data_exploration.qmd`, `_data_wrangling.qmd`, `_limitations.qmd`.

### R06: Cell Counts Across Patterns

**Responds to:** Michael's question about differing counts and zero activities.

Clarified that counts differ between conditions, not between patterns: every
pattern includes 12,315 control and 12,358 stimulated cells. Added R/Python
checks for missing values, unique identifiers, six patterns, and agreement with
metadata counts. Zero activities are retained without an invented explanation
for their origin. This is a clarified and checked table, not a repaired filtering
bug.

**File:** `_data_wrangling.qmd`.

### R07: Diagnostics and Uncertainty

**Responds to:** Michael and Sarah; related to full-review finding F10.

Kept the full trace and added a separately scaled sampling-only view. Exported
phase/iteration labels only when saved metadata and the verified CoGAPS 3.22.0
schedule support reconstruction; unknown schedules remain unresolved. Original
trace measurements are unchanged. Removed plateau-as-convergence wording,
identified the snapshot overlap as an average across patterns, and distinguished
posterior activity SD from gene-loading CV. These changes do not certify
convergence or equal stability of all patterns.

**Files:** `_data_analysis.qmd`, `_full_reproduction_guide.qmd`,
`scripts/export_cogaps_diagnostics_r.R`, `scripts/reproduction/test_trace_timing_r.R`,
selected R trace/diagnostic metadata and README, rendered diagnostic plots.

### R08: Biological Wording and References

**Responds to:** Michael, Genevieve, Elana, and Sarah.

Clarified transcriptional regulation, cell-identity-associated structure,
transformed model input, run-specific pattern numbers, and the experiment's
clinical limits. Standardized SLE terminology while preserving article titles.
Retained Fertig 2010, placed the existing single-cell CoGAPS and Kang citations
where used, and added Johnson 2023 and one shared Stein-O'Brien 2017 reference
for the overlapping literature requests. Mentioning PatternMarkers does not
mean a new PatternMarkers analysis was performed.

**Files:** `_context.qmd`, `_motivation.qmd`, `_data_import.qmd`,
`_limitations.qmd`, `_cite.qmd`, terminology in related modules and glossary,
`resources/references.bib`.

### R09: Rank Choice and Optional Code Explanations

**Responds to:** Hunter and Ryan; author-directed presentation decisions.

Defined K and the actual tested grids; explained interpretability as well as
IFN recovery and the dataset-specific iteration trade-off. Kept model selection
in Data Import, rather than moving it to Data Analysis. Used the existing
thinking-question/disclaimer styles. Preserved the programming explanations in
optional language-specific foldouts across three subsections; code and results
remain visible, and no required object exists only inside a foldout.

**Files:** `_data_description.qmd`, `_data_import.qmd`, `_data_analysis.qmd`.

### R10: Contents of the Saved Objects

**Responds to:** Ryan.

Added a matrix-to-file reference, real R/Python inspection screenshots, text
inventories, and explanations of full-object contents. Explicitly distinguished
software/mathematical matrix orientation, uncertainty, diagnostics, inherited
inputs, and empty fields. Added a disabled, illustration-only R fitting call
checked against the installed package. The author declined a direct full-object
exercise because those large files are not part of the learner repository.
Learners continue with included compact outputs, without a new model run or
additional download.

**Files:** `_data_exploration.qmd`, `assets/object_overviews/` and rendered copies.

### R11: Motivation and Opening Explanation

**Responds to:** Hunter, Sarah, and Genevieve.

Explained the opening atlas's run-specific Pattern4 emphasis in Motivation,
defined PBMCs as a mixture, clarified the research purpose and clinical boundary,
and shortened repeated introductions. Kept the bulk-versus-single-cell
comparison and all substantive Context material. The opening image remains
standalone: no heading or caption was added around it.

Clarified that self nucleic acids discussed in the text are not explicitly
labeled in the reproduced SLE pathway diagram. Kept its artwork and attribution
unchanged and moved the redraw note out of the learner page. This adapts the
missing-label concern; it does not add labels to the artwork. The pathway
adaptation and opening-image internal-title decisions remain separate.

**Files:** `_motivation.qmd`, `_context.qmd`.

### AI Attribution

**Requested by:** the author.

Added the approved acknowledgment describing human-directed, reviewed Codex
and ChatGPT assistance, the stated model list, and the distinct use of ChatGPT
for illustrations and code rather than earlier narrative drafts.

**File:** `_acknowledgements.qmd`.

## Validation

### PR Preparation: Local-Only RStudio Access

During PR review, the author approved restricting both Docker launch examples
and the alternate-port troubleshooting example to `127.0.0.1`. RStudio remains
available at the same localhost URL, but the published example credentials are
not exposed through all host interfaces. This is a setup safety correction, not
R12 or an analysis change.

**Files:** `_packages.qmd`, `_full_reproduction_guide.qmd`.

Each checkpoint was rendered and checked locally. For this PR, validation uses
the documented Docker runtime and does not fit CoGAPS models or download the
optional full model files.

- 36 Python reproduction-script unit tests pass.
- 27 R trace-reconstruction checks pass, including missing/mismatched metadata.
- Final clean-archive render, Python terminal workflow, asset/link checks, and
  repository packaging checks are being completed before submission.

The unit tests can be run from the project root in the validated container:

```bash
/opt/cogaps-venv/bin/python -m unittest discover -s scripts/reproduction -p 'test_*.py' -v
Rscript scripts/reproduction/test_trace_timing_r.R
```

Page rendering is a maintainer validation step, not an added learner instruction.

## Not Resolved by This PR

- **R12:** the UMAP/PCA explanation and additional metadata views have not started.
- **R13/F1:** directionality still uses the 3,000-HVG denominator. Its normalization
  and affected scientific interpretations require correction before final release.
- **F3:** upstream cell filtering, annotations, embedding settings, and the cell-count
  discrepancy are not resolved by distributing an identified source object.
- **Python seed provenance:** the saved object's internal seed is 0 while its
  filename/run record says 2. The discrepancy is disclosed, not silently fixed.
- **Broader science:** all-pattern stability, technical-covariate assessment, and
  later enrichment/PatternMarkers revisions remain separate work.
- **Site QA:** Motivation and the top-level Context introduction are missing from
  the current generated search index. Fresh interactive desktop/mobile browser
  checks remain pending after the local-file browser-access limitation.
- **Artwork:** no opening-title change or pathway-diagram adaptation is included.
- **Jessica's review:** it will be addressed separately, as requested by the author.

Review readiness and successful rendering are not scientific release approval.
