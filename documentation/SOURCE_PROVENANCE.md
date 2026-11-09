# Source provenance and limitations

Source scripts and the three reference CSV files are taken **verbatim** from the small code-only ZIP uploaded by the researcher, `JPR_CODE_ONLY_FOR_REPLICATION.zip`. No regression specifications were changed in those archived scripts. Original script versions can refer to absolute Windows paths and not all are runnable in a clean environment.

- `code/raw_upstream/` contains the original R00--R99 global build pipeline. The original source names are **R02_ingest_airbases.R** and **R03_ingest_ged.R**, not the incorrectly guessed names in earlier runners.
- `code/afghanistan/` contains the original Afghanistan-specific construction/estimation pipeline.
- `code/published_H2/` contains the original preferred H2 (joint ACLED-LWJ) and diagnostics code; many scripts require the saved FINAL30 panel and external CSVs.
- `reference/` holds the 30-airfield coordinates/whitelist and the original H1/H3 canonical results file.

## Unresolved critical link

No uploaded source R script writes `11_final_analysis_panel.rds` from scratch. The FINAL30 airfield construction step and its exact screening/decomposition algorithms have not been established from these files. Even an upstream raw-data re-run can therefore **not** currently be certified as regenerating the published FINAL30 sample.

`RUN_PUBLICATION.R` re-estimates the main published H1, joint H2, original ACLED-only H2 and H3 on the **published FINAL30 analysis panel**, thus testing the estimation tier. `RUN_RAW_UPSTREAM.R` attempts the older original raw-data pipeline and explicitly halts if the FINAL30 generation stage is unavailable. These are two separate tiers; do not conflate them.

Do not publicly share raw ACLED or other datasets without permission, nor claim full raw-to-final replication until the missing FINAL30 writer and all figure-generation steps are fully identified and validated.
