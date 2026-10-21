---
name: qdi-texture-qaqc
description: Read-only QA/QC of the MGS Quaternary Data Index texture/grain-count database (QDI "Texture Data" layer, view mgs_qdi.qdi.vw_qdtx). Finds mistyped or unknown stratunit codes (Quaternary and bedrock), rare typo-like codes, percentages out of range or not summing to 100, values that disagree with their lab weights or grain counts, bad depths/elevations, duplicates and orphan rows, and writes proposed changes to CSVs without touching the database. Use when asked to check, audit, QA/QC or find errors in QDI texture data, strat codes, or the Texture Data .lyr.
---

# QDI texture QA/QC

## Ground rules
- **Never write to the database.** Read only through `arcpy.da.SearchCursor`. Never use update or insert cursors or edit sessions, and never call `LayerFile.save()`. Every finding goes into CSVs of *proposed* changes.
- **Source:** `L:\gdrs\data\pub\edu_umn_mngs\geos_quaternary_data_index\Texture Data.lyr`. It is a PostgreSQL SDE query layer, logged in as the read-only user `qdiviewer`, on the view `mgs_qdi.qdi.vw_qdtx` (view SQL: `...\metadata\SQL text for QDI_Views.txt`).
- **The published layer hides rows.** Its object-ID field is `relateid`, which isn't unique, so ArcGIS shows only 71,396 of the 71,903 rows (2026-10-15). The hidden rows are the ones with a NULL relateid. `scripts/qdi_load.py` switches the in-memory layer to `objectid` so every row is read. Report this every time.
- **Units:** elevation and depths are in **feet**, and weights are in grams. Column meanings are in `references/field_dictionary.csv`, as supplied by the user.
- **Strat codes:** check **every** stratunit code, bedrock included, against `references/stratunit_codes.tsv`, the full MGS table. A blank or NULL stratunit means "not classified yet" and is skipped by the strat checks.
- **A strat "outlier" is a rare code that looks like a typo, not a composition outlier.** Do not use texture or clast composition to judge a stratunit.

## Run
Use ArcGIS Pro's Python, because it has arcpy and pandas.  It can typically be found at:
```
"C:\Program Files\ArcGIS\Pro\bin\Python\envs\arcgispro-py3\python.exe" ~/.claude/skills/qdi-texture-qaqc/scripts/run_qaqc.py ^
    --out PROMPT_USER_FOR_OUTPUT_PATH --refresh
```
- `--refresh` re-reads the database, which takes about 30 s. Without it, the script uses the pickle cache (`--cache`).
- `--lyr` points at another copy of the layer, and `--codes` at another code table.

Supplemental CWI data used for identifying missing data can be found in the Microsoft access database at `T:/db_prod/cwi/data/cwidata.mdb`.  The 
relevant tables are 
- `c4ix` Main borehole data table containing summary and spatial information
- `c4c1` Well construction information

## Process
1. **Load.** Read through the layer file with the in-memory OID fix, cache the data, and record the row counts and max(date_modify) in `run_info.csv`.
2. **Profile before you trust a rule.** For every derived percentage, `run_qaqc.py` measures how often the formula reproduces the stored value. That agreement is written to `run_info.csv`:
   - Below 97%: the formula is not used.
   - 97% to 99.5%: mismatches are reported as `low`.
   - 99.5% or more: mismatches are reported as `medium`, with the recomputed value proposed.

   If you add a rule, verify it the same way first.
3. **Strat checks** (`strat_codes.csv`):
   - `STRAT_UNKNOWN_CODE`: the code is not in the table. The proposal lists the valid codes one keystroke away, with letter swaps first (QSHU→QHSU) and then the most-used.
   - `STRAT_RARE_TYPO`: a valid code used on ≤ 10 rows that is one keystroke from a code of a *different formation* used at least 10× as often. It is `medium` when the same borehole uses the suggested code, otherwise `low`. Codes that differ only in the last letter (member or provenance of the same formation) are never typo suspects.
   - `STRAT_CONTEXT_TYPO` / `STRAT_CONTEXT_MEMBER`: a sample whose code is one keystroke from the code shared by the samples directly above and below it in the same borehole.
   - `STRAT_WHITESPACE` and `STRAT_FORMAT`: stray spaces or wrong case.
   - `STRATFORMATION_*` and `STRATMEMBER_UNKNOWN`: invalid codes, or a stratformation that doesn't match the stratunit's formation.
4. **Range and sum checks** (`ranges_sums.csv`):
   - `PCT_OUT_OF_RANGE`: a percentage outside 0 to 100.
   - `GROUP_SUM_NOT_100`: the groups in `SUM_GROUPS` don't add up to 100. The allowed slack is k/2 for k rounded parts.
   - `DERIVED_PCT_MISMATCH`: a stored % differs from the % recomputed from weights or counts.
   - `TEXTURE_ZERO_*`: "no data" stored as zeros.
   - Weights: `NEGATIVE_WEIGHT`, `GRAVEL_GT_TOTAL`, `HEAVY_SAMPLE` (> 500 g), `ADJUSTWEIGHT_MISMATCH` (≠ totalweight − gravel), `SSC_WEIGHTS_NE_ADJDRY`, `SAND_CLAY_WEIGHT_MISMATCH`.
   - Counts: `COUNT_ROLLUP_MISMATCH` (a total ≠ the sum of its parts) and `COUNT_COMPONENTS_MISSING`.

   Recomputed values are proposed only when the input weights are sane.
5. **Value domains** (`domains.csv`):
   - deposittype: digits or empty strings. The user said "T or L", but the data also uses T*, O, F, C, G, … Those are listed in `field_value_counts.csv` for the user to rule on and are not flagged.
   - leached must be Y or N.
   - type: wrong case.
   - Dry and wet colors that don't parse as Munsell.
6. **Depth and elevation** (`depth_elev.csv`):
   - Sampledepth of −99 (placeholder for null) or negative.
   - Depths over 1,000 ft.
   - Depth_top below depth_bottom (a swap is proposed).
   - Sampledepth outside its interval, checked only where depth_top > 0 because 0 is a placeholder.
   - depth_unit is `m` or NULL.
   - Borehole elevation outside 600–2,310 ft, reported once per relateid in table `qdix`. A value of 0 or −1 means missing; values of 270–360 look like metres.
   - Use `c4ix` to try and extract the correct elevation value via shared RelateID.
7. **Integrity** (`integrity.csv`):
   - `NULL_RELATEID`.
   - `RELATEID_PLACEHOLDER` (00Q00 and similar).
   - `RELATEID_FORMAT`: the id is not 10 characters. Zero-padding is proposed only when the padded id exists.
   - `NO_LOCATION`: the relateid has no qdix record.  Search the `c4ix` table to see if an entry with the same RelateID 
      has location data.  If it does, suggest that location data as values.
   - `DUPLICATE_EXACT` (proposed: delete) and `DUPLICATE_CONFLICT` (lists the columns that differ). Duplicates share the same relateid, sampledepth and sample_num.
8. **Review.** Before reporting, read about 5 flags from each check. If a check returns thousands of hits, look for a placeholder or convention behind them, as with depth_top = 0 and sample_num = borehole id, and tighten the check instead of dumping the noise on the user.
9. **Verify read-only.** After the run, re-read the layer and confirm the row count and max(date_modify) are unchanged and the .lyr file's timestamp is unchanged.
10. **Identify missing entries.**  Use `c4c1` to find boreholes that were drilled using Rotary Sonic (coded as `V`) that do not have entries in QDI; they will share a RelateID.

## Output columns
`check, category, severity (high/medium/low), objectid, relateid, sample_num, sampledepth, table (qdtx|qdix), field, current_value, proposed_value, reason`

- `all_proposed_changes.csv` holds everything.
- `summary.csv` has counts per check.
- `stratunit_code_counts.csv` lists every code used, with its row count and whether it is in the table.
- `field_value_counts.csv` lists the value domains.
- `run_info.csv` has the layer facts and formula agreement.

An empty `proposed_value` means the right value can't be inferred, so a person has to look at the lab sheet or log.
