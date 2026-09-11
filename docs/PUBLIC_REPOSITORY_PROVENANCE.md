# Public repository provenance

This GitHub-facing tree was derived from the final pre-submission reproducibility code package prepared for the 2026 pentagon-array manuscript. The scientific source files were preserved while the packaging was reorganized for public readability.

Public-facing changes are intentionally non-scientific:

- removed `.DS_Store` / `__MACOSX` metadata;
- excluded raw DAQ data (not present in the code package in the first place);
- excluded large generated `.npz` and raw shower `.parquet` products;
- excluded the frozen full-byte third-party catalog ZIP while retaining its provenance manifest;
- added portable path arguments to a public orchestration driver;
- added README, environment/requirements files, data/licensing notes, issue guidance, static QA, and CI metadata.

The authoritative scientific analysis remains `experimental/analysis/Pentagon_Array_Directional_Analysis_V12.py`; the reorganization does not alter that file.
