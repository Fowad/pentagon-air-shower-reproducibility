# Repository QA summary

This public-safe repository was prepared from the final pre-submission reproducibility code package used for the pentagon-array manuscript.

The source package had already passed the following static checks before this GitHub-facing reorganization:

- every packaged Python source parsed/compiled;
- every `.command` launcher passed `bash -n` syntax validation;
- packaged JSON files parsed successfully;
- frozen catalog ZIP integrity was checked;
- the scientific branches were designed to fail closed rather than silently substitute historical data products.

For this public repository, the checks are repeated by `scripts/repository_smoke_check.py` and by `.github/workflows/static-checks.yml`. Raw DAQ files, raw shower products, large Monte Carlo arrays, macOS metadata, and third-party catalog byte snapshots have been excluded.

A static check does **not** replace the full scientific rerun, which requires the original data and substantial compute time.
