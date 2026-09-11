# Contributing

This repository is primarily a scientific-reproducibility record. Contributions that improve portability, documentation, tests, numerical diagnostics, or bug fixes are welcome.

## Issues

When reporting a bug, please include:

- the exact command used;
- Python and operating-system versions;
- the relevant traceback/log excerpt;
- whether the run used live catalog downloads or a frozen catalog cache;
- for stochastic calculations, the seed/configuration and realization count.

Do not upload raw DAQ files or other collaboration-controlled data to an issue.

## Pull requests

Please keep scientific-scope changes separate from software-maintenance changes. A PR that changes event selection, calibration, trial definitions, source families, or physics-response assumptions should explain the scientific rationale and expected numerical impact explicitly.

Before opening a PR, run:

```bash
python scripts/repository_smoke_check.py
```

The GitHub Actions workflow repeats the static checks on pushes and pull requests.
