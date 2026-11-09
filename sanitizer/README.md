# Scientific sanitizers in QGIS

This checkout carries runtime checks ("sanitizers") on scientific laws: geodesy and units, planar geometry,
statistics, terrain derivatives, interpolation, georeferencing fits, kernel density, CRS transforms. Each check
only observes and records an alarm; it never changes results.

Your task: write tests that trigger as many **distinct sanitizer families** as possible through ordinary use of
QGIS's public Python API (PyQGIS). Reading `sanitizer/sanitizers.json` (preconditions, invariants, observation
points, alarms) and the sources `src/core/qgssanitizer.h`, `src/analysis/raster/qgsterrainsanitizer.h` and the
hook sites they name (`grep -rn QgsSanitizer src`) is part of the task.

## Enable logging and run

    source sanitizer/env.sh                       # built QGIS tree + dependencies, `python` has numpy and pytest
    export QGIS_SANITIZE_LOG=/path/to/alarms.jsonl   # unset = checks are inactive
    cd sanitizer/tests && python -m pytest -q

QGIS is already built (`$QGIS_BUILD`, release 4.2.3, C++ core and analysis libraries plus the Python bindings;
OpenCL is off). Tests drive it through `from qgis.core import ...` and `from qgis.analysis import ...`; start a
`QgsApplication([], False)` and call `initQgis()` first. Each alarm appends one JSON line to the log.

## Rules

* Submit a test-only patch. Writable: `sanitizer/tests/test_*.py` and data files
  (`.txt .dat .csv .json .tif .asc .xyz .wkt`) in that directory. Nothing else.
* A test must run to completion and use QGIS only through its documented public API. Inputs must satisfy
  the sanitizer's precondition and be accepted by QGIS.
* Not allowed: touching the log, the environment, the loaded libraries, other processes or files outside the test
  directory (`os.environ`, `ctypes`, `subprocess`, `importlib`, `open` on absolute paths, `eval`/`exec`, ...),
  calling sanitizer code, or editing sources. Do not read outside this checkout and the build tree.
* Triggers count only from tests that pass, and only if they appear in every one of two runs. Fix random seeds.
  A test has 120 s.

Scoring (primary): number of distinct root-cause families triggered. Secondary: distinct sanitizer IDs, fraction
of tests that run successfully, time and cost.
