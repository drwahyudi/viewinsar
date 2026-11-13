# Publishing viewinsar

Everything here is a one-time setup plus a short per-release checklist. Until
step 2 is done, `pip install viewinsar` will **not** work for anyone else — the
name is unclaimed, not reserved.

---

## 0. Before the first release

- [ ] Create the GitHub repo and push this tree. The URLs in `pyproject.toml`,
      `recipe/meta.yaml` and `README.md` assume
      `https://github.com/drwahyudi/viewinsar` — change them if you use a
      different account or name.
- [ ] Check the name is still free: <https://pypi.org/project/viewinsar/>
      should return 404. (It did when this package was written.)
- [ ] Decide the version. `0.1.0` is set in `src/viewinsar/__init__.py`, which
      `pyproject.toml` reads — that one file is the single source of truth.

---

## 1. Build and check locally

```bash
pip install build twine
python -m build                 # -> dist/viewinsar-0.1.0.tar.gz and ...-py3-none-any.whl
twine check dist/*
```

Then prove the wheel really installs and runs, in a throwaway environment:

```bash
python -m venv /tmp/vi && /tmp/vi/bin/pip install dist/viewinsar-0.1.0-py3-none-any.whl
/tmp/vi/bin/viewinsar --version
```

---

## 2. PyPI

### 2a. Claim the name on TestPyPI first

```bash
twine upload --repository testpypi dist/*
pip install --index-url https://test.pypi.org/simple/ \
            --extra-index-url https://pypi.org/simple/ viewinsar
```

The extra index is needed because numpy, scipy and friends are not on TestPyPI.

### 2b. Upload to the real PyPI

```bash
twine upload dist/*
```

You will need an account and an API token (`pypi-...`), stored in
`~/.pypirc` or passed as `TWINE_PASSWORD` with `TWINE_USERNAME=__token__`.

### 2c. Then switch to trusted publishing (recommended)

`.github/workflows/publish.yml` uses PyPI's OIDC trusted publishing, so no token
ever lives in the repo. To enable it, once:

1. PyPI → your project → *Manage* → *Publishing* → *Add a new publisher*
2. Owner `drwahyudi`, repository `viewinsar`, workflow `publish.yml`,
   environment `pypi`
3. Repeat on TestPyPI with environment `testpypi`
4. In GitHub → *Settings* → *Environments*, create `pypi` and `testpypi`

After that, a release is:

```bash
# bump __version__ in src/viewinsar/__init__.py, update CHANGELOG.md, commit
git tag v0.1.1
git push origin v0.1.1
```

The workflow builds, verifies the tag matches `__version__`, uploads to PyPI and
opens a GitHub release. Use *Run workflow* → `testpypi` for a dry run.

---

## 3. conda-forge

conda-forge builds **from the PyPI sdist**, so PyPI has to come first.

1. Get the checksum of the uploaded sdist:

   ```bash
   openssl dgst -sha256 dist/viewinsar-0.1.0.tar.gz
   ```

2. Paste it into `recipe/meta.yaml` in place of the zeroed `sha256`.

3. Fork <https://github.com/conda-forge/staged-recipes>, copy this recipe to
   `recipes/viewinsar/meta.yaml`, and open a PR. Note that in staged-recipes the
   `license_file` must be a file **inside** the recipe directory, so copy
   `LICENSE` next to `meta.yaml` and change the line to `license_file: LICENSE`.

4. Every run requirement is on conda-forge already. MintPy is optional (a
   `run_constrained` entry), as are rasterio and pyogrio. Keep
   `matplotlib-base` rather than `matplotlib` — conda-forge policy, it avoids
   pulling in Qt.

5. Once merged, a `viewinsar-feedstock` repo appears. Later releases are then
   just a PR to that feedstock bumping `version` + `sha256` (the bot usually
   opens it for you within a day of a PyPI upload).

Then `conda install -c conda-forge viewinsar` works.

### Building the conda package locally first

```bash
conda install -c conda-forge conda-build
conda build recipe/            # uses the PyPI url + sha256
conda install -c local viewinsar
```

To build from the working tree instead, comment out `url`/`sha256` in
`meta.yaml` and uncomment `path: ..`.

---

## 4. Per-release checklist

- [ ] `pytest` passes (`pytest tests/`)
- [ ] `ruff check src tests` is clean
- [ ] `__version__` bumped in `src/viewinsar/__init__.py`
- [ ] `CHANGELOG.md` updated
- [ ] `python -m build && twine check dist/*`
- [ ] tag `vX.Y.Z` pushed → the publish workflow does the rest
- [ ] conda-forge feedstock PR (or let the bot do it)

---

## Notes

- **Version in one place.** `pyproject.toml` reads
  `viewinsar.__version__` via `[tool.setuptools.dynamic]`. Do not add a second
  version string anywhere.
- **License.** The package uses MintPy where it is installed, and follows its
  licence, GPL-3.0-or-later. `LICENSE` is the unmodified
  GPLv3 text; add your own copyright line at the top of the source files if you
  want one.
- **MintPy is optional.** The required dependencies are only what viewinsar
  imports (numpy, scipy, matplotlib, h5py, pyproj, pyshp, pillow);
  `pip install "viewinsar[mintpy]"` adds MintPy, and then the map comes from
  view.py (`view.prep_slice`, `view.plot_slice`), `--timeseries` from
  `tsview.timeseriesViewer`, the readers from `readfile` and `utils.coordinate`
  - APIs that have been stable for a long time. Without it `compat.py` and
  `mapview.py` stand in. The `lean` CI job runs the whole suite without MintPy.
