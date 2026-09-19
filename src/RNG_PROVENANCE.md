# Manuscript random-number compatibility

`manuscript_rng.c` adapts NumPy's BSD-licensed binomial inversion/BTPE code and implements PCG64 and bounded integer sampling for the manuscript's fixed seed. The R runner compiles this helper locally; no Python executable, Python package, or saved random draws are required to run the analysis.

Source reference: NumPy v2.2.3 `numpy/random/src/distributions/distributions.c` and `numpy/random/src/pcg64/pcg64.h`, retrieved from the NumPy GitHub repository. Compatibility was checked against the available NumPy 2.3.5 runtime. See `NUMPY_LICENSE.txt` for the required notice. The implementation is adapted, not an unmodified NumPy distribution.

The initialized PCG64 state for `default_rng(20260829)` is:

- State: 226717667518859928713350147869998893023
- Increment: 76530840036586822202556446989362783643
- Buffered 32-bit value: absent

These generator initialization constants are derived from NumPy's `SeedSequence(20260829)` initialization. They contain no study data and are not bootstrap draws or fitted results. The fixed-seed API rejects other seeds. Each new R pointer starts this stream; integer and multinomial calls consume it continuously, including buffered 32-bit values.

For paired biomarker bootstrap ordering, the helper also reconstructs CPython 3.12's positive-integer `set(defaultdict)` iteration algorithm (incremental insertion, table growth, linear probing, perturbation, ascending slot iteration). Input keys come from the R reconstruction of pregnancy candidate insertion order; no reference identifier list is used by the analysis. This matches the original script's unsorted factorization order before woman-level bootstrap sampling.
