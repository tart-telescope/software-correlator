# Changelog

## [0.4.0] — 2026-06-13

### Added
- Rayon parallelism: baseband conversion, PFB channelizer, correlator, and
  means computation now use all available CPU cores.
  - `to_baseband()`: 15× speedup in release mode (30s → 2.0s).
  - `channelize()`: 11× speedup (629ms → 56ms).
  - `correlate_channel()`: 17× speedup (287ms → 17ms).
  - Integration tests: 7× faster in debug mode (110s → 15s).
- `--info` now prints observation duration in milliseconds.

### Dependencies
- Added `rayon` for data-parallelism.

## [0.3.0] — 2026-06-13

### Added
- Antenna correlation (`--correlate`): computes complex visibilities for all
  baselines by direct complex correlation V_ij = ⟨x_i · conj(x_j)⟩.
- `--integration-time` argument (seconds) to control the correlation integration
  window.
- `Visibility` struct: carries baseline indices, complex value, and sample count.
- `baselines()` helper: generates all N(N-1)/2 antenna pairs.
- `van_vleck_correction()`: corrects 1-bit quantization bias via sin(π/2 · R).
- Unit tests for baselines, correlation (in-phase and quadrature), and van Vleck.
- Works in both single-channel and PFB-channelized modes.

## [0.2.0] — 2026-06-13

### Added
- Complex baseband conversion (`--baseband`): converts real 1-bit IF samples at
  16.368 MHz to complex I/Q baseband at 4.092 MHz via:
  - Unipolar (0,1) → bipolar (-1,+1) conversion
  - Complex mix-down by the IF center frequency
  - Anti-aliasing FIR filter (64-tap Hann-windowed sinc)
  - 4× decimation
- Polyphase filterbank channelizer (`--channel-width`): splits the complex
  baseband into equally-spaced frequency channels using a windowed-sinc
  prototype filter and an FFT-based PFB.
- `baseband_frequency` field added to `Settings` (was in JSON but not deserialized).

### Changed
- `lat`, `lon`, `alt` in `Settings` changed from `String` to `f64` to match the
  actual JSON format in HDF5 observation files.
- VLEN HDF5 dataset reading now uses `hdf5-reader` crate's built-in
  `resolve_vlen_reference_bytes()` API instead of manual global-heap parsing.
- Main binary uses library crate (`tart_correlator`) so modules are accessible
  from integration tests.

### Dependencies
- Added `num-complex` for complex number arithmetic.
- Added `rustfft` for FFT in the polyphase filterbank.
- Replaced `hdf5` with `hdf5-reader` (pure Rust, no C build dependency).

## [0.1.0] — 2026-06-12

### Added
- Initial release: load TART radio-astronomy observations from HDF5 files.
- `--data <file.h5>` argument to specify an observation file.
- `--info` argument to print observation metadata (timestamp, location,
  frequencies, antenna count, sample count).
- `--means` argument to print per-antenna signal means.
- `Observation` struct with HDF5 loading, bit-unpacking, Julian date, and MJD.
- `Settings` struct for telescope configuration deserialization.
- Integration test against `test-data/mg-tana-raw.hdf`.
