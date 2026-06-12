# TART Software Correlator

A Rust-based software correlator for the [TART radio telescope](https://github.com/tmolteno/TART). Loads raw observation data from HDF5 files, converts to complex baseband, channelizes with a polyphase filterbank, and correlates all antenna pairs to produce complex visibilities.

## Installation

Requires Rust 1.85+ (edition 2024). No C build dependencies — uses `hdf5-reader` (pure Rust) for HDF5 I/O.

```bash
cargo build --release
```

## Usage

All features operate on a TART observation file produced by `obs.to_hdf5()`.

```
tart-correlator [OPTIONS]
```

### Options

| Flag | Description |
|------|-------------|
| `-d`, `--data <FILE>` | Path to an HDF5 observation file |
| `-i`, `--info` | Print observation metadata (timestamp, location, frequencies, antenna/sample counts) |
| `-m`, `--means` | Print the mean of each antenna's signal (bipolar -1..+1) |
| `-b`, `--baseband` | Convert real IF samples to complex baseband I/Q |
| `--channel-width <HZ>` | Polyphase filterbank channel width in Hz (requires `--baseband`) |
| `-c`, `--correlate` | Correlate all antenna pairs and print visibilities (requires `--baseband`) |
| `--integration-time <SECS>` | Integration window for correlation in seconds (requires `--correlate`; default: all samples) |
| `-h`, `--help` | Print help |
| `-V`, `--version` | Print version |

### Examples

**Inspect an observation:**

```bash
tart-correlator --data observation.h5 --info
```

Output:
```
Timestamp:            2026-06-12 11:20:40 UTC
Julian Date:          2461203.972685
MJD:                  61203.472685
Telescope name:       Madagascar - Antananarivo
Location:             lat=-18.8971446, lon=47.5551677, alt=1280.0m
Operating frequency:  1575.420 MHz
Bandwidth:            2.500 MHz
Sampling frequency:   16.368 MHz
Number of antennas:   24
Samples per antenna:  1048576
```

**Antenna signal means:**

```bash
tart-correlator --data observation.h5 --means
```

**Convert to complex baseband (single channel):**

```bash
tart-correlator --data observation.h5 --baseband
```

Converts 1-bit real IF samples (16.368 MHz, centered at 4.092 MHz) to complex I/Q baseband at 4.092 MHz via:
1. Unipolar (0,1) → bipolar (-1,+1)
2. Complex mix-down by the IF center frequency
3. 64-tap Hann-windowed FIR anti-aliasing filter
4. 4× decimation

**Polyphase filterbank channelizer:**

```bash
tart-correlator --data observation.h5 --baseband --channel-width 100000
```

Splits the 4.092 MHz band into 64 channels of ~63.9 kHz each (next power-of-two above 4.092 MHz / 100 kHz). Uses a 4-tap windowed-sinc prototype filter and an FFT-based PFB.

```bash
tart-correlator --data observation.h5 --baseband --channel-width 10000
```

~512 channels of ~8 kHz each for high-resolution spectroscopy.

**Correlate all antenna pairs:**

```bash
# Full bandwidth, 10 ms integration
tart-correlator --data observation.h5 --baseband --correlate --integration-time 0.01

# Per-channel (channel 0 only), 1 s integration
tart-correlator --data observation.h5 --baseband --channel-width 100000 --correlate --integration-time 1.0
```

Computes complex visibilities V_ij = ⟨x_i · conj(x_j)⟩ for all 276 baselines (24 antennas). Output includes raw amplitude/phase and van Vleck-corrected values.

### Combined examples

```bash
# Full pipeline: info + baseband + channelizer + correlation
tart-correlator --data observation.h5 --info --baseband --channel-width 50000 --correlate --integration-time 0.1
```

## Signal processing pipeline

```
HDF5 file
  │
  ├─ unpack 1-bit packed data → unipolar (0,1) samples at 16.368 MHz
  │
  ├─ --info: metadata display
  ├─ --means: per-antenna DC offsets
  │
  └─ --baseband:
       │
       ├─ unipolar → bipolar (-1,+1)
       ├─ complex mix-down: exp(-j 2π fc n / fs), fc = 4.092 MHz
       ├─ 64-tap FIR low-pass filter (Hann-windowed sinc)
       └─ 4× decimate → 4.092 MHz complex I/Q
            │
            ├─ (default): I/Q statistics per antenna
            │
            ├─ --channel-width: PFB channelizer
            │     ├─ prototype filter: M×4-tap windowed sinc
            │     ├─ polyphase decomposition
            │     ├─ M-point FFT per time step
            │     └─ per-channel power spectrum or --correlate
            │
            └─ --correlate: complex visibilities
                  ├─ V_ij = ⟨x_i · conj(x_j)⟩ / N
                  ├─ van Vleck correction: ρ = sin(π/2 · R)
                  └─ amplitude + phase per baseline
```

## Library API

The crate can also be used as a library:

```rust
use tart_correlator::observation::Observation;
use tart_correlator::correlator;

let obs = Observation::from_hdf5("observation.h5")?;
println!("{}", obs.info_string());

let bb = obs.to_baseband();
let channelized = bb.channelize(100_000.0); // 100 kHz channels

// Correlate channel 0 with 1 second integration
let ch0: Vec<Vec<_>> = channelized.iter().map(|ant| ant.channels[0].clone()).collect();
let vis = correlator::correlate_channel(&ch0, channelized[0].channel_width, 1.0);
```

## Test data

The repository includes `test-data/mg-tana-raw.hdf` — a TART observation from the Madagascar-Antananarivo telescope (24 antennas, 1048576 samples at 16.368 MHz).

Run tests:
```bash
cargo test                          # all tests (unit + integration)
cargo test -p tart-correlator --lib # unit tests only (fast)
```

## License

See [LICENSE](LICENSE).
