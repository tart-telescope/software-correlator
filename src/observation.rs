use chrono::{DateTime, Utc};
use num_complex::Complex64;

use crate::settings::Settings;

/// Complex baseband output from converting a real IF observation.
pub struct ComplexBaseband {
    /// Per-antenna complex baseband samples (I + jQ).
    pub data: Vec<Vec<Complex64>>,
    /// Sample rate in Hz of the baseband signal.
    pub sample_rate: f64,
    /// Baseband center frequency that was mixed down from (Hz).
    pub center_frequency: f64,
}

/// Output of the polyphase filterbank channelizer.
///
/// `channels[ch][t]` = complex sample in channel `ch` at time step `t`.
pub struct Channelized {
    /// Per-channel time series.  Shape: [num_channels][num_time_steps]
    pub channels: Vec<Vec<Complex64>>,
    /// Bandwidth of each channel in Hz.
    pub channel_width: f64,
    /// Number of channels.
    pub num_channels: usize,
    /// Number of PFB taps (prototype filter length = num_channels * taps).
    pub taps: usize,
}

impl Channelized {
    pub fn num_time_steps(&self) -> usize {
        self.channels.first().map(|c| c.len()).unwrap_or(0)
    }
}

impl ComplexBaseband {
    pub fn num_antenna(&self) -> usize {
        self.data.len()
    }

    pub fn num_samples(&self) -> usize {
        self.data.first().map(|a| a.len()).unwrap_or(0)
    }

    /// Channelize all antennas using a polyphase filterbank.
    ///
    /// `channel_width_hz` is the desired width of each channel in Hz.
    /// The actual number of channels is the next power-of-two such that
    /// the channel width does not exceed `channel_width_hz`.
    /// Uses 4 PFB taps.
    pub fn channelize(&self, channel_width_hz: f64) -> Vec<Channelized> {
        let taps = 4;
        let sample_rate = self.sample_rate;
        let raw_channels = (sample_rate / channel_width_hz).ceil() as usize;
        let num_channels = next_power_of_two(raw_channels).max(2);
        let actual_width = sample_rate / num_channels as f64;

        self.data
            .iter()
            .map(|ant| {
                let mut ch = Self::pfb_channelize(ant, num_channels, taps);
                ch.channel_width = actual_width;
                ch
            })
            .collect()
    }

    /// Apply a polyphase filterbank to channelize a single antenna's baseband.
    ///
    /// Splits the band into `num_channels` equally-spaced channels using a
    /// windowed-sinc prototype filter with `taps` polyphase taps per channel.
    ///
    /// The prototype filter length is `num_channels * taps`.
    ///
    /// Returns `num_channels` complex time series, one per channel.
    pub fn pfb_channelize(
        antenna_data: &[Complex64],
        num_channels: usize,
        taps: usize,
    ) -> Channelized {
        assert!(num_channels >= 2, "need at least 2 channels");
        assert!(num_channels.is_power_of_two(), "num_channels must be a power of two");

        let proto_len = num_channels * taps;
        let proto_filter = design_pfb_prototype(num_channels, taps);

        // Polyphase decomposition: sub-filters[i] has taps at indices i, i+M, i+2M, ...
        // The PFB processes M (num_channels) samples per output time step.
        // We reverse the filter coefficients and conjugate for the standard PFB.
        let n_in = antenna_data.len();
        let n_time = n_in / num_channels;

        let mut channels: Vec<Vec<Complex64>> = (0..num_channels)
            .map(|_| Vec::with_capacity(n_time))
            .collect();

        let mut planner = rustfft::FftPlanner::new();
        let fft = planner.plan_fft_forward(num_channels);

        let mut fft_buf = vec![Complex64::new(0.0, 0.0); num_channels];

        for t in 0..n_time {
            let base = t * num_channels;

            // For each polyphase branch (channel), convolve with its sub-filter
            for ch in 0..num_channels {
                let mut acc = Complex64::new(0.0, 0.0);
                for p in 0..taps {
                    let sample_idx = base as isize + ch as isize - p as isize * num_channels as isize;
                    if sample_idx >= 0 && (sample_idx as usize) < n_in {
                        let coeff_idx = p * num_channels + ch;
                        // Standard PFB: use conjugated, time-reversed prototype filter
                        let coeff = proto_filter[proto_len - 1 - coeff_idx].conj();
                        acc += antenna_data[sample_idx as usize] * coeff;
                    }
                }
                fft_buf[ch] = acc;
            }

            // FFT over the polyphase outputs → channelized spectrum
            fft.process(&mut fft_buf);

            for ch in 0..num_channels {
                channels[ch].push(fft_buf[ch]);
            }
        }

        Channelized {
            channels,
            channel_width: 0.0, // caller sets this
            num_channels,
            taps,
        }
    }
}

/// Represents one observation: timestamp, telescope config, and antenna signal data.
///
/// Each antenna's signal is stored as a `Vec<u8>` of unipolar samples (0 or 1).
/// The HDF5 file stores these as bit-packed `uint8` arrays which are unpacked on load.
pub struct Observation {
    pub timestamp: DateTime<Utc>,
    pub config: Settings,
    /// Per-antenna unipolar binary data (0 or 1 values, one `u8` per sample).
    pub data: Vec<Vec<u8>>,
}

impl Observation {
    /// Retrn the number of antennas.
    pub fn num_antenna(&self) -> usize {
        self.data.len()
    }

    /// Retrn the number of samples per antenna.
    pub fn num_samples(&self) -> usize {
        self.data.first().map(|a| a.len()).unwrap_or(0)
    }

    /// Calculate and retrn the mean of each antenna's signal.
    ///
    /// The mean is computed as `(sum(data[i]) / len(data[i])) * 2 - 1`,
    /// mapping unipolar 0→-1, 1→+1 before averaging, matching the Python `boolean_mean`
    /// convention.
    pub fn means(&self) -> Vec<f64> {
        self.data
            .iter()
            .map(|ant| {
                let sum: f64 = ant.iter().map(|&v| v as f64).sum();
                let n = ant.len() as f64;
                if n == 0.0 {
                    0.0
                } else {
                    (sum / n) * 2.0 - 1.0
                }
            })
            .collect()
    }

    /// Load an `Observation` from an HDF5 file.
    ///
    /// The file must contain:
    /// - `config`    — VLEN byte dataset (JSON string)
    /// - `timestamp` — VLEN byte dataset (ISO-8601 string)
    /// - `data`      — packed uint8 arrays, one per antenna (2D regular or 1D VLEN)
    pub fn from_hdf5(path: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let h5 = hdf5_reader::Hdf5File::open(path)?;

        // --- config ---
        let config_json = read_vlen_bytes_dataset(&h5, "config")?;
        let config: Settings = serde_json::from_str(&config_json)?;

        // --- timestamp ---
        let ts_str = read_vlen_bytes_dataset(&h5, "timestamp")?;
        let timestamp = DateTime::parse_from_rfc3339(&ts_str)
            .or_else(|_| DateTime::parse_from_str(&ts_str, "%Y-%m-%dT%H:%M:%S%.f%:z"))?;
        let timestamp = timestamp.with_timezone(&Utc);

        // --- data ---
        let data_ds = h5.dataset("data")?;
        let shape = data_ds.shape();
        let num_antennas: usize = shape[0].try_into()?;

        let mut data: Vec<Vec<u8>> = Vec::with_capacity(num_antennas);

        if shape.len() == 2 {
            // Regular 2D dataset — all rows have the same length
            let array = data_ds.read_array::<u8>()?;
            let row_len: usize = shape[1].try_into()?;
            let flat = array
                .as_slice_memory_order()
                .ok_or("data array not contiguous")?;
            for i in 0..num_antennas {
                let start = i * row_len;
                let end = start + row_len;
                let unpacked = unpack_bits(&flat[start..end]);
                data.push(unpacked);
            }
        } else {
            // 1D VLEN dataset — each element is a variable-length u8 sequence
            let raw = data_ds.read_raw_bytes()?;
            let ref_size = data_ds.vlen_reference_size();

            for i in 0..num_antennas {
                let off = i * ref_size;
                if off + ref_size > raw.len() {
                    return Err(format!(
                        "data VLEN reference at index {i} out of bounds: need {ref_size} bytes at offset {off}, have {}",
                        raw.len()
                    ).into());
                }
                let ref_bytes = &raw[off..off + ref_size];
                // base element size is 1 (u8)
                let packed = data_ds.resolve_vlen_reference_bytes(ref_bytes, 1)?;
                data.push(unpack_bits(&packed));
            }
        }

        Ok(Observation {
            timestamp,
            config,
            data,
        })
    }

    /// Retrn basic information about this observation as a multi-line string.
    pub fn info_string(&self) -> String {
        let mut lines = Vec::new();
        lines.push(format!("Timestamp:            {}", self.timestamp));
        lines.push(format!(
            "Julian Date:          {:.6}",
            self.julian_date()
        ));
        lines.push(format!("MJD:                  {:.6}", self.mjd()));
        lines.push(format!("Telescope name:       {}", &self.config.name));
        lines.push(format!(
            "Location:             lat={:.7}, lon={:.7}, alt={:.1}m",
            self.config.lat, self.config.lon, self.config.alt
        ));
        lines.push(format!(
            "Operating frequency:  {:.3} MHz",
            self.config.frequency / 1e6
        ));
        lines.push(format!(
            "Bandwidth:            {:.3} MHz",
            self.config.bandwidth / 1e6
        ));
        lines.push(format!(
            "Sampling frequency:   {:.3} MHz",
            self.config.sampling_frequency / 1e6
        ));
        lines.push(format!("Number of antennas:   {}", self.config.num_antenna));
        lines.push(format!(
            "Samples per antenna:  {}",
            self.num_samples()
        ));
        lines.join("\n")
    }

    /// Julian Date for the observation timestamp.
    pub fn julian_date(&self) -> f64 {
        // JD = (timestamp - Unix epoch).days + 2440587.5
        let unix_epoch = DateTime::from_timestamp(0, 0).unwrap();
        let duration = self.timestamp.signed_duration_since(unix_epoch);
        let days = duration.num_milliseconds() as f64 / 86_400_000.0;
        days + 2_440_587.5
    }

    /// Modified Julian Date for the observation timestamp.
    pub fn mjd(&self) -> f64 {
        self.julian_date() - 2_400_000.5
    }

    /// Convert from real IF samples to complex baseband.
    ///
    /// The observation contains 1-bit real samples at `sampling_frequency` Hz,
    /// centered at `baseband_frequency` Hz.  This method:
    ///
    /// 1. Converts unipolar (0,1) to bipolar (-1,+1)
    /// 2. Complex-mixes down by the center frequency: `exp(-j 2π fc n / fs)`
    /// 3. Applies a low-pass FIR filter to prevent aliasing
    /// 4. Decimates to the baseband rate
    ///
    /// The output sample rate equals `baseband_frequency`.
    pub fn to_baseband(&self) -> ComplexBaseband {
        let fs_in = self.config.sampling_frequency; // 16.368 MHz
        let fc = self.config.baseband_frequency; // 4.092 MHz
        let fs_out = fc; // decimate to 4.092 MHz

        let decim_factor = (fs_in / fs_out).round() as usize;
        assert!(decim_factor >= 1, "decimation factor must be >= 1");

        // Design a low-pass FIR filter with cutoff just below fs_out/2.
        // We use a 64-tap windowed-sinc filter.
        let filter = design_lowpass_fir(decim_factor, 64);

        let baseband_data: Vec<Vec<Complex64>> = self
            .data
            .iter()
            .map(|ant| {
                convert_one_antenna(ant, fc, fs_in, decim_factor, &filter)
            })
            .collect();

        ComplexBaseband {
            data: baseband_data,
            sample_rate: fs_out,
            center_frequency: fc,
        }
    }
}

/// Read a scalar VLEN byte dataset and retrn its contents as a String.
fn read_vlen_bytes_dataset(
    h5: &hdf5_reader::Hdf5File,
    name: &str,
) -> Result<String, Box<dyn std::error::Error>> {
    let ds = h5.dataset(name)?;
    let raw = ds.read_raw_bytes()?;
    let ref_size = ds.vlen_reference_size();
    let num_elements = ds.num_elements()? as usize;

    if raw.len() < ref_size * num_elements {
        return Err(format!(
            "dataset '{name}' too short for VLEN: have {}, need {} ({} refs × {} bytes)",
            raw.len(),
            ref_size * num_elements,
            num_elements,
            ref_size
        )
        .into());
    }

    // The base element size for a VLEN byte dataset (base=u8) is 1.
    let data = ds.resolve_vlen_reference_bytes(&raw[..ref_size], 1)?;
    Ok(String::from_utf8(data)?)
}

/// Unpack an MSB-first bit-packed `&[u8]` into a `Vec<u8>` of 0/1 values.
///
/// Each byte yields 8 bits, MSB first.  The caller is responsible for knowing
/// the true number of valid samples (the last byte may be partially padded in
/// the original — numpy's `packbits` always pads to a full byte).  We return
/// all 8*len bits.
fn unpack_bits(packed: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(packed.len() * 8);
    for &byte in packed {
        for shift in (0..8).rev() {
            out.push((byte >> shift) & 1);
        }
    }
    out
}

/// Design a low-pass FIR filter for anti-aliasing before decimation.
///
/// Returns `ntaps` coefficients for a windowed-sinc filter with cutoff
/// at `fs_out / 2 = fs_in / (2 * decim_factor)`, i.e. the Nyquist frequency
/// of the decimated output.  Uses a Hann window.
fn design_lowpass_fir(decim_factor: usize, ntaps: usize) -> Vec<f64> {
    let cutoff = 1.0 / (2.0 * decim_factor as f64); // normalized cutoff (Nyquist = 0.5)
    let half = (ntaps - 1) as f64 / 2.0;
    let mut taps = Vec::with_capacity(ntaps);

    for i in 0..ntaps {
        let n = i as f64 - half;
        if n.abs() < 1e-12 {
            // sinc(0) = 1
            taps.push(2.0 * cutoff);
        } else {
            let sinc = (2.0 * std::f64::consts::PI * cutoff * n).sin()
                / (std::f64::consts::PI * n);
            // Hann window
            let window = 0.5 * (1.0 - (2.0 * std::f64::consts::PI * i as f64 / (ntaps - 1) as f64).cos());
            taps.push(sinc * window);
        }
    }

    // Normalize to unity gain at DC
    let sum: f64 = taps.iter().sum();
    for t in &mut taps {
        *t /= sum;
    }

    taps
}

/// Design the prototype low-pass filter for a polyphase filterbank.
///
/// This is a windowed-sinc filter of length `num_channels * taps` with
/// cutoff at the channel width (i.e., 1/M in normalized frequency where
/// M = num_channels).  Uses a Hann window.
fn design_pfb_prototype(num_channels: usize, taps: usize) -> Vec<Complex64> {
    let proto_len = num_channels * taps;
    let mut coeffs = Vec::with_capacity(proto_len);
    let half = (proto_len - 1) as f64 / 2.0;

    // Cutoff at half the channel spacing: 1/(2*M) normalized
    let cutoff = 1.0 / (2.0 * num_channels as f64);

    for i in 0..proto_len {
        let n = i as f64 - half;
        let val = if n.abs() < 1e-12 {
            2.0 * cutoff
        } else {
            let sinc = (2.0 * std::f64::consts::PI * cutoff * n).sin()
                / (std::f64::consts::PI * n);
            let window = 0.5
                * (1.0
                    - (2.0 * std::f64::consts::PI * i as f64 / (proto_len - 1) as f64)
                        .cos());
            sinc * window
        };
        coeffs.push(Complex64::new(val, 0.0));
    }

    // Normalize
    let sum: f64 = coeffs.iter().map(|c| c.re).sum();
    for c in &mut coeffs {
        c.re /= sum;
    }

    coeffs
}

/// Round up to the next power of two.
fn next_power_of_two(n: usize) -> usize {
    if n == 0 {
        return 1;
    }
    let mut p = 1;
    while p < n {
        p <<= 1;
    }
    p
}

/// Convert a single antenna's real unipolar samples to complex baseband.
fn convert_one_antenna(
    unipolar: &[u8],
    fc: f64,
    fs: f64,
    decim_factor: usize,
    filter: &[f64],
) -> Vec<Complex64> {
    let n_in = unipolar.len();
    let n_out = n_in / decim_factor;
    let mut out = Vec::with_capacity(n_out);

    let ntaps = filter.len();
    let phase_per_sample = -2.0 * std::f64::consts::PI * fc / fs;

    for i_out in 0..n_out {
        let i_center = i_out * decim_factor;

        // Apply FIR filter around the decimated sample, then mix down
        let mut acc = Complex64::new(0.0, 0.0);
        let tap_start = if i_center >= ntaps / 2 {
            i_center - ntaps / 2
        } else {
            0
        };

        for (tap_idx, &coeff) in filter.iter().enumerate() {
            let sample_idx = tap_start + tap_idx;
            if sample_idx >= n_in {
                break;
            }
            // Convert unipolar (0,1) → bipolar (-1,+1)
            let bipolar = (unipolar[sample_idx] as f64) * 2.0 - 1.0;

            // Complex mix down: exp(j * phase_per_sample * sample_idx)
            let phase = phase_per_sample * sample_idx as f64;
            let mixer = Complex64::new(phase.cos(), phase.sin());

            acc += (bipolar * coeff) * mixer;
        }

        out.push(acc);
    }

    out
}
