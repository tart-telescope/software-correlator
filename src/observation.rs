use chrono::{DateTime, Utc};

use crate::settings::Settings;

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
