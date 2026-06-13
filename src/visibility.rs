use crate::settings::Settings;
use num_complex::Complex64;

/// Compute UVW coordinates (meters) for each baseline given antenna ENU positions.
///
/// Phase center is zenith (el=90°, az=0°). UVW convention follows CASA MS v2:
/// w toward source, v north, u east. Baseline direction is from ANTENNA2 to
/// ANTENNA1 (pos[i] - pos[j]).
pub fn compute_uvw(
    baselines: &[(usize, usize)],
    antenna_positions: &[[f64; 3]],
) -> Vec<[f64; 3]> {
    baselines
        .iter()
        .map(|&(i, j)| {
            let pi = antenna_positions[i];
            let pj = antenna_positions[j];
            // UVW from j to i: pos_i - pos_j
            // Zenith phase center → u=east, v=north, w=up
            [pi[0] - pj[0], pi[1] - pj[1], pi[2] - pj[2]]
        })
        .collect()
}

/// Compute channel center frequencies (Hz) for N equally-spaced channels.
///
/// Channels span [0, sample_rate) with the first channel at channel_width/2.
pub fn channel_frequencies(num_channels: usize, channel_width_hz: f64) -> Vec<f64> {
    (0..num_channels)
        .map(|ch| (ch as f64 + 0.5) * channel_width_hz)
        .collect()
}

/// Write visibilities and metadata to a CASA-style HDF5 file.
///
/// `vis` is a 3D array: [N_channel, N_integration, N_baseline].
pub fn write_visibilities_hdf5(
    path: &str,
    config: &Settings,
    timestamp: &str,
    baselines: &[(usize, usize)],
    vis_3d: &[Vec<Vec<Complex64>>],
    channel_width_hz: f64,
    antenna_positions: &[[f64; 3]],
) -> Result<(), Box<dyn std::error::Error>> {
    let mut builder = hdf5_pure::FileBuilder::new();
    let n_ant = antenna_positions.len();
    let n_ch = vis_3d.len();
    let n_int = vis_3d.first().map(|c| c.len()).unwrap_or(0);
    let n_bl = vis_3d
        .first()
        .and_then(|c| c.first().map(|t| t.len()))
        .unwrap_or(0);

    // config — JSON string (as bytes)
    let config_json = serde_json::to_string(config)?;
    builder
        .create_dataset("config")
        .with_u8_data(config_json.as_bytes());

    // phase_elaz
    builder
        .create_dataset("phase_elaz")
        .with_f64_data(&[90.0, 0.0]);

    // baselines [n_bl, 2]
    let bl_flat: Vec<i64> = baselines
        .iter()
        .flat_map(|&(i, j)| [i as i64, j as i64])
        .collect();
    builder
        .create_dataset("baselines")
        .with_shape(&[n_bl as u64, 2])
        .with_i64_data(&bl_flat);

    // vis [n_ch, n_int, n_bl] complex32
    let mut vis_flat: Vec<(f32, f32)> =
        Vec::with_capacity(n_ch * n_int * n_bl);
    for ch in 0..n_ch {
        for t in 0..n_int {
            for bl in 0..n_bl {
                let v = vis_3d[ch][t][bl];
                vis_flat.push((v.re as f32, v.im as f32));
            }
        }
    }
    builder
        .create_dataset("vis")
        .with_shape(&[n_ch as u64, n_int as u64, n_bl as u64])
        .with_complex32_data(&vis_flat);

    // SPECTRAL_WINDOW / chan_freq [n_ch]
    let chan_freq = channel_frequencies(n_ch, channel_width_hz);
    builder
        .create_dataset("chan_freq")
        .with_f64_data(&chan_freq);

    // SPECTRAL_WINDOW / chan_width
    let chan_widths = vec![channel_width_hz; n_ch];
    builder
        .create_dataset("chan_width")
        .with_f64_data(&chan_widths);

    // UVW [n_bl, 3]
    let uvw = compute_uvw(baselines, antenna_positions);
    let uvw_flat: Vec<f64> = uvw.iter().flat_map(|v| [v[0], v[1], v[2]]).collect();
    builder
        .create_dataset("uvw")
        .with_shape(&[n_bl as u64, 3])
        .with_f64_data(&uvw_flat);

    // gains (unity)
    builder
        .create_dataset("gains")
        .with_f32_data(&vec![1.0f32; n_ant]);

    // phases (zero)
    builder
        .create_dataset("phases")
        .with_f32_data(&vec![0.0f32; n_ant]);

    // antenna_positions [n_ant, 3]
    let ant_pos_flat: Vec<f32> = antenna_positions
        .iter()
        .flat_map(|p| [p[0] as f32, p[1] as f32, p[2] as f32])
        .collect();
    builder
        .create_dataset("antenna_positions")
        .with_shape(&[n_ant as u64, 3])
        .with_f32_data(&ant_pos_flat);

    // timestamp
    builder
        .create_dataset("timestamp")
        .with_u8_data(timestamp.as_bytes());

    builder.write(path)?;
    Ok(())
}

/// Load antenna positions from a JSON file.
pub fn load_antenna_positions(path: &str) -> Result<Vec<[f64; 3]>, Box<dyn std::error::Error>> {
    let raw = std::fs::read_to_string(path)?;
    let parsed: serde_json::Value = serde_json::from_str(&raw)?;
    let positions = parsed["antenna_positions"]
        .as_array()
        .ok_or("missing 'antenna_positions' array")?;

    let mut result = Vec::with_capacity(positions.len());
    for entry in positions {
        let coords: Vec<f64> = entry
            .as_array()
            .ok_or("antenna position is not an array")?
            .iter()
            .map(|v| v.as_f64().unwrap_or(0.0))
            .collect();
        if coords.len() != 3 {
            return Err(format!("expected 3 coordinates, got {}", coords.len()).into());
        }
        result.push([coords[0], coords[1], coords[2]]);
    }

    Ok(result)
}
