use crate::settings::Settings;
use num_complex::Complex64;

/// Write visibilities and metadata to a TART-format HDF5 file.
pub fn write_visibilities_hdf5(
    path: &str,
    config: &Settings,
    timestamp: &str,
    baselines: &[(usize, usize)],
    vis_values: &[Complex64],
    antenna_positions: &[[f64; 3]],
) -> Result<(), Box<dyn std::error::Error>> {
    let mut builder = hdf5_pure::FileBuilder::new();
    let n_ant = antenna_positions.len();
    let n_bl = vis_values.len();

    // config — JSON string
    let config_json = serde_json::to_string(config)?;
    let config_bytes: Vec<u8> = config_json.into_bytes();
    builder
        .create_dataset("config")
        .with_u8_data(&config_bytes);

    // phase_elaz
    builder
        .create_dataset("phase_elaz")
        .with_f64_data(&[90.0, 0.0]);

    // baselines
    let bl_flat: Vec<i64> = baselines
        .iter()
        .flat_map(|&(i, j)| [i as i64, j as i64])
        .collect();
    builder
        .create_dataset("baselines")
        .with_shape(&[n_bl as u64, 2])
        .with_i64_data(&bl_flat);

    // vis (complex64 as (f32,f32) pairs)
    let vis_c32: Vec<(f32, f32)> = vis_values
        .iter()
        .map(|v| (v.re as f32, v.im as f32))
        .collect();
    builder
        .create_dataset("vis")
        .with_complex32_data(&vis_c32);

    // gains (unity)
    let gains = vec![1.0f32; n_ant];
    builder
        .create_dataset("gains")
        .with_f32_data(&gains);

    // phases (zero)
    let phases = vec![0.0f32; n_ant];
    builder
        .create_dataset("phases")
        .with_f32_data(&phases);

    // antenna_positions
    let ant_pos_flat: Vec<f32> = antenna_positions
        .iter()
        .flat_map(|p| [p[0] as f32, p[1] as f32, p[2] as f32])
        .collect();
    builder
        .create_dataset("antenna_positions")
        .with_shape(&[n_ant as u64, 3])
        .with_f32_data(&ant_pos_flat);

    // timestamp
    let ts_bytes = timestamp.as_bytes().to_vec();
    builder
        .create_dataset("timestamp")
        .with_u8_data(&ts_bytes);

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
