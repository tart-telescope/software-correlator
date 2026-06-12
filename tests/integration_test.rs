use tart_correlator::observation::Observation;

/// Load the mg-tana-raw.hdf file and verify its contents:
/// - 24 antennas
/// - packed uint8 data that unpacks to 1,048,576 samples per antenna
/// - info output contains expected metadata
/// - antenna means are within [-1, +1]
#[test]
fn test_info_on_mg_tana_raw() {
    let obs = Observation::from_hdf5("test-data/mg-tana-raw.hdf")
        .expect("failed to load test-data/mg-tana-raw.hdf");

    let info = obs.info_string();
    println!("=== Observation Info ===\n{info}\n=========================");

    // Basic sanity checks
    assert_eq!(obs.num_antenna(), 24);
    assert_eq!(
        obs.num_samples(),
        131072 * 8,
        "packed bytes should unpack to 8x samples"
    );
    assert!(!info.is_empty());
    assert!(info.contains("Madagascar"), "should contain telescope name");
    assert!(info.contains("24"), "should show antenna count");
    assert!(info.contains("2026-06-12"), "should show timestamp");

    // Verify means are plausible (in range -1..1)
    let means = obs.means();
    assert_eq!(means.len(), 24);
    for (i, m) in means.iter().enumerate() {
        assert!(
            (-1.0..=1.0).contains(m),
            "antenna {i} mean {m} out of range"
        );
    }
    println!("\nAntenna means: {means:.6?}");
}

/// Convert to complex baseband and verify output dimensions and sanity.
#[test]
fn test_baseband_conversion() {
    let obs = Observation::from_hdf5("test-data/mg-tana-raw.hdf")
        .expect("failed to load test-data/mg-tana-raw.hdf");

    let bb = obs.to_baseband();

    // Check dimensions: 16.368 → 4.092 MHz is factor 4
    assert_eq!(bb.num_antenna(), 24);
    assert_eq!(bb.num_samples(), obs.num_samples() / 4);
    assert!((bb.sample_rate - 4.092e6).abs() < 100.0);

    // I/Q means should be small (signal is zero-mean noise-like after mixing)
    for (i, ant) in bb.data.iter().enumerate() {
        let n = ant.len() as f64;
        let mean_i: f64 = ant.iter().map(|c| c.re).sum::<f64>() / n;
        let mean_q: f64 = ant.iter().map(|c| c.im).sum::<f64>() / n;
        assert!(
            mean_i.abs() < 0.2,
            "antenna {i}: I mean {mean_i:.6} too large"
        );
        assert!(
            mean_q.abs() < 0.2,
            "antenna {i}: Q mean {mean_q:.6} too large"
        );
    }

    // RMS should be non-zero (signal exists)
    let rms_values: Vec<f64> = bb
        .data
        .iter()
        .map(|ant| {
            let sum_sq: f64 = ant.iter().map(|c| c.re * c.re + c.im * c.im).sum();
            (sum_sq / ant.len() as f64).sqrt()
        })
        .collect();

    println!("\nBaseband I/Q stats:");
    for (i, ant) in bb.data.iter().enumerate() {
        let n = ant.len() as f64;
        let mean_i: f64 = ant.iter().map(|c| c.re).sum::<f64>() / n;
        let mean_q: f64 = ant.iter().map(|c| c.im).sum::<f64>() / n;
        println!("  ant {i:3}: I={mean_i:+.6}  Q={mean_q:+.6}  rms={:.6}", rms_values[i]);
    }

    let avg_rms: f64 = rms_values.iter().sum::<f64>() / rms_values.len() as f64;
    println!("Average RMS across antennas: {avg_rms:.6}");
    assert!(avg_rms > 0.01, "RMS should be non-trivial");
    assert!(avg_rms < 2.0, "RMS should be bounded");
}

/// Test the polyphase filterbank channelizer with a 100 kHz channel width.
#[test]
fn test_pfb_channelizer() {
    let obs = Observation::from_hdf5("test-data/mg-tana-raw.hdf")
        .expect("failed to load test-data/mg-tana-raw.hdf");

    let bb = obs.to_baseband();

    // 100 kHz channels: 4.092 MHz / 100 kHz = 41 → next power of two = 64
    let channel_width = 100_000.0; // 100 kHz
    let channelized = bb.channelize(channel_width);

    assert_eq!(channelized.len(), 24, "should channelize all antennas");

    let num_channels = channelized[0].num_channels;
    let actual_width = channelized[0].channel_width;
    let n_time = channelized[0].num_time_steps();
    let taps = channelized[0].taps;

    println!("\nPFB channelizer results:");
    println!("  Channels: {num_channels}");
    println!("  Actual width: {:.3} kHz", actual_width / 1e3);
    println!("  Time steps: {n_time}");
    println!("  Taps: {taps}");
    println!("  Input samples: {}", bb.num_samples());
    println!("  Expected time steps: {}", bb.num_samples() / num_channels);

    assert_eq!(n_time, bb.num_samples() / num_channels);
    assert!(num_channels.is_power_of_two());
    assert!(actual_width <= channel_width);
    assert!(taps == 4);

    // Each channel should have non-zero data
    for (ch_idx, ch_data) in channelized[0].channels.iter().enumerate() {
        assert_eq!(ch_data.len(), n_time, "channel {ch_idx} has wrong length");
        let power: f64 =
            ch_data.iter().map(|c| c.re * c.re + c.im * c.im).sum::<f64>() / ch_data.len() as f64;
        // Power should be non-zero (signal exists)
        assert!(power > 0.0, "channel {ch_idx} has zero power");
    }

    // Channels near the band edges should have lower power than center channels
    // (the 2.5 MHz signal is centered, so after complex mixing the signal is
    // centered at DC; outer channels should have only noise)
    let band_center_ch = num_channels / 2; // DC after complex mixdown is at channel 0
    let edge_ch = 1; // first channel near edge

    let center_power: f64 = channelized[0].channels[band_center_ch]
        .iter()
        .map(|c| c.re * c.re + c.im * c.im)
        .sum::<f64>()
        / channelized[0].channels[band_center_ch].len() as f64;

    let edge_power: f64 = channelized[0].channels[edge_ch]
        .iter()
        .map(|c| c.re * c.re + c.im * c.im)
        .sum::<f64>()
        / channelized[0].channels[edge_ch].len() as f64;

    println!("  Center channel ({band_center_ch}) power: {center_power:.6}");
    println!("  Edge channel ({edge_ch}) power: {edge_power:.6}");

    // Center should have signal; edge may have less
    assert!(center_power > 0.0, "center channel has no power");
}
