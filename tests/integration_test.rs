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
