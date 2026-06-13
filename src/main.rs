use clap::Parser;
use num_complex::Complex64;
use rayon::prelude::*;
use tart_correlator::correlator;
use tart_correlator::observation::Observation;
use tart_correlator::visibility;

/// TART software correlator — load and inspect radio-observation HDF5 files.
#[derive(Parser)]
#[command(name = "tart-correlator", version)]
struct Cli {
    /// Path to an HDF5 observation file.
    #[arg(short, long = "data")]
    data: Option<String>,

    /// Print the mean of each antenna signal.
    #[arg(short, long)]
    means: bool,

    /// Print information about the observation.
    #[arg(short, long)]
    info: bool,

    /// Convert to complex baseband and print I/Q stats for each antenna.
    #[arg(short, long)]
    baseband: bool,

    /// Channel width in Hz for polyphase filterbank (requires --baseband).
    /// Default: single channel (full bandwidth).
    #[arg(long = "channel-width")]
    channel_width: Option<f64>,

    /// Correlate all antennas and print visibilities (requires --baseband).
    #[arg(short, long)]
    correlate: bool,

    /// Integration time in seconds for correlation (requires --correlate).
    /// Default: use all available samples.
    #[arg(long = "integration-time", default_value = "0.0")]
    integration_time: f64,

    /// Path to antenna positions JSON file (required for --save-vis).
    #[arg(long = "antenna-positions")]
    antenna_positions: Option<String>,

    /// Save visibilities to an HDF5 file (requires --correlate, --antenna-positions).
    #[arg(long = "save-vis")]
    save_vis: Option<String>,
}

fn main() {
    let cli = Cli::parse();

    if cli.means && cli.data.is_none() {
        eprintln!("Error: --means requires --data <file.h5>");
        std::process::exit(1);
    }
    if cli.info && cli.data.is_none() {
        eprintln!("Error: --info requires --data <file.h5>");
        std::process::exit(1);
    }
    if cli.baseband && cli.data.is_none() {
        eprintln!("Error: --baseband requires --data <file.h5>");
        std::process::exit(1);
    }
    if cli.channel_width.is_some() && !cli.baseband {
        eprintln!("Error: --channel-width requires --baseband");
        std::process::exit(1);
    }
    if cli.correlate && !cli.baseband {
        eprintln!("Error: --correlate requires --baseband");
        std::process::exit(1);
    }
    if cli.integration_time != 0.0 && !cli.correlate {
        eprintln!("Error: --integration-time requires --correlate");
        std::process::exit(1);
    }
    if cli.save_vis.is_some() && !cli.correlate {
        eprintln!("Error: --save-vis requires --correlate");
        std::process::exit(1);
    }
    if cli.save_vis.is_some() && cli.antenna_positions.is_none() {
        eprintln!("Error: --save-vis requires --antenna-positions <file.json>");
        std::process::exit(1);
    }

    if cli.data.is_none() {
        println!("TART software correlator v{}", env!("CARGO_PKG_VERSION"));
        println!("Use --help for usage.");
        return;
    }

    let path = cli.data.as_ref().unwrap();
    let obs = match Observation::from_hdf5(path) {
        Ok(o) => o,
        Err(e) => {
            eprintln!("Failed to load observation from `{path}`: {e}");
            std::process::exit(1);
        }
    };

    if cli.info {
        println!("{}", obs.info_string());
    }

    if cli.means {
        let means = obs.means();
        println!("Antenna means (bipolar -1..+1):");
        for (i, mean) in means.iter().enumerate() {
            println!("  antenna {:3}: {:+.6}", i, mean);
        }
    }

    if cli.baseband {
        use std::time::Instant;
        let t0 = Instant::now();
        let bb = obs.to_baseband();
        let elapsed = t0.elapsed();
        println!(
            "Complex baseband: {} antennas × {} samples, rate={:.3} MHz",
            bb.num_antenna(),
            bb.num_samples(),
            bb.sample_rate / 1e6
        );
        println!("Conversion took {:?}", elapsed);

        if let Some(cw) = cli.channel_width {
            // PFB channelize
            let t0 = Instant::now();
            let channelized = bb.channelize(cw);
            let ch_elapsed = t0.elapsed();

            let num_ch = channelized[0].num_channels;
            let actual_width = channelized[0].channel_width;
            let n_time = channelized[0].num_time_steps();

            println!(
                "\nPFB channelizer: requested width {:.3} kHz → {} channels × {:.3} kHz, {} taps, {} time steps per channel",
                cw / 1e3,
                num_ch,
                actual_width / 1e3,
                channelized[0].taps,
                n_time,
            );
            println!("PFB took {:?}", ch_elapsed);

            // If correlating, correlate ALL channels, ALL integrations
            if cli.correlate {
                let t0 = Instant::now();
                let int_time = if cli.integration_time > 0.0 {
                    cli.integration_time
                } else {
                    n_time as f64 / actual_width
                };

                // Correlate each channel independently (parallel over channels)
                // Returns [channel][integration][baseline]
                let all_vis: Vec<Vec<Vec<Complex64>>> = (0..num_ch)
                    .into_par_iter()
                    .map(|ch_idx| {
                        let ch_data: Vec<Vec<_>> = channelized
                            .iter()
                            .map(|ant| ant.channels[ch_idx].clone())
                            .collect();
                        let multi = correlator::correlate_channel_multi(&ch_data, actual_width, int_time);
                        multi.iter()
                            .map(|vis| vis.iter().map(|v| v.value).collect())
                            .collect()
                    })
                    .collect();
                let corr_elapsed = t0.elapsed();

                let n_int = all_vis.first().map(|c| c.len()).unwrap_or(0);
                let n_bl = all_vis.first().and_then(|c| c.first().map(|t| t.len())).unwrap_or(0);

                println!("\nCorrelation ({num_ch} channels × {n_int} integrations, {:.3} s each):", int_time);
                println!("Correlation took {:?}", corr_elapsed);
                println!("Baselines per integration: {n_bl}");

                // Save to HDF5 if requested
                if let Some(ref save_path) = cli.save_vis {
                    let ant_pos_path = cli.antenna_positions.as_ref().unwrap();
                    let ant_pos = visibility::load_antenna_positions(ant_pos_path)
                        .unwrap_or_else(|e| {
                            eprintln!("Failed to load antenna positions: {e}");
                            std::process::exit(1);
                        });
                    let ts = obs.timestamp.to_rfc3339();

                    // Baseline pairs from first channel, first integration
                    let ch0_data: Vec<Vec<_>> = channelized
                        .iter()
                        .map(|ant| ant.channels[0].clone())
                        .collect();
                    let multi_ref = correlator::correlate_channel_multi(&ch0_data, actual_width, int_time);
                    let bl_pairs: Vec<_> = multi_ref[0].iter().map(|v| (v.i, v.j)).collect();

                    if let Err(e) = visibility::write_visibilities_hdf5(
                        save_path, &obs.config, &ts, &bl_pairs, &all_vis, actual_width, &ant_pos,
                    ) {
                        eprintln!("Failed to save visibilities: {e}");
                    } else {
                        println!("Visibilities saved to {save_path}");
                        print_h5_summary(save_path);
                    }
                } else {
                    // Print per-channel amplitude summary (first integration only)
                    println!("\nPer-channel amplitude summary (integration 0):");
                    println!("{:>5}  {:>10}  {:>10}  {:>10}", "ch", "mean|V|", "max|V|", "min|V|");
                    for (ch_idx, ch_vis) in all_vis.iter().enumerate() {
                        let amps: Vec<f64> = ch_vis[0].iter().map(|v| v.norm()).collect();
                        let mean: f64 = amps.iter().sum::<f64>() / amps.len() as f64;
                        let max = amps.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                        let min = amps.iter().cloned().fold(f64::INFINITY, f64::min);
                        println!("  {:>3}  {:10.6}  {:10.6}  {:10.6}", ch_idx, mean, max, min);
                    }
                }
            } else {
                // Print per-channel power for the first antenna
                println!("\nPer-channel power (antenna 0, integrated over all time steps, dBFS):");
                let ant0 = &channelized[0];
                let mut powers: Vec<(usize, f64)> = ant0
                    .channels
                    .iter()
                    .enumerate()
                    .map(|(ch, data)| {
                        let power: f64 =
                            data.iter().map(|c| c.re * c.re + c.im * c.im).sum::<f64>()
                                / data.len() as f64;
                        (ch, power)
                    })
                    .collect();

                powers.sort_by_key(|(ch, _)| *ch);

                let max_power = powers.iter().map(|(_, p)| *p).fold(0.0f64, f64::max);
                for (ch, power) in &powers {
                    let db = if *power > 0.0 {
                        10.0 * power.log10()
                    } else {
                        f64::NEG_INFINITY
                    };
                    let db_rel = if max_power > 0.0 {
                        10.0 * (*power / max_power).log10()
                    } else {
                        0.0
                    };
                    let freq_mhz = *ch as f64 * actual_width / 1e6;
                    println!(
                        "  ch {:4}  {:8.3} MHz:  power={:+.2} dB  rel={:+.2} dB",
                        ch, freq_mhz, db, db_rel
                    );
                }
            }
        } else if cli.correlate {
            // Single-channel correlation with multi-integration support
            let t0 = Instant::now();
            let int_time = if cli.integration_time > 0.0 {
                cli.integration_time
            } else {
                bb.num_samples() as f64 / bb.sample_rate
            };
            let multi = correlator::correlate_channel_multi(&bb.data, bb.sample_rate, int_time);
            let corr_elapsed = t0.elapsed();

            let n_int = multi.len();
            let n_bl = multi.first().map(|v| v.len()).unwrap_or(0);
            println!("\nCorrelation (single channel × {n_int} integrations, {:.3} s each):", int_time);
            println!("Correlation took {:?}", corr_elapsed);
            println!("Baselines per integration: {n_bl}");

            // Save to HDF5 if requested
            if let Some(ref save_path) = cli.save_vis {
                save_vis_multi(save_path, &obs, &multi, bb.sample_rate, cli.antenna_positions.as_ref().unwrap());
                print_h5_summary(save_path);
            } else {
                println!("\nVisibilities (amplitude, phase) — integration 0:");
                for v in &multi[0] {
                    let amp = v.value.norm();
                    let phase = v.value.arg();
                    let vv_re = correlator::van_vleck_correction(v.value.re);
                    let vv_im = correlator::van_vleck_correction(v.value.im);
                    let vv_amp = (vv_re * vv_re + vv_im * vv_im).sqrt();
                    let vv_phase = vv_im.atan2(vv_re);
                    println!(
                        "  baseline ({:2},{:2}):  amp={:.6}  phase={:+.4} rad   (VV-corrected: amp={:.6}  phase={:+.4} rad)",
                        v.i, v.j, amp, phase, vv_amp, vv_phase,
                    );
                }
            }
        } else {
            // Single-channel stats
            println!("\nI/Q statistics per antenna (mean_I, mean_Q, rms):");
            for (i, ant) in bb.data.iter().enumerate() {
                let n = ant.len() as f64;
                let sum_i: f64 = ant.iter().map(|c| c.re).sum();
                let sum_q: f64 = ant.iter().map(|c| c.im).sum();
                let sum_sq: f64 = ant.iter().map(|c| c.re * c.re + c.im * c.im).sum();
                let mean_i = sum_i / n;
                let mean_q = sum_q / n;
                let rms = (sum_sq / n).sqrt();
                println!(
                    "  antenna {:3}: I={:+.6}  Q={:+.6}  rms={:.6}",
                    i, mean_i, mean_q, rms
                );
            }
        }
    }

    if !cli.info && !cli.means && !cli.baseband {
        // Default: just print a one-line summary
        println!(
            "Loaded observation: {} antennas, {} samples each, timestamp {}",
            obs.num_antenna(),
            obs.num_samples(),
            obs.timestamp
        );
    }
}

/// Helper: save multi-integration visibilities to HDF5.
fn save_vis_multi(
    path: &str,
    obs: &Observation,
    multi: &[Vec<tart_correlator::correlator::Visibility>],
    channel_width_hz: f64,
    ant_pos_path: &str,
) {
    let ant_pos = visibility::load_antenna_positions(ant_pos_path)
        .unwrap_or_else(|e| {
            eprintln!("Failed to load antenna positions: {e}");
            std::process::exit(1);
        });
    let ts = obs.timestamp.to_rfc3339();

    let bl_pairs: Vec<_> = multi[0].iter().map(|v| (v.i, v.j)).collect();

    // Build 3D: [1 channel][N_int integrations][N_bl baselines]
    let vis_3d: Vec<Vec<Vec<Complex64>>> = vec![multi
        .iter()
        .map(|vis| vis.iter().map(|v| v.value).collect())
        .collect()];

    if let Err(e) = visibility::write_visibilities_hdf5(
        path, &obs.config, &ts, &bl_pairs, &vis_3d, channel_width_hz, &ant_pos,
    ) {
        eprintln!("Failed to save visibilities: {e}");
    } else {
        println!("Visibilities saved to {path}");
    }
}

/// Print a summary of all datasets in an HDF5 file.
fn print_h5_summary(path: &str) {
    let h5 = match hdf5_reader::Hdf5File::open(path) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("Could not read back HDF5 file: {e}");
            return;
        }
    };
    println!("\nHDF5 file summary:");
    for name in &[
        "config", "timestamp", "phase_elaz", "baselines",
        "uvw", "antenna_positions", "gains", "phases",
        "chan_freq", "chan_width", "vis",
    ] {
        if let Ok(ds) = h5.dataset(name) {
            let shape: Vec<String> = ds.shape().iter().map(|s| s.to_string()).collect();
            let elem_size = hdf5_reader::dtype_element_size(ds.dtype()).unwrap_or(0);
            let n_elems: u64 = ds.shape().iter().product();
            let total_bytes = n_elems * elem_size as u64;
            println!(
                "  {:<20}  shape=({})  size={}",
                name,
                shape.join(", "),
                format_bytes(total_bytes),
            );
        }
    }
}

fn format_bytes(bytes: u64) -> String {
    if bytes >= 1_048_576 {
        format!("{:.1} MiB", bytes as f64 / 1_048_576.0)
    } else if bytes >= 1024 {
        format!("{:.1} KiB", bytes as f64 / 1024.0)
    } else {
        format!("{bytes} B")
    }
}
