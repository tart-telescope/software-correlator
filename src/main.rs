use clap::Parser;
use tart_correlator::observation::Observation;

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

            // Sort by channel index for display
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
