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

    if !cli.info && !cli.means {
        // Default: just print a one-line summary
        println!(
            "Loaded observation: {} antennas, {} samples each, timestamp {}",
            obs.num_antenna(),
            obs.num_samples(),
            obs.timestamp
        );
    }
}
