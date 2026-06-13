use num_complex::Complex64;
use rayon::prelude::*;

/// A single complex visibility for one baseline.
#[derive(Debug, Clone)]
pub struct Visibility {
    /// Index of first antenna.
    pub i: usize,
    /// Index of second antenna.
    pub j: usize,
    /// Complex visibility: V = ⟨x_i · conj(x_j)⟩
    pub value: Complex64,
    /// Number of samples integrated.
    pub n_samples: usize,
}

/// Generate all baseline pairs (i, j) for N antennas where i < j.
/// Returns N(N-1)/2 pairs.
pub fn baselines(num_antennas: usize) -> Vec<(usize, usize)> {
    let n_bl = num_antennas * (num_antennas - 1) / 2;
    let mut bl = Vec::with_capacity(n_bl);
    for i in 0..num_antennas {
        for j in i + 1..num_antennas {
            bl.push((i, j));
        }
    }
    bl
}

/// Correlate all antennas over multiple sequential integration windows.
///
/// Divides the data into integration windows of `integration_time_s` seconds
/// each and returns one `Vec<Visibility>` per window.
///
/// Returns `Vec<Vec<Visibility>>` where the outer Vec is per-integration
/// (length N_int) and the inner Vec is per-baseline (length N_bl).
pub fn correlate_channel_multi(
    antenna_data: &[Vec<Complex64>],
    sample_rate: f64,
    integration_time_s: f64,
) -> Vec<Vec<Visibility>> {
    let num_antennas = antenna_data.len();
    if num_antennas < 2 {
        return Vec::new();
    }

    let n_total = antenna_data[0].len();
    let n_per = (sample_rate * integration_time_s).round() as usize;
    let n_per = n_per.min(n_total).max(1);
    let n_int = n_total / n_per;

    let bl_pairs = baselines(num_antennas);

    (0..n_int)
        .into_par_iter()
        .map(|i_win| {
            let start = i_win * n_per;
            bl_pairs
                .iter()
                .map(|&(i, j)| {
                    let xi = &antenna_data[i][start..][..n_per];
                    let xj = &antenna_data[j][start..][..n_per];

                    let mut sum = Complex64::new(0.0, 0.0);
                    for k in 0..n_per {
                        sum += xi[k] * xj[k].conj();
                    }
                    let value = sum / n_per as f64;

                    Visibility {
                        i,
                        j,
                        value,
                        n_samples: n_per,
                    }
                })
                .collect()
        })
        .collect()
}

/// Apply van Vleck correction for 1-bit quantization.
///
/// The 1-bit correlator output R is corrected to the true correlation ρ via:
///   ρ = sin(π/2 · R)
///
/// Reference: https://arxiv.org/abs/1608.04367
pub fn van_vleck_correction(r: f64) -> f64 {
    (std::f64::consts::PI / 2.0 * r).sin()
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;

    #[test]
    fn test_baselines() {
        let bl = baselines(4);
        assert_eq!(bl.len(), 6);
        assert_eq!(bl, vec![
            (0, 1), (0, 2), (0, 3),
            (1, 2), (1, 3),
            (2, 3),
        ]);
    }

    #[test]
    fn test_correlate_channel_basic() {
        // Two antennas with identical signals → correlation should be real, positive
        let n = 1000;
        let data = vec![
            (0..n).map(|i| Complex64::new((i as f64).sin(), 0.0)).collect::<Vec<_>>(),
            (0..n).map(|i| Complex64::new((i as f64).sin(), 0.0)).collect::<Vec<_>>(),
        ];

        let multi = correlate_channel_multi(&data, 100.0, 10.0);
        assert_eq!(multi.len(), 1);
        let vis = &multi[0];
        assert_eq!(vis.len(), 1);
        assert!(vis[0].value.re > 0.0);
        assert!(vis[0].value.im.abs() < 0.1);
    }

    #[test]
    fn test_correlate_channel_multi_windows() {
        // 2000 samples at 100 Hz = 20s; 5s integration = 4 windows
        let n = 2000;
        let data = vec![
            (0..n).map(|i| Complex64::new((i as f64).sin(), 0.0)).collect::<Vec<_>>(),
            (0..n).map(|i| Complex64::new((i as f64).sin(), 0.0)).collect::<Vec<_>>(),
        ];

        let multi = correlate_channel_multi(&data, 100.0, 5.0);
        assert_eq!(multi.len(), 4, "should have 4 integration windows");
        for vis in &multi {
            assert_eq!(vis.len(), 1);
            assert!(vis[0].value.re > 0.0);
        }
    }

    #[test]
    fn test_correlate_channel_quadrature() {
        let n = 4000;
        let data = vec![
            (0..n).map(|i| {
                let phase = 2.0 * std::f64::consts::PI * i as f64 / 100.0;
                Complex64::new(phase.cos(), phase.sin())
            }).collect::<Vec<_>>(),
            (0..n).map(|i| {
                let phase = 2.0 * std::f64::consts::PI * i as f64 / 100.0 + std::f64::consts::PI / 2.0;
                Complex64::new(phase.cos(), phase.sin())
            }).collect::<Vec<_>>(),
        ];

        let multi = correlate_channel_multi(&data, 100.0, 40.0);
        assert_eq!(multi.len(), 1);
        let vis = &multi[0];
        assert_eq!(vis.len(), 1);
        assert!(vis[0].value.im < -0.9);
        assert!(vis[0].value.re.abs() < 0.1);
    }

    #[test]
    fn test_van_vleck() {
        assert!((van_vleck_correction(0.0) - 0.0).abs() < 1e-10);
        assert!((van_vleck_correction(1.0) - 1.0).abs() < 1e-10);
        assert!((van_vleck_correction(-1.0) + 1.0).abs() < 1e-10);
    }
}
