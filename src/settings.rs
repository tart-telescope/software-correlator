use serde::{Deserialize, Serialize};

/// Telescope config settings, deserialized from the JSON stored in the HDF5 `config` dataset.
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct Settings {
    pub num_antenna: u32,
    pub sampling_frequency: f64,
    #[serde(default)]
    pub frequency: f64,
    #[serde(default)]
    pub bandwidth: f64,
    #[serde(default)]
    pub baseband_frequency: f64,
    #[serde(default)]
    pub name: String,
    #[serde(default)]
    pub lat: f64,
    #[serde(default)]
    pub lon: f64,
    #[serde(default)]
    pub alt: f64,
    #[serde(default)]
    #[allow(dead_code)]
    pub array_orientation: f64,
}

#[allow(dead_code)]
impl Settings {
    pub fn get_num_antenna(&self) -> u32 {
        self.num_antenna
    }
    pub fn get_sampling_frequency(&self) -> f64 {
        self.sampling_frequency
    }
    pub fn get_operating_frequency(&self) -> f64 {
        self.frequency
    }
    pub fn get_bandwidth(&self) -> f64 {
        self.bandwidth
    }
    pub fn get_baseband_frequency(&self) -> f64 {
        self.baseband_frequency
    }
    pub fn get_name(&self) -> &str {
        &self.name
    }
    pub fn get_lat(&self) -> f64 {
        self.lat
    }
    pub fn get_lon(&self) -> f64 {
        self.lon
    }
    pub fn get_alt(&self) -> f64 {
        self.alt
    }
}
