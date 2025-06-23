
# GPS Anti-Spoofing via Spatial Subspace Processing

This MATLAB project implements a GPS anti-spoofing algorithm based on spatial subspace estimation. The method uses a 2×2 antenna array to detect and mitigate spoofed GPS signals in complex environments where spoofed and authentic signals are both near the noise floor.

##  Overview

Traditional GPS spoofing detection relies on time-frequency domain correlation to identify PRN codes and Doppler shifts, which is computationally expensive. This project bypasses that need by estimating the **spoofing signal subspace vector (SSV)** using spatial information. It relies on the fact that spoofers transmit multiple PRNs from a single physical location, resulting in **constructive spatial energy accumulation**. In contrast, authentic signals arrive from multiple directions with lower aggregate spatial energy.

This technique:

* Avoids full 2D time-frequency acquisition for spoofed signals.
* Exploits spatial dominance of spoofing sources.
* Leverages PRN code periodicity for signal alignment and enhancement.



## 📁 Project Structure

```
.
├── array sig gen/                          # Antenna array signal generation utilities
├── Library/                                # Supporting functions and definitions
├── AS_Core1_withoutPMax.m                  # Core 1: Null Steering only
├── AS_core2_with_pmax.m                    # Core 2: Null Steering + Power Maximization
├── AS_results.png                          # Summary of acquisition results (pre/post AS)
├── main_GPS_AntiSpoofing_arraybased.m      # Main simulation script
├── plot_all_data.m                         # Visualization script
├── simulate_gps_signal_noActualData.m      # GNSS signal simulator (no real data)
└── README.md                               # This file

```


## 🚀 Processing Pipeline

1. **Signal Simulation**

   * GPS signal with spoofed PRNs and jamming is generated over 3 segments.
2. **Subspace Estimation**

   * Spatial information of spoofed signals is extracted using SSV estimation.
3. **Anti-Spoofing Core 1**

   * Projects signals orthogonally to spoofed/jammed subspace.
4. **Anti-Spoofing Core 2**

   * Refines with Pmax filtering using detected Doppler/PRNs.
5. **Signal Acquisition**

   * GPS PRNs acquired on original, Core 1, and Core 2 outputs.
6. **Beam Pattern Analysis**

   * Visualizes spatial filtering effect via beamforming gain plots.
7. **DOA Estimation**

   * Estimates direction of arrival for spoofing, jamming, and filtered signals.


## Requirements

* **MATLAB** (R2020 or newer recommended)
* Signal Processing Toolbox (for plotting and basic DSP)

## 🧠 Theory Reference

This project is based on techniques described in:

> Borio, D., Closas, P., & Humphreys, T. (2012). *A new approach to GPS spoofing detection based on spatial processing*. ION GNSS+ 2012, Session B3, Nashville, TN.

Key insight: spoofed signals from a single source exhibit dominant spatial energy, enabling their isolation without full code/Doppler acquisition.


## 📃 License

This project is intended for academic and research use.

