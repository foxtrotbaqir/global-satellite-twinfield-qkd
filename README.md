# Satellite Twin-Field QKD Simulation

![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)

**Keywords:** TF-QKD, Twin-Field QKD, Satellite QKD, Quantum Communication, MATLAB, SKR, Geodesic, Path Search, Simulation, Distributed Trust

---

## Project Overview
This repository provides a **MATLAB simulation framework** for evaluating **global, distributed-trust satellite Twin-Field Quantum Key Distribution (TF-QKD)**.  

It models realistic **orbital dynamics**, **inter-satellite links (ISLs)**, and **ground-to-satellite uplinks**, including key physical impairments:

- Geometric losses  
- Atmospheric extinction (Beer-Lambert)  
- Telescope aperture sizing  
- Beam pointing errors  
- Synchronization penalties  
- Sun/Moon background noise  

A **graph-based routing engine** calculates shortest paths and k-disjoint paths between ground stations. Each path's **secret key rate (SKR)** is computed using an **asymmetric decoy-state, finite-key TF-QKD model**. End-to-end secure throughput is evaluated using **XOR-based key forwarding** across multiple simultaneous paths.

Simulations quantify **multi-path availability, SKR continuity, outages, delivered keys per pass, and global throughput** over month-long orbital cycles.

---

## Purpose
This repository is intended for **academic, research, and learning purposes only**.  
Commercial use is **prohibited**, and proper attribution must be given when using or referencing this work.

---

## Usage
1. Open the MATLAB project or script files.  
2. Run the main simulation script (`main_simulation.m` or equivalent).  
3. Outputs include:
   - CSV dataset of paths and SKR values: `tfqkd_paths_xor_dataset.csv`  
   - Bar plot of XOR SKR over time: `tfqkd_skr_barplot.png`

---

## License
This project is licensed under the **Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)** license.  

- You may **share and adapt** the material **for non-commercial purposes only**.  
- **Attribution is required**: credit the author when using, modifying, or referencing the work.  
- **Commercial use is prohibited**.  

[Full license text](https://creativecommons.org/licenses/by-nc/4.0/legalcode)

---

## References
- The Thesis: *“Globally Secure Satellite QKD with Distributed Trust”*  
- Relevant literature on TF-QKD, satellite constellations, and secure quantum communication.
