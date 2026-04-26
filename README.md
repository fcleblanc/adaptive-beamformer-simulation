# Adaptive Beamformer Simulation

## Overview

This project implements and evaluates an adaptive beamforming algorithm for a linear antenna array using MATLAB. The goal is to isolate a desired signal arriving from a known direction while suppressing strong interference and random noise. The beamformer is developed and tested across three scenarios of increasing complexity, demonstrating both fixed and adaptive behavior in static and dynamic signal environments.

---

## Key Features

- Linear antenna array modeling using steering vectors  
- Linearly Constrained Minimum Variance (LCMV) beamformer  
- Adaptive covariance estimation using exponential weighting  
- Suppression of interference signals  
- Tracking and mitigation of time-varying interference  
- Visualization of spatial response and beam patterns over time  

---

## Project Structure

The project is divided into three parts, each exploring a different beamforming scenario:

### Part 1 – Fixed Beamformer

Implements an LCMV beamformer with known steering vectors for the desired signal and interference sources.

- Static signal environment  
- Two fixed interferers  
- Spatial response visualization  
- Comparison of raw vs. beamformed signals  

---

### Part 2 – Adaptive Beamformer with Emerging Interference

Extends the beamformer to adapt to changes in the signal environment.

- Covariance matrix updated using an exponential forgetting factor  
- Introduction of a new interference source during the simulation  
- Analysis of beam pattern evolution over time  

---

### Part 3 – Adaptive Beamformer with Time-Varying Interference

Further extends the adaptive beamformer to handle continuously changing interference directions.

- Interferer angles vary over time  
- Beamformer dynamically tracks and suppresses moving interference  
- Evaluation of response magnitude in desired and interference directions  
- Visualization of performance over time  

---

## Results

The simulations demonstrate:

- Effective suppression of interference signals  
- Preservation of the desired signal direction  
- Successful adaptation to:
  - Newly introduced interferers  
  - Continuously moving interference sources  

Beam patterns and response plots illustrate how the array adjusts its spatial filtering over time.

---

## Notes

- This project focuses on simulation and analysis of beamforming algorithms  
- The implementation is intended to demonstrate core concepts in array processing and adaptive filtering
- The simulation expects an `ABF.mat` file containing audio signals stored as 1-dimensional vectors. Signal length should be sufficient for the chosen number of simulation samples. The dataset is not     included and can be replaced with user-provided signals.
