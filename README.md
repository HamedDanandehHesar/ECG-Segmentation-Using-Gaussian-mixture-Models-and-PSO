<img width="1366" height="643" alt="untitled1" src="https://github.com/user-attachments/assets/d1fdb1f7-0eed-4de5-8eec-f406922520b7" />![<img width="1366" height="643" alt="untitled2" src="https://github.com/user-attachments/assets/8613ccbc-1f79-4465-9dec-52295f4cb220" /><img width="1366" height="643" alt="untitled" src="https://github.com/user-attachments/assets/783ee007-d3b7-456b-8fa9-fbb0177e3fc7" />
Uploading untitled1.svg…]()


:::
# ECG Wave Segmentation using Phase-Domain Gaussian Modeling

MATLAB implementation for **automatic segmentation of ECG waves (P, QRS, and T)** using a **phase‑domain Gaussian mixture model**.

The algorithm first learns a **parametric model of ECG morphology** from the mean ECG waveform and then uses the resulting **synthetic ECG components** to perform **robust segmentation of P, QRS, and T waves** across the entire signal.

This framework provides:

- Detection of **P peaks, T peaks, and R peaks**
- Estimation of **P onset / P offset**
- Estimation of **QRS onset / QRS offset**
- Estimation of **T onset / T offset**

The approach combines **phase-domain ECG alignment**, **Gaussian kernel modeling**, and **energy-based segmentation**.

---

# Overview

Electrocardiogram (ECG) signals consist of repeating cardiac cycles composed of:

- **P wave** – atrial depolarization  
- **QRS complex** – ventricular depolarization  
- **T wave** – ventricular repolarization  

Reliable detection of these waves is a fundamental problem in biomedical signal processing.

Traditional ECG segmentation methods operate directly in the **time domain**, which can be sensitive to:

- heart rate variability  
- morphological differences between beats  
- noise and baseline wander

This project addresses these issues by transforming the ECG signal into the **cardiac phase domain**, where each heartbeat is normalized to a phase interval:

```
[-π , π]
```

with

```
θ = 0  →  R peak
```

Using this representation, the algorithm extracts the **mean ECG morphology**, models it with **Gaussian kernels**, and generates **synthetic ECG components** that guide the segmentation of the original signal.

---

# Main Capabilities

The code automatically performs:

- R‑peak detection using **Pan–Tompkins algorithm**
- ECG **phase transformation**
- **Dynamic Time Warping (DTW)** alignment of heartbeats
- Extraction of **mean ECG morphology**
- Modeling of **P wave, QRS complex, and T wave** using Gaussian mixtures
- Generation of **synthetic ECG signals**
- Detection of **P peak and T peak**
- Detection of **P onset / P offset**
- Detection of **QRS onset / QRS offset**
- Detection of **T onset / T offset**

---

# Input Data

The code expects a MATLAB `.mat` file containing:

```
x   → ECG signal matrix
fs  → sampling frequency
```

Example loading:

```matlab
[file, path] = uigetfile('*.mat','Select ECG mat file');
data = load(file);

fs = data.fs;
x = data.x;

ecg = x(1,2000:8000);
```

Only the **first channel** of the signal is used.

---

# Algorithm Pipeline

The algorithm consists of several main stages.

---

# 1. ECG Preprocessing

Baseline wander is removed using median filtering:

```matlab
ecg = ecg - medfilt1(ecg,round(fs/3));
```

This improves the stability of the segmentation process.

---

# 2. R‑Peak Detection

R peaks are detected using the **Pan–Tompkins QRS detection algorithm**:

```matlab
[qrs_positions] = pantompkins_qrs(abs(ecg),fs);
```

These peaks define the boundaries of cardiac cycles.

---

# 3. ECG Phase Calculation

Each ECG sample is mapped to a **cardiac phase** between consecutive R peaks.

The phase is defined in the interval:

```
[-π , π]
```

with

```
θ = 0  at R peak
```

The phase is computed using:

```
calculate_linear_phase_ver2()
```

Optionally, **nonlinear phase alignment** using **Dynamic Time Warping (DTW)** can be applied to improve alignment between beats.

---

# 4. Mean ECG Morphology Extraction

The ECG signal is averaged in the phase domain using **phase bins**.

Example configuration:

```
Number of phase bins = 200
```

For each phase bin the algorithm computes:

- mean ECG amplitude
- standard deviation

Function used:

```
ecgsd_extractor_ver1()
```

Output variables:

```
ECGmean
ECGmean_nonlinear_phase
ECGsd
```

---

# 5. ECG Component Separation

The mean ECG waveform is divided into physiological regions.

### P wave region

```
-π/6 - π/2 < θ < -π/6
```

### QRS region

```
-π/6 < θ < π/6
```

### T wave region

```
π/6 < θ < π/6 + 3π/4
```

These intervals isolate the corresponding ECG components.

---

# 6. Gaussian Mixture Modeling of ECG Waves

Each ECG component is modeled using a **Gaussian mixture representation**.

A Gaussian kernel is defined as:

```
G(θ) = a_i * exp( - (θ − θ_i)^2 / (2 b_i^2) )
```

where

- `a_i` → amplitude  
- `b_i` → width  
- `θ_i` → phase location

Number of Gaussian kernels used:

```
P wave  → 3 Gaussians
QRS     → 5 Gaussians
T wave  → 3 Gaussians
```

---

# 7. Parameter Optimization

Gaussian parameters are estimated using **Particle Swarm Optimization (PSO)**.

Example configuration:

```
SwarmSize = 500 – 1300
MaxIterations = 200
```

The optimization minimizes the reconstruction error between the **mean ECG component** and the **Gaussian model**.

Estimated parameters include:

```
ai_Pwave, bi_Pwave, tetai_Pwave
ai_QRSwave, bi_QRSwave, tetai_QRSwave
ai_Twave, bi_Twave, tetai_Twave
```

---

# 8. Synthetic ECG Generation

Using the estimated Gaussian parameters, the algorithm reconstructs:

```
Synthetic_P
Synthetic_QRS
Synthetic_T
```

These are combined to generate the **synthetic ECG waveform**:

```
Synthetic_ECG
```

This synthetic signal approximates the original ECG morphology.

---

# 9. Peak Detection for P and T Waves

The synthetic P and T components are used to detect peaks across the ECG.

Steps:

1. Apply phase masks
2. Correlate the signal with a triangular window
3. Remove baseline
4. Detect peaks using `findpeaks`

Outputs:

```
idx_P_peak_final
idx_T_peak_final
```

---

# 10. Estimation of Wave Boundaries

The algorithm determines the onset and offset of each ECG wave using **energy distribution** of the synthetic components.

### P wave boundaries

```
P_on
P_off
```

Computed using cumulative energy of `Synthetic_P`.

---

### T wave boundaries

```
T_on
T_off
```

Computed using cumulative energy of `Synthetic_T`.

---

### `Synthetic_QRS` within a search window around each R_off
```

Computed using energy of `Synthetic_QRS` within a search window around each R peak.

Typical search window:

```
±50 ms around R peak
```

Energy thresholds:

```
0.1%  → onset
99%   → offset
```

---

# Visualization Outputs

The code produces multiple figures showing intermediate and final results:

- ECG signal with **R peaks**
- Mean ECG morphology
- Gaussian reconstruction of **P wave**
- Gaussian reconstruction of **QRS complex**
- Gaussian reconstruction of **T wave**
- Synthetic ECG vs original ECG
- Detected **P peaks and T peaks**
- Detected **P onset / offset**
- Detected **T onset / offset**
- Detected **QRS onset / offset**

These plots allow visual verification of the segmentation.

---

# Required MATLAB Toolboxes

The code requires:

- Signal Processing Toolbox
- Global Optimization Toolbox
- Statistics Toolbox

---

External dependency:

```
pantompkins_qrs.m
```

---

# Applications

This framework can be used in:

- ECG waveform segmentation
- biomedical signal analysis
- cardiac morphology modeling
- ECG simulation
- ECG compression research
- cardiac feature extraction

---

# Summary

This project presents a **model‑based ECG segmentation framework** that:

- aligns ECG beats in the **phase domain**
- models ECG morphology using **Gaussian mixtures**
- reconstructs **synthetic ECG components**
- uses these components to perform **robust segmentation of P, QRS, and T waves**

The result is a **physiology‑aware segmentation algorithm** capable of accurately detecting ECG wave boundaries across multiple heartbeats.
:::

