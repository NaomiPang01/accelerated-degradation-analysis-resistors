# accelerated-degradation-analysis-resistors
Reliability Engineering_Final Project

# 🔬 Reliability Prediction of Carbon-Film Resistors using Accelerated Degradation Models

[![R](https://img.shields.io/badge/Language-R-blue.svg)](https://www.r-project.org/)
[![Status](https://img.shields.io/badge/Status-Completed-success.svg)]()

This repository contains the final project for the **Reliability Analysis** course, focusing on the two-stage accelerated degradation testing (ADT) and lifetime prediction of carbon-film resistors.

---

## 📂 Project Overview

* **Data Source**: Table C.3 (Percent Increase in Resistance Over Time of Carbon-Film Resistors) from Meeker & Escobar (1998) via `SMRD.data`.
* **Accelerating Factor**: Temperature ($83^\circ\text{C}$, $133^\circ\text{C}$, $173^\circ\text{C}$).
* **Soft Failure Threshold**: $5\%$ resistance increase.
* **Authors**: 潘皓心、黃芯鈺 — National Cheng Kung University (NCKU).

---

## 🛠️ Methodology & Workflow

The analysis follows a robust **Two-Stage Degradation Modeling** approach:

1. **Exploratory Non-Parametric Analysis**:
   - Utilized **Kaplan-Meier (KM)** estimator to visualize degradation trends and verify the significant accelerating effect of temperature.
2. **Stage 1: Degradation Path Modeling**:
   - Fitted a **Power Law degradation model** ($\log D(t) = \beta_0 + \beta_1 \log t + \varepsilon$) for individual units.
   - Extrapolated to find the **Pseudo-Failure Times ($T_f$)** at the $5\%$ failure threshold.
3. **Stage 2: Accelerated Life Testing (ALT) Modeling**:
   - Fitted a **Weibull–Arrhenius model** using Maximum Likelihood Estimation (MLE) via `survreg`.
   - Estimated activation energy ($E_a$) and Weibull shape/scale parameters.
4. **Uncertainty Analysis & Predictions**:
   - Evaluated life metrics at normal use condition ($50^\circ\text{C}$), including **MTTF**, **B10 Life**, and failure probabilities.
   - Compared **Normal Approximation (Delta Method)** and **Bootstrap (BCa/Percentile)** confidence intervals.
   - Performed **Reverse Prediction** to calculate the maximum allowable operating temperature for a 5-year B10 design requirement.

---

## 🚀 Code Structure (`scripts/TableC3.R`)

* **Step 1**: Data loading, cleaning, and structure preparation.
* **Step 2**: Spaghetti plots and pseudo-failure time extraction (Power Law).
* **Step 3**: Weibull-Arrhenius MLE modeling (`survreg`).
* **Step 4**: Use-condition reliability metrics prediction ($50^\circ\text{C}$).
* **Step 5**: Asymptotic Normal Approximation CI (Delta Method).
* **Step 6**: Two-stage Bootstrap uncertainty and bias correction.
* **Step 7**: Reverse temperature optimization ($B_{10} \ge 5\text{ years}$).

---

## 📦 Required R Packages

Ensure the following packages are installed before running the script:

```R
install.packages(c("tidyverse", "survival", "SMRD.data", "ggfortify", "boot"))
