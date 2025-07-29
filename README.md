# MATLAB Code for Maximum-Likelihood-Estimator-Based Differentiator (MLEdiff) Using Quadratic and Zero-Order Optimal Splines

This repository contains MATLAB code supporting the results in:

**"Derivative Estimation from Coarse, Irregular, Noisy Samples: An MLE–Spline Approach"**  
by *Konstantin E. Avrachenkov* and *Leonid B. Freidovich*, 2025.

📧 Contact: [leonid.freidovich@umu.se](mailto:leonid.freidovich@umu.se)

---

## 📄 Description

This code implements a **maximum-likelihood-based numerical differentiation method** using **optimal smoothing quadratic and zero-order splines**. It is designed to estimate the derivative of a signal measured under **coarse, non-uniform sampling** and corrupted by **Gaussian noise**.

The implementation reproduces the key results and figures described in the above article, including comparisons with classical differentiators and performance under various noise and sampling conditions.

_Last updated: **2025-07-29**_

---

## 📦 Files

### 🔹 Primary Scripts
- ml_quadratic_spline_v16.m — Full implementation of the quadratic-spline-based nd zero-order-spline-based differentiators, including data generation, spline fitting, and comparison with two other methods (generates figures for the paper referenced above).
- ml_zero_spline_v1.m — Implementation of the zero-order spline-based differentiator.
- rec_mle_splines_v2.m — Full implementation of the recursive quadratic-spline-based and zero-order-spline-based differentiators, including data generation, spline fitting, and comparison with ttwo other methods.

### 🔹 Signal Generation
- make_data.m — Generates synthetic signals with non-uniform sampling and Gaussian noise.
- x_fun.m — Defines the test signal and its analytical derivative.

### 🔹 Spline-Based Differentiation
- quadratic_spline_step_QC.m — Constructs matrices for quadratic spline smoothing.
- quadratic_spline_step_new.m — Alternative formulation of the quadratic spline method.
- zero_order_spline_step_QC.m — Constructs matrices for zero-order spline smoothing.
- z_from_p.m — Computes derivative estimates from spline parameters.
- interp_z_from_z_K.m — Interpolates spline-based derivative estimates over time.

### 🔹 Other Differentiation Methods for Compsrision
- levant_step.m — Implements Levant’s first-order super-twisting differentiator (requires upper bound on second derivative), discretized using Euler's method.
- hgo_step.m — Implements a high-gain observer (HGO) for derivative estimation (requires tuning of a single small parameter, uses upper bound on the estimated first derivative for a protactive saturation), discretized using Euler's method.

### 🔹 Recursive and Update Utilities
- update_quadratic_spline.m — Efficient (optimized for speed in Matlab) recursive update of quadratic spline estimates using matrix identities.
- updateQC.m — Computes updates to second-order spline matrices for recursive estimation (for update of zero-order spline matrices see rec_mle_splines_v2.m).

---

## 🚀 How to Run

1. Clone or download this repository.
2. Open MATLAB and set the current folder to the repository root.
3. Run the main script:   ml_quadratic_spline_v16.m
  

---

## 🖋 Citation

If you use this code in your research or publications, please cite either the article, once published, or this depository: https://github.com/fleonid/MLEdiff/

---
