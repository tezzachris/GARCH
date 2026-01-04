# GARCH2F Model  
*A GARCH Model with Two Volatility Components and Two Driving Factors*

---

## Overview

This repository contains MATLAB code for the estimation, simulation, and option pricing of a **GARCH model with two innovation factors and two volatility components** (GARCH2F).

The model extends standard GARCH-type specifications by allowing for multiple volatility dynamics and driving shocks, providing increased flexibility in capturing empirical features of financial return series.

---

## Model Description

The GARCH2F model includes:

- Two innovation (shock) factors  
- Two conditional volatility components  

This structure allows the model to better capture volatility persistence, clustering, and heterogeneous shock transmission mechanisms observed in financial markets.

---

## Model Estimation

The estimation of model parameters is non-trivial due to the high dimensionality and nonlinearity of the likelihood function.

The current implementation is written in **MATLAB** and relies on built-in nonlinear optimization routines:

- `fmincon`
- `fminunc`

To improve convergence and mitigate sensitivity to initial conditions, a **starting parameter routine** is provided. This routine performs a standard **grid search** to identify suitable initial parameter values prior to likelihood maximization.

**Notes:**

- Ideally, the estimation procedure should be repeated multiple times using different starting points to gain confidence in the results.
- In our experience:
  - `fminunc` often delivers more accurate estimates
  - However, it requires additional care to enforce parameter constraints explicitly

---

## Model Simulation

The simulation module allows for generating synthetic return paths from the GARCH2F model.

Simulation requires:

- An initial value for the conditional volatility components  
- An initial value for the return process (optional)  

These simulations can be used for validation, stress testing, and option pricing applications.

---

## Option Pricing via Characteristic Function

European call and put option prices are computed using the **characteristic function approach**.

Specifically, option prices are obtained via the inversion formula of **Gil-Pelaez**, following Equation (11) in Reference [1]. This approach allows for semi-closed-form option valuation under the GARCH2F dynamics.

---

## Data

- **Risk-free rate**:  
  Freely available from the Federal Reserve Economic Data (FRED):  
  https://fred.stlouisfed.org/series/TB3MS

- **Equity and option data**:  
  S&P 500 returns and option data were obtained via **Refinitiv Eikon Datastream** under a University license.

---

## File Structure

- `start_param.m % Randomized starting values`
- `loglik.m % Kalman filter log-likelihood`
- `RungeKuttaFuture.m % Runge–Kutta solver for affine coefficients used in log-future price formula`
- `main.m % Main estimation routine`

---

## Reference

If you use this code, please cite:
   ```matlab
      Ballestra, L. V., D’Innocenzo, E., & Tezza, C.  
      "A GARCH Model with Two Volatility Components and Two Driving Factors"
      Journal of Empirical Finance**,  
      Volume 85, February 2026, Article 101671
