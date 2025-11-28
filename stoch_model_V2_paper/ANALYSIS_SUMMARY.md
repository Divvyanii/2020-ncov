# Analysis Summary: stoch_model_V2_paper

## Figure 2 Description

**Figure 2** (saved as `plots/Figure_2.pdf`) is a multi-panel figure showing the results of the stochastic SEIR model fitted to early COVID-19 outbreak data. The figure contains **7 panels** arranged in a 3×3 grid layout:

### Panel A: Time-Varying Reproduction Number (Rₜ)
- **Y-axis**: Rₜ (reproduction number) ranging from 0 to 8
- **X-axis**: Date from December 15, 2019 to February 5, 2020
- Shows:
  - Median Rₜ (blue line)
  - 50% credible interval (IQR) - darker blue shaded region
  - 95% credible interval - lighter blue shaded region
  - Red vertical line at January 23, 2020 (Wuhan travel restrictions)
  - Horizontal dashed line at Rₜ = 1 (epidemic threshold)
- **Key finding**: Rₜ declines from ~2.5 before travel restrictions to ~1.05 after restrictions

### Panel B: New Onsets in Wuhan
- **Y-axis**: Number of new onsets (0-25)
- **X-axis**: Date from December 15, 2019 onwards
- Shows:
  - Model estimates (blue lines and shaded regions)
  - Observed data points:
    - Triangles (Wuhan data)
    - Diamonds (China data)
  - Red vertical line marking travel restrictions
- Model fits to local case onset data

### Panel C: New International Onsets
- **Y-axis**: Number of new international onsets (0-10)
- **X-axis**: Date from December 15, 2019 onwards
- Shows:
  - Model estimates (blue)
  - Observed international case data (black points)
  - Red vertical line for travel restrictions
- Validates model predictions against exported cases

### Panel D: Prevalence of Pre-symptomatic Infections in Wuhan
- **Y-axis**: Prevalence proportion (0-15%)
- **X-axis**: Date from December 15, 2020 onwards
- Shows:
  - Model estimates of pre-symptomatic prevalence
  - Validation data from evacuation flights:
    - Japan flights (3 data points with confidence intervals)
    - Germany/Korea flight
    - Singapore flight
    - Italy flight
    - Malaysia flight
- Each flight data point shows proportion infected with binomial confidence intervals

### Panel E: New Cases in Wuhan
- **Y-axis (left)**: Model estimates of new cases (0-50,000)
- **Y-axis (right)**: Confirmed cases (0-5,000)
- **X-axis**: Date from December 15, 2019 onwards
- Shows:
  - Model estimates (blue)
  - Observed confirmed cases (black points, right axis)
- Dual y-axis allows comparison of model predictions with confirmed case data

### Panel F: New International Exports Confirmed
- **Y-axis**: Number of new international exports confirmed (0-10)
- **X-axis**: Date from December 15, 2019 onwards
- Shows:
  - Model estimates (blue)
  - Observed data (black points, labeled as "non-fitted" for validation)
  - Red vertical line for travel restrictions
- Used for model validation (not fitted)

### Panel G: Expected vs Observed International Exports
- **X-axis**: Expected international exports from Wuhan
- **Y-axis**: Confirmed international exports
- Shows:
  - Scatter plot of expected vs observed
  - 1:1 line (dashed) for reference
  - Points represent different countries
- Validates model's ability to predict country-specific export patterns

---

## Monte Carlo Simulation Summary

### Method: Sequential Monte Carlo (SMC) Particle Filtering

The analysis uses **Sequential Monte Carlo (SMC)**, also known as **particle filtering**, to estimate the transmission dynamics and latent states of the COVID-19 outbreak.

### Key Components

#### 1. **Particle System**
- **Number of particles**: 1,000-2,000 particles per simulation
- **Bootstrap replicates**: 100 independent runs for uncertainty quantification
- Each particle represents a complete trajectory of the epidemic

#### 2. **Process Model (SEIR Structure)**
The stochastic SEIR model includes:
- **Compartments**:
  - Susceptible (S)
  - Exposed (E) - split into E1 and E2 for two-stage incubation
  - Infectious (I) - split into I1 and I2 for two-stage recovery
  - Recovered (R)
  - Travel compartments (tr_exp1, tr_exp2, tr_waiting) for exported cases
  - Local reporting compartments (waiting_local, cases_local, reports_local)

- **Time step**: dt = 0.25 days (4 steps per day)
- **Stochastic transmission**: Uses Euler-Maruyama numerical integration
- **Time-varying transmission**: Random walk on transmission rate (β) with volatility parameter

#### 3. **SMC Algorithm Steps**

For each time step t = 2 to T:

**Step 1: Prediction**
- Evolve each particle forward using the process model
- Add random walk noise to transmission rate: `β[t] = β[t-1] × exp(ε[t])` where ε ~ N(0, σ_vol)
- Account for travel restrictions (set travel fraction to 0 after Jan 23, 2020)

**Step 2: Weight Assignment**
Calculate likelihood weights for each particle based on:
- **Local case onsets** (Wuhan): Poisson likelihood
- **International case onsets**: Poisson likelihood (summed across countries)
- **Local confirmed cases**: Negative binomial likelihood (accounts for overdispersion)
- **Evacuation flight data**: Binomial likelihood (proportion infected on flights)
- **Country-specific exports**: Poisson likelihood for each destination country

Weight formula: `w[i,t] = exp(log-likelihood[i,t])`

**Step 3: Normalization**
- Normalize weights: `W[i,t] = w[i,t] / Σw[j,t]`

**Step 4: Resampling**
- Sample parent particles according to weights: `A[i,t] ~ Multinomial(W[1:nn,t])`
- Resample state variables: `X[i,t] = X[A[i,t],t]`
- This prevents particle degeneracy (few particles with high weight)

**Step 5: Backward Sampling**
- After forward pass, sample a single trajectory backward through the particle ancestry
- This provides a single representative trajectory for each bootstrap run

#### 4. **Likelihood Calculation**

The likelihood combines multiple data sources:

1. **Local onset data**:
   ```
   λ_local = cases_local × confirmed_prop × onset_prop × local_rep_prop
   log-lik = dpois(observed, λ_local)
   ```

2. **International onset data**:
   ```
   λ_int = cases_exported × confirmed_prop × onset_prop_int × country_risk
   log-lik = dpois(observed_total, sum(λ_int))
   ```

3. **Local confirmed cases** (negative binomial for overdispersion):
   ```
   μ = reports_local × confirmed_prop × local_rep_prop
   log-lik = dnbinom(observed, mu=μ, size=1/rep_local_var)
   ```

4. **Evacuation flights**:
   ```
   p_inf = prevalence / population_travel
   log-lik = dbinom(infected, total_passengers, p_inf)
   ```

5. **Country-specific exports**:
   ```
   λ_country = reports_exported × confirmed_prop × country_risk
   log-lik = Σ dpois(observed_country, λ_country)
   ```

#### 5. **Bootstrap Procedure**

To quantify uncertainty:
1. Run 100 independent SMC simulations
2. Each simulation uses different random seeds
3. Extract trajectories for:
   - Susceptible population (S)
   - Exposed population (E)
   - Infectious population (I)
   - Local cases (C_local)
   - International cases (C)
   - Reproduction number (R₀ = β / recovery_rate)
4. Calculate quantiles (2.5%, 25%, 50%, 75%, 97.5%) across bootstrap runs
5. Plot credible intervals and median trajectories

#### 6. **Key Parameters Estimated**

- **R₀ (time-varying)**: Basic reproduction number, estimated daily
- **Transmission rate (β)**: Time-varying, follows random walk
- **Latent states**: S, E, I, R populations over time
- **Case counts**: Local and international cases (with reporting delays)

#### 7. **Model Assumptions**

- **Initial conditions**: 1 infectious case on Nov 22, 2019
- **Population at risk**: 11 million (Wuhan travel population)
- **Incubation period**: Mean 5.2 days (gamma distribution, shape=2)
- **Infectious period**: Mean 2.9 days (gamma distribution, shape=2)
- **Reporting delay**: Mean 6.1 days (exponential)
- **Reporting proportions**:
  - Local reporting: 0.66% of cases
  - Onset known: 16% locally, 47% internationally
  - Confirmed proportion: 100% (all symptomatic cases)

#### 8. **Computational Details**

- **Parallel processing**: Uses `foreach` and `doMC` with 4 CPU cores
- **Time period**: Nov 22, 2019 to end of February 2020
- **Total time steps**: ~100 days × 4 steps/day = 400 steps
- **Total particles**: 1,000-2,000 per run
- **Total computations**: 100 bootstrap runs × 1,000 particles × 400 steps = 40 million particle-time steps

### Results Summary

From the bootstrap analysis:

**R₀ Estimates**:
- **Before travel restrictions** (Jan 1-23, 2020): Median = 2.35 (95% CI: 1.15-4.77)
- **After travel restrictions** (Jan 31, 2020): Median = 1.05 (95% CI: 0.41-2.39)
- **Overall period** (Jan 1-23): Median = 2.15 (95% CI: 1.56-2.57)

**Key Findings**:
1. R₀ declined significantly after travel restrictions
2. Model successfully fits both local and international case data
3. Pre-symptomatic prevalence estimates validated by evacuation flight data
4. Model captures country-specific export patterns

### Advantages of SMC Approach

1. **Handles non-linear dynamics**: Can model complex stochastic processes
2. **Time-varying parameters**: Estimates R₀ and transmission rate over time
3. **Multiple data sources**: Combines local cases, international exports, and flight data
4. **Uncertainty quantification**: Bootstrap provides credible intervals
5. **Latent state estimation**: Estimates unobserved populations (E, I)
6. **Real-time inference**: Can update estimates as new data arrives

### Technical Implementation

- **Language**: R
- **Key packages**: `foreach`, `doMC` (parallel processing)
- **Numerical method**: Euler-Maruyama for stochastic differential equations
- **Resampling**: Multinomial resampling
- **Likelihood**: Poisson, negative binomial, and binomial distributions

---

## Files Generated

- `plots/Figure_2.pdf`: Main results figure
- `outputs/bootstrap_fit_1.RData`: Bootstrap results (100 runs)
- `outputs/case_model.csv`: Model-estimated case counts
- `outputs/before_after_R.csv`: R₀ estimates before/after restrictions
- `out_R0.csv`: Time series of R₀ estimates
- `out_date.csv`: Corresponding dates

---

## Reference

Kucharski AJ, Russell TW, Diamond C et al. Early dynamics of transmission and control of 2019-nCoV: a mathematical modelling study. Lancet Infectious Diseases, 2020.



