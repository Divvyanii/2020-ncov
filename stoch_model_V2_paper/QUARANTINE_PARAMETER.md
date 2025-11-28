# Quarantine Parameter Addition

## Overview

A new **quarantine** parameter has been added to the stochastic SEIR model to account for the effect of quarantine/isolation measures on disease transmission.

## Parameters Added

Two new parameters have been added to `inputs/theta_initial_conditions.csv`:

1. **`quarantine_effectiveness`** (default: 0.5)
   - Type: Numeric (0-1)
   - Description: Proportion of infectious individuals that are effectively quarantined/isolated
   - Effect: Reduces the effective infectious population contributing to transmission
   - Example: 0.5 means 50% of infectious individuals are effectively isolated and don't transmit

2. **`quarantine_start_date`** (default: "2020-01-23")
   - Type: Date string (YYYY-MM-DD format)
   - Description: Date when quarantine measures begin
   - Effect: Quarantine only applies from this date onwards
   - Default: Same as travel restrictions date

## How It Works

### Mathematical Implementation

The quarantine parameter reduces the effective infectious population in the transmission equation:

**Before quarantine:**
```
S_to_E1 = susceptible × (infectious_t1 + infectious_t2) × inf_rate
```

**After quarantine:**
```
effective_inf1 = infectious_t1 × (1 - quarantine_effectiveness)
effective_inf2 = infectious_t2 × (1 - quarantine_effectiveness)
effective_tr_exp2 = tr_exposed_t2 × (1 - quarantine_effectiveness)  # pre-symptomatic

S_to_E1 = susceptible × (effective_inf1 + effective_inf2 + pre_symp×effective_tr_exp2) × inf_rate
```

### Time-Varying Implementation

- Quarantine is **inactive** before `quarantine_start_date`
- Quarantine becomes **active** on and after `quarantine_start_date`
- The effectiveness is constant once active (can be modified to be time-varying if needed)

## Files Modified

1. **`inputs/theta_initial_conditions.csv`**
   - Added `quarantine_effectiveness` and `quarantine_start_date` parameters

2. **`scripts/main_model.R`**
   - Added parameter loading from CSV
   - Added calculation of `quarantine_start_time` after data loading

3. **`R/model_functions.R`**
   - Modified `process_model()` function signature to accept `quarantineF` parameter
   - Updated transmission equation to use effective infectious population
   - Modified `smc_model()` to check quarantine timing and pass parameter to process model

4. **`R/load_timeseries_data.R`**
   - Added calculation of `quarantine_start_time` (with fallback default)

## Usage Example

To modify quarantine parameters, edit `inputs/theta_initial_conditions.csv`:

```csv
quarantine_effectiveness,0.7
quarantine_start_date,2020-01-20
```

This would mean:
- 70% of infectious individuals are effectively quarantined
- Quarantine measures begin on January 20, 2020

## Interpretation

- **quarantine_effectiveness = 0**: No quarantine effect (baseline model)
- **quarantine_effectiveness = 0.5**: 50% reduction in effective transmission
- **quarantine_effectiveness = 1.0**: Complete isolation (100% reduction, though this is unrealistic)

## Relationship to Other Parameters

The quarantine parameter works alongside:
- **Travel restrictions**: Both reduce transmission but through different mechanisms
- **R₀ decline**: Quarantine contributes to the overall reduction in transmission
- **Time-varying transmission**: Quarantine adds an additional layer of intervention

## Notes

- Quarantine affects both symptomatic (I1, I2) and pre-symptomatic (tr_exp2) transmission
- The parameter can be estimated via maximum likelihood or set based on intervention data
- Future enhancements could include:
  - Time-varying quarantine effectiveness
  - Different effectiveness for different compartments
  - Quarantine compliance/adherence parameters

