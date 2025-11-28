#!/usr/bin/env python3
"""
Generate quarantine effect visualization using Python
This creates a conceptual graph showing how quarantine affects transmission
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import os

# Set up the plot style
plt.style.use('seaborn-v0_8-darkgrid' if 'seaborn-v0_8-darkgrid' in plt.style.available else 'default')
fig, axes = plt.subplots(3, 1, figsize=(10, 12))
fig.suptitle('Effect of Quarantine on COVID-19 Transmission Dynamics', fontsize=14, fontweight='bold')

# Create date range
start_date = datetime(2019, 12, 15)
end_date = datetime(2020, 2, 5)
date_range = pd.date_range(start=start_date, end=end_date, freq='D')
travel_restriction_date = datetime(2020, 1, 23)
quarantine_start_date = datetime(2020, 1, 23)

# Read existing R0 data if available
try:
    r0_data = pd.read_csv('out_R0.csv', header=None)
    # Use median of first few columns as baseline
    baseline_r0 = r0_data.iloc[:, :10].median(axis=1).values
    # Interpolate to match date range length
    if len(baseline_r0) > len(date_range):
        baseline_r0 = baseline_r0[:len(date_range)]
    elif len(baseline_r0) < len(date_range):
        # Extend with last value
        baseline_r0 = np.append(baseline_r0, [baseline_r0[-1]] * (len(date_range) - len(baseline_r0)))
except:
    # Create synthetic baseline R0 if file doesn't exist
    baseline_r0 = np.ones(len(date_range)) * 2.5
    # Add some variation
    for i in range(len(baseline_r0)):
        if date_range[i] >= travel_restriction_date:
            baseline_r0[i] = 2.5 * np.exp(-0.02 * (i - np.where(date_range == travel_restriction_date)[0][0]))
        else:
            baseline_r0[i] = 2.5 + 0.3 * np.sin(i / 10)

# Quarantine levels to compare
quarantine_levels = [0, 0.3, 0.5, 0.7]
colors = ['#FF4444', '#FF8800', '#FFAA00', '#00AA00']  # Red, Orange, Yellow, Green
labels = ['0%', '30%', '50%', '70%']

# Calculate R0 for each quarantine scenario
r0_scenarios = {}
for q_level, color, label in zip(quarantine_levels, colors, labels):
    r0_values = baseline_r0.copy()
    # Apply quarantine effect after quarantine start date
    quarantine_start_idx = np.where(date_range >= quarantine_start_date)[0]
    if len(quarantine_start_idx) > 0:
        start_idx = quarantine_start_idx[0]
        # Reduce R0 by quarantine effectiveness after start date
        r0_values[start_idx:] = r0_values[start_idx:] * (1 - q_level)
    r0_scenarios[label] = r0_values

# Panel 1: Reproduction Number (R_t)
ax1 = axes[0]
for label, r0_vals in r0_scenarios.items():
    ax1.plot(date_range, r0_vals, label=f'Quarantine: {label}', linewidth=2.5, 
             color=colors[labels.index(label)])

ax1.axvline(travel_restriction_date, color='red', linestyle='--', linewidth=2, label='Travel restrictions')
ax1.axhline(1, color='gray', linestyle=':', linewidth=1.5, alpha=0.7)
ax1.set_ylabel('Reproduction Number (Rₜ)', fontsize=11, fontweight='bold')
ax1.set_title('Effect of Quarantine on Reproduction Number', fontsize=12, fontweight='bold')
ax1.set_ylim(0, 8)
ax1.legend(loc='upper right', fontsize=9)
ax1.grid(True, alpha=0.3)
ax1.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
ax1.xaxis.set_major_locator(mdates.WeekdayLocator(interval=2))

# Panel 2: Daily New Cases (simulated based on R0)
ax2 = axes[1]
case_scenarios = {}
initial_cases = 1
for label, r0_vals in r0_scenarios.items():
    cases = [initial_cases]
    for i in range(1, len(r0_vals)):
        # Simple exponential growth model
        growth_rate = (r0_vals[i-1] - 1) / 5.2  # Assuming ~5 day generation time
        new_cases = cases[-1] * np.exp(growth_rate)
        cases.append(new_cases)
    # Convert to daily new cases
    daily_cases = np.diff([0] + cases)
    daily_cases = np.maximum(daily_cases, 0)  # Ensure non-negative
    case_scenarios[label] = daily_cases
    ax2.plot(date_range, daily_cases, label=f'Quarantine: {label}', linewidth=2.5,
             color=colors[labels.index(label)])

ax2.axvline(travel_restriction_date, color='red', linestyle='--', linewidth=2)
ax2.axvline(quarantine_start_date, color='blue', linestyle='--', linewidth=2, label='Quarantine starts')
ax2.set_ylabel('Daily New Cases in Wuhan', fontsize=11, fontweight='bold')
ax2.set_title('Effect of Quarantine on Daily Case Incidence', fontsize=12, fontweight='bold')
ax2.legend(loc='upper left', fontsize=9)
ax2.grid(True, alpha=0.3)
ax2.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
ax2.xaxis.set_major_locator(mdates.WeekdayLocator(interval=2))

# Panel 3: Infection Prevalence (cumulative)
ax3 = axes[2]
prevalence_scenarios = {}
for label, daily_cases in case_scenarios.items():
    # Cumulative prevalence (simplified)
    prevalence = np.cumsum(daily_cases) * 0.1  # Scale factor for prevalence
    prevalence_scenarios[label] = prevalence
    ax3.plot(date_range, prevalence, label=f'Quarantine: {label}', linewidth=2.5,
             color=colors[labels.index(label)])

ax3.axvline(travel_restriction_date, color='red', linestyle='--', linewidth=2)
ax3.axvline(quarantine_start_date, color='blue', linestyle='--', linewidth=2)
ax3.set_ylabel('Infection Prevalence (E+I)', fontsize=11, fontweight='bold')
ax3.set_xlabel('Date', fontsize=11, fontweight='bold')
ax3.set_title('Effect of Quarantine on Infection Prevalence', fontsize=12, fontweight='bold')
ax3.legend(loc='upper left', fontsize=9)
ax3.grid(True, alpha=0.3)
ax3.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
ax3.xaxis.set_major_locator(mdates.WeekdayLocator(interval=2))

# Rotate x-axis labels
for ax in axes:
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')

plt.tight_layout()

# Ensure plots directory exists
os.makedirs('plots', exist_ok=True)

# Save as PDF
output_file = 'plots/quarantine_effect_quick.pdf'
plt.savefig(output_file, format='pdf', dpi=300, bbox_inches='tight')
print(f"✓ Graph saved to: {output_file}")

# Also save as PNG for quick viewing
output_png = 'plots/quarantine_effect_quick.png'
plt.savefig(output_png, format='png', dpi=300, bbox_inches='tight')
print(f"✓ Also saved as PNG: {output_png}")

plt.close()

print("\nGraph generation complete!")
print("The PDF shows how different quarantine effectiveness levels (0%, 30%, 50%, 70%)")
print("affect transmission dynamics, case incidence, and infection prevalence.")

