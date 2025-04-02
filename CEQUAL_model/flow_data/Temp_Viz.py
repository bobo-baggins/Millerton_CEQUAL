import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# Read the CSV files
k2p_temp = pd.read_csv('K2P_Temp.csv')
k2p_weighted = pd.read_csv('K2P_Weighted_Temp.csv')
sja_temp = pd.read_csv('SJA_Temp.csv')

# Convert first column to datetime and handle NaT values
for df in [k2p_temp, k2p_weighted, sja_temp]:
    df.iloc[:, 0] = pd.to_datetime(df.iloc[:, 0])
    # Remove any rows with NaT values
    df.dropna(subset=[df.columns[0]], inplace=True)

# Create the plot
plt.figure(figsize=(12, 6))

# Plot each temperature series
plt.plot(k2p_temp.iloc[:, 0].values, k2p_temp.iloc[:, 1].values, label='K2P Temperature')
plt.plot(k2p_weighted.iloc[:, 0].values, k2p_weighted.iloc[:, 1].values, label='Weighted Temperature')
plt.plot(sja_temp.iloc[:, 0].values, sja_temp.iloc[:, 1].values, label='SJA Temperature')

# Format x-axis to show dates nicely
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
plt.gca().xaxis.set_major_locator(mdates.MonthLocator())

# Add vertical dashed line at July 30th with label
july_30 = pd.to_datetime('2024-07-30').date()
plt.axvline(x=july_30, color='black', linestyle='--', alpha=0.7)

# Add boxed label with arrow
plt.annotate('K2P Data Available\nbeginning July 30th', 
            xy=(july_30, 19), # Point on the line
            xytext=(july_30 - pd.Timedelta(days=60), 20), # Text position moved further right
            bbox=dict(facecolor='white', edgecolor='black', alpha=0.8),
            horizontalalignment='left',
            verticalalignment='top',
            arrowprops=dict(arrowstyle='-', color='black', linewidth=1))

# Customize the plot
plt.title('Inflow Temperature Comparison')
plt.xlabel('Date')
plt.ylabel('Temperature (*C)')
plt.legend()
plt.grid(True)

# Rotate x-axis labels for better readability
plt.xticks(rotation=45)

# Adjust layout to prevent label cutoff
plt.tight_layout()

# Show the plot
plt.show()
