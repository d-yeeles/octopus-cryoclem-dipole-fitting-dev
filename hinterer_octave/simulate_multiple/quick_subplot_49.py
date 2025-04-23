import matplotlib.pyplot as plt

# Example data
inclination = [0.00, 1.84, 3.67, 5.51, 7.35, 9.18, 11.02, 12.86, 14.69, 16.53, 18.37, 20.20, 22.04, 23.88, 25.71, 27.55, 29.39, 31.22, 33.06, 34.90, 36.73, 38.57, 40.41, 42.24, 44.08, 45.92, 47.76, 49.59, 51.43, 53.27, 55.10, 56.94, 58.78, 60.61, 62.45, 64.29, 66.12, 67.96, 69.80, 71.63, 73.47, 75.31, 77.14, 78.98, 80.82, 82.65, 84.49, 86.33, 88.16]
inc_err = [0.00, 88.16, 76.28, 65.31, 71.86, 66.45, 78.98, 12.86, 14.69, 16.53, 18.37, 57.57, 57.03, 56.44, 25.71, 48.59, 47.46, 47.80, 45.87, 41.64, 42.70, 44.03, 39.28, 40.72, 34.23, 35.97, 33.09, 29.95, 30.02, 29.70, 27.29, 22.93, 24.30, 19.66, 18.15, 18.64, 16.41, 14.52, 12.34, 12.17, 10.29, 9.82, 8.60, 5.65, 7.21, 5.50, 2.77, 1.93, 0.49]
loc_err = [1.93, 9.55, 18.24, 23.14, 19.94, 17.13, 17.88, 0.68, 0.09, 2.81, 3.98, 15.93, 21.19, 10.09, 5.71, 17.71, 18.50, 20.90, 17.92, 16.44, 10.69, 12.01, 15.35, 8.72, 12.32, 14.77, 11.63, 8.36, 13.14, 20.28, 11.77, 13.31, 9.28, 4.66, 4.68, 9.35, 8.89, 10.06, 8.03, 5.02, 7.84, 4.76, 0.70, 3.82, 3.34, 3.78, 1.29, 1.24, 1.99]

# Create a figure with 3 subplots
fig, axs = plt.subplots(1, 3, figsize=(15, 5))

# Plot x vs y
axs[0].scatter(inclination, inc_err)
axs[0].set_xlabel('Inclination, °')
axs[0].set_ylabel('Inclination error, °')
axs[0].set_title('')

# Plot y vs z
axs[1].scatter(inclination, loc_err)
axs[1].set_xlabel('Inclination, °')
axs[1].set_ylabel('Localisation Error = mean(|X_err|, |Y_err|), nm')
axs[1].set_title('')

# Plot x vs z
axs[2].scatter(inc_err, loc_err)
axs[2].set_xlabel('Inclination error, °')
axs[2].set_ylabel('Localisation Error = mean(|X_err|, |Y_err|), nm')
axs[2].set_title('')

# Adjust layout to prevent overlap
plt.tight_layout()

# Show the plot
plt.show()

