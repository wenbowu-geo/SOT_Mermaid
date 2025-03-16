import numpy as np
import matplotlib.pyplot as plt

# Define file paths
model_files = {
    "PREM_noCutoff": "results/H03N1/mode_1_2.50Hz.txt",
    "PREM_Cutoff": "../Modes_MultipleStations_GivenSeafloorCutoff/results_PREM/H03N1/mode_1_2.50Hz.txt",
    "PREM_Cutoff_Vs1000": "../Modes_MultipleStations_GivenSeafloorCutoff/results_PREMVpRho_Vs1000/H03N1/mode_1_2.50Hz.txt",
    "Vp4000Vs2000Rho2": "../Modes_MultipleStations_GivenSeafloorCutoff/results_Vp4000Vs2000Rho2/H03N1/mode_1_2.50Hz.txt"
}

def load_mode_data(file_path, complex_first_col=False):
    """
    Load mode data from a given file. If complex_first_col is True,
    the first column is treated as complex numbers (with 'i' replaced by 'j').
    """
    if complex_first_col:
        raw_data = np.genfromtxt(file_path, dtype=str)
        complex_col = np.array([complex(s.replace('i', 'j')) for s in raw_data[:, 0]])
        real_col = raw_data[:, 1].astype(float)
        return np.column_stack((np.abs(complex_col), real_col))  # Use absolute values of complex numbers
    else:
        return np.loadtxt(file_path)

# Load datasets
data_model_PREM_noCutoff = load_mode_data(model_files["PREM_noCutoff"])
data_model_PREM_Cutoff = load_mode_data(model_files["PREM_Cutoff"])
data_model_PREM_Cutoff_Vs1000 = load_mode_data(model_files["PREM_Cutoff_Vs1000"], complex_first_col=True)
data_model_Cutoff_Vp4000Vs1000Rho2 = load_mode_data(model_files["Vp4000Vs2000Rho2"], complex_first_col=True)

# Plot the data
plt.figure(figsize=(8, 5))
plt.plot(data_model_PREM_noCutoff[:, 0], data_model_PREM_noCutoff[:, 1], color='tab:blue', label="PREM")
plt.plot(data_model_PREM_Cutoff[:, 0], data_model_PREM_Cutoff[:, 1], color='tab:green', label="PREM_Cutoff")
plt.plot(data_model_PREM_Cutoff_Vs1000[:, 0], data_model_PREM_Cutoff_Vs1000[:, 1], color='tab:orange', label="PREM_Cutoff_Vs1000")
plt.plot(data_model_Cutoff_Vp4000Vs1000Rho2[:, 0], data_model_Cutoff_Vp4000Vs1000Rho2[:, 1], color='gray', label="Vp4000Vs2000Rho2")


# Determine global min and max amplitude values for horizontal line placement
all_x_values = np.concatenate([
    data_model_PREM_noCutoff[:, 0],
    data_model_PREM_Cutoff[:, 0],
    data_model_PREM_Cutoff_Vs1000[:, 0],
    data_model_Cutoff_Vp4000Vs1000Rho2[:,0]
])
max_x, min_x = np.max(all_x_values), np.min(all_x_values)
print(max_x, min_x )

# Plot a horizontal line at the last depth value of the PREM_noCutoff dataset
plt.plot([min_x,max_x],[data_model_PREM_Cutoff[-1, 1],data_model_PREM_Cutoff[-1, 1]], color='black', linestyle='--')

# Customize the plot
plt.xlabel("Normalized Amplitude")
plt.ylabel("Depth (m)")
plt.title("Ocean Acoustic Modes")
plt.legend()
plt.gca().invert_yaxis()  # Flip y-axis so depth increases downward

# Show the plot
plt.show()
#plt.savefig("modes.png",dpi=400)

