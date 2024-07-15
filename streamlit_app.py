import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from scipy.signal import savgol_filter, windows
from scipy.fft import fft, fftfreq
import os

# Define checkbox keys
checkbox_keys = [
    "Fluorescence (Normalized)",
    "Derivative of Normalized Fluorescence",
    "Fluorescence",
    "Transmission Sample",
    "Derivative of Transmission Sample",
    "Transmission Foil",
    "Derivative of Transmission Foil",
    "k-space",
    "R-space",
]

# Function to remember checkbox state using session state
def remember_checkbox_state():
    for key in checkbox_keys:
        if key not in st.session_state:
            st.session_state[key] = False

# Define plots for all graphs
def plot_graph(x, y, x_label, y_label, title, legend_labels):
    fig = go.Figure()
    for i, y_data in enumerate(y):
        fig.add_trace(go.Scatter(x=x[i], y=y_data, mode='lines', name=legend_labels[i]))
    fig.update_layout(title=title, xaxis_title=x_label, yaxis_title=y_label)
    return fig

# Function to open, read, and modify files
def read_modify_files(file_paths, choices):
    all_dfs = []
    legend_labels = []
    combined_plots = {}
    
    for path in file_paths:
        with open(path, "r") as f:
            lines = f.readlines()
        lines[2]

    for path in file_paths:
        try:
            # Read and process file
            with open(path, "r") as f:
                lines = f.readlines()
            
            lines[14] = lines[14].replace("#C", "")
            
            temp_path = path + "_temp"
            with open(temp_path, "w") as temp_file:
                temp_file.writelines(lines)
            
            data_start = 0
            for i in range(len(lines)):
                if not lines[i].startswith("#"):
                    data_start = i
                    break
            
            df = pd.read_csv(temp_path, skiprows=data_start, delimiter="\s+")
            df["Normalization"] = df["10_flatot"] / df["3_i0"]

            energy_diff = np.diff(df['1_Energy'])
            normalization_diff = np.diff(df['Normalization'])
            slopes = normalization_diff / energy_diff
            mean_slope = np.mean(slopes)
            std_slope = np.std(slopes)

            significant_change_indices = np.where(slopes > mean_slope + 3 * std_slope)[0]

            if significant_change_indices.size > 0:
                shoot_up_index = significant_change_indices[0]
                shoot_up_energy = df['1_Energy'].iloc[shoot_up_index]
                ##st.write(f"Edge begins at energy level: {shoot_up_energy:.2f} eV")
            else:
                shoot_up_energy = None

            df['Smoothed_Normalization'] = savgol_filter(df['Normalization'], window_length=60, polyorder=3)
            
            if 1 in choices:
                if "Fluorescence (Normalized)" not in combined_plots:
                    combined_plots["Fluorescence (Normalized)"] = {"x": [], "y": [], "labels": []}
                combined_plots["Fluorescence (Normalized)"]["x"].append(df["1_Energy"])
                combined_plots["Fluorescence (Normalized)"]["y"].append(df["Normalization"])
                combined_plots["Fluorescence (Normalized)"]["labels"].append(path)

            if 2 in choices:
                derivative = normalization_diff / energy_diff
                df['derivative'] = np.append(derivative, np.nan)
                if "Derivative of Normalized Fluorescence" not in combined_plots:
                    combined_plots["Derivative of Normalized Fluorescence"] = {"x": [], "y": [], "labels": []}
                combined_plots["Derivative of Normalized Fluorescence"]["x"].append(df["1_Energy"])
                combined_plots["Derivative of Normalized Fluorescence"]["y"].append(df["derivative"])
                combined_plots["Derivative of Normalized Fluorescence"]["labels"].append(path)
                
            if 3 in choices:
                if "Fluorescence" not in combined_plots:
                    combined_plots["Fluorescence"] = {"x": [], "y": [], "labels": []}
                combined_plots["Fluorescence"]["x"].append(df["1_Energy"])
                combined_plots["Fluorescence"]["y"].append(df["10_flatot"])
                combined_plots["Fluorescence"]["labels"].append(path)
                
            if 4 in choices:
                df['data'] = df["7_mu01"]
                if "Transmission Sample" not in combined_plots:
                    combined_plots["Transmission Sample"] = {"x": [], "y": [], "labels": []}
                combined_plots["Transmission Sample"]["x"].append(df["1_Energy"])
                combined_plots["Transmission Sample"]["y"].append(df["data"])
                combined_plots["Transmission Sample"]["labels"].append(path)

            if 5 in choices:
                data_diff = np.diff(df['data'])
                derivative2 = data_diff / energy_diff
                df['derivative2'] = np.append(derivative2, np.nan)
                if "Derivative of Transmission Sample" not in combined_plots:
                    combined_plots["Derivative of Transmission Sample"] = {"x": [], "y": [], "labels": []}
                combined_plots["Derivative of Transmission Sample"]["x"].append(df["1_Energy"])
                combined_plots["Derivative of Transmission Sample"]["y"].append(df["derivative2"])
                combined_plots["Derivative of Transmission Sample"]["labels"].append(path)

            if 6 in choices:
                df['reference_foil'] = df["8_mu12"]
                if "Transmission Foil" not in combined_plots:
                    combined_plots["Transmission Foil"] = {"x": [], "y": [], "labels": []}
                combined_plots["Transmission Foil"]["x"].append(df["1_Energy"])
                combined_plots["Transmission Foil"]["y"].append(df["reference_foil"])
                combined_plots["Transmission Foil"]["labels"].append(path)

            if 7 in choices:
                reference_foil_diff = np.diff(df['reference_foil'])
                derivative3 = reference_foil_diff / energy_diff
                df['derivative3'] = np.append(derivative3, np.nan)
                if "Derivative of Transmission Foil" not in combined_plots:
                    combined_plots["Derivative of Transmission Foil"] = {"x": [], "y": [], "labels": []}
                combined_plots["Derivative of Transmission Foil"]["x"].append(df["1_Energy"])
                combined_plots["Derivative of Transmission Foil"]["y"].append(df["derivative3"])
                combined_plots["Derivative of Transmission Foil"]["labels"].append(path)

            all_dfs.append(df)
            legend_labels.append(path)

        except Exception as e:
            st.write(f"An error occurred while processing {path}: {e}")
        finally:
            if os.path.exists(temp_path):
                os.remove(temp_path)

    for key, data in combined_plots.items():
        fig = plot_graph(data["x"], data["y"], "Energy Level (eV)", key, key, data["labels"])
        st.plotly_chart(fig)
        
    for path, df in zip(file_paths, all_dfs):
        try:
            if 8 in choices:
                if shoot_up_energy is not None:
                    df['Oscillatory_Part'] = df['Normalization'] - df['Smoothed_Normalization']
                    df['Normalized_Oscillatory_Part'] = df['Oscillatory_Part'] / df['Smoothed_Normalization']

                    E0 = shoot_up_energy
                    m_e = 9.10938356e-31  # mass of an electron in kg
                    hbar = 1.0545718e-34  # planks constant in Js
                    eV_to_J = 1.60218e-19  # eV to Joules
                    
                    df['k_space'] = np.sqrt((2 * m_e * (df['1_Energy'] - E0) * eV_to_J) / hbar)
                    df = df[df['k_space'].notna() & (df['k_space'] >= 0)]  # Filter invalid values
                    df['k2_Xk'] = df['k_space']**2 * df['Normalized_Oscillatory_Part']
                    
                    fig = plot_graph([df['k_space']], [df['k2_Xk']], "wavenumber (Å^-1)", "k^2 * X(k) (Å^-2)", f"k-space - {path}", [path])
                    st.plotly_chart(fig)
                else:
                    st.write(f"Cannot compute k-space for {path} because shoot_up_energy is not defined.")
            
            if 9 in choices and 8 in choices:
                if 'k_space' in df.columns and 'k2_Xk' in df.columns:
                    k_space = df['k_space'].values
                    k2_Xk = df['k2_Xk'].values
                    N = len(k_space)
                    T = k_space[1] - k_space[0]

                    # Windowed FFT
                    w = windows.hann(N)
                    ywf = fft(k2_Xk * w)
                    xf = fftfreq(N, T)[:N//2]
                    
                    fig = go.Figure()
                    fig.add_trace(go.Scatter(x=xf[1:N//2], y=2.0/N * np.abs(ywf[1:N//2]), mode='lines', name='|X(R)|'))
                    fig.add_trace(go.Scatter(x=xf[1:N//2], y=2.0/N * ywf.real[1:N//2], mode='lines', name='Real part of X(R)'))
                    fig.add_trace(go.Scatter(x=xf[1:N//2], y=2.0/N * ywf.imag[1:N//2], mode='lines', name='Imaginary part of X(R)'))
                    fig.update_layout(title=f'R-space - {path}', xaxis_title='Radial distance (Å)', yaxis_title='|X(R)|, Real part, and Imaginary part (Å^-3)')
                    st.plotly_chart(fig)
                else:
                    st.write(f"Cannot compute R-space for {path} because k_space and/or k2_Xk are not defined.")

        except Exception as e:
            st.write(f"An error occurred while processing {path} for k-space or r-space: {e}")

#streamlit app 
def main():
    st.title("XAS Data Visualization")
    try:
        uploaded_files = st.file_uploader("Upload your raw data file(s)", accept_multiple_files=True)
        
    except Exception as e: 
        st.write(f"Incorrect file type uploaded. Please try another file. {e}")

    if uploaded_files:
        choices = [i + 1 for i, key in enumerate(checkbox_keys) if st.session_state.get(key, False)]
        file_paths = []
        for uploaded_file in uploaded_files:
            file_path = f"temp_{uploaded_file.name}"
            with open(file_path, "wb") as f:
                f.write(uploaded_file.getbuffer())
            file_paths.append(file_path)
            
        remember_checkbox_state()
        
        # Display checkboxes and add "All of the above" option
        st.write("Select which graphs to generate:")
        st.caption("When choosing a derivaitive graph, please also choose its non-derivative counterpart.")
        st.caption("When choosing R-space, please also choose k-space.")
        if st.checkbox("All graphs"):
            for key in checkbox_keys:
                st.session_state[key] = True
        for key in checkbox_keys:
            st.session_state[key] = st.checkbox(key, st.session_state.get(key, False))

        # Generate graphs button
        if st.button("Generate Graphs"):
            read_modify_files(file_paths, choices)

            for file_path in file_paths:
                if os.path.exists(file_path):
                    os.remove(file_path)

if __name__ == "__main__":
    main()
