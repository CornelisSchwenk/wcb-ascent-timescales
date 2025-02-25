# WCB Moisture Transport Analysis

This repository contains the code and instructions to reproduce the results from the **2024 paper** titled:  
**"The Role of Ascent Timescale for WCB Moisture Transport into the UTLS"**  
by **Cornelis Schwenk**.

## Overview
This project involves **data analysis and visualization** using **Julia** and **Python** to study moisture transport by **Warm Conveyor Belts (WCBs)** into the **Upper Troposphere and Lower Stratosphere (UTLS)**.


![Visualization](assets/visualization.gif)


---

## Prerequisites
- **Julia Version:** `1.10.0`
- **Python Version:** `3.x`
- **Minimum System Requirements:**  
  - **64 GB RAM**  
  - **Estimated Runtime:** ~2 hours for plotting  

---

## Directory Structure
📁 **`EXECUTE_WCB_SELECTION/`** → Scripts for WCB selection  
📁 **`EXECUTE_WIND_SHEAR_CALCULATION/`** → Scripts for wind shear calculation  
📁 **`julia_code/`** → Julia scripts for data processing and plotting  
📁 **`python_code/`** → Python scripts for additional visualization  
📁 **`nest_data/`** → Folder for required NetCDF data files  
📁 **`plots/`** → Directory for storing generated plots  

---

## Setup & Execution Instructions

### **Step 0: Download Required Data**
Before running any scripts, you need to **download the required NetCDF files** from **Zenodo**:

📥 **Download Link:** [10.5281/zenodo.13913563](https://doi.org/10.5281/zenodo.13913563)  

Required files:
- `PECR.nc`
- `WCB_tau_vars_realtime.nc`
- `WCB_tau_vars_normed.nc`

After downloading, **place these files inside the `nest_data/` directory**.

💡 **Need Help?**  
For assistance with data access, contact:  
📧 **Cornelis Schwenk** - coschwen@uni-mainz.de  
📧 **Annette Miltenberger** - amiltenb@uni-mainz.de  

---

### **Step 1: Install Julia Dependencies**
Navigate to the `julia_code/` directory and install the necessary Julia packages:
```bash
julia --project=. install_packages.jl
```

---

### **Step 2: Generate Normalized WCB Variables**
Run the following command to create the `WCB_tau_vars_normed.nc` file:
```bash
julia --project=. make_normed_file.jl
```

---

### **Step 3: Generate Plots**
Execute the following scripts to generate the necessary plots:
```bash
julia --project=. make_plots.jl
julia --project=. pyplots.jl
```

---

### **Step 4: Generate Figure 3 (Python)**
Navigate to the `python_code/` directory, install dependencies, and generate **Figure 3**:
```bash
pip install -r requirements.txt
python plot_Fig3.py
```

---

## **(Optional) Full Workflow: Reproduce Entire Analysis**
⚠️ **The following steps require access to the full simulation data and should be run on the MOGON-II supercomputer.**

### **Step 1: WCB Selection**
Navigate to the `EXECUTE_WCB_SELECTION/` directory and submit the SLURM job:
```bash
sbatch mogon_merge_and_WCB_sel.slurm
```

### **Step 2: Wind Shear Calculation**
Navigate to the `EXECUTE_WIND_SHEAR_CALCULATION/` directory and execute:
```bash
sbatch do_calc_mean_wind_shear_mogon.slurm
```

### **Step 3: Write `PECR.nc` File**
Navigate to the `julia_code/` directory and run:
```bash
julia --project=. write_PECR.jl
```

---

## **Acknowledgments**
This research is conducted at **Johannes Gutenberg University Mainz** and is associated with the **MOGON-II Supercomputing Cluster**.

For questions, suggestions, or collaboration inquiries, please contact:  
📧 **Cornelis Schwenk** - coschwen@uni-mainz.de

