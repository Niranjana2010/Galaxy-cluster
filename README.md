# Galaxy-cluster
# Estimating the Dynamical Mass of a Galaxy Cluster

This project analyzes spectroscopic redshift data of galaxies in a given field to estimate the **dynamical mass** of a galaxy cluster. The analysis includes data cleaning, visualization, velocity calculation, and mass estimation using astrophysical principles and observational data.

## 📁 Files Used
- `Skyserver_SQL6_22_2025 6_29_58 AM.csv`: Raw catalog of galaxy data retrieved from SDSS SkyServer.
- `code.py`: Main Python script that performs the full analysis pipeline.

## 🔧 Key Tools & Libraries
- `pandas` for data manipulation  
- `matplotlib.pyplot` for plotting  
- `numpy` for numerical operations  
- `astropy.constants` and `astropy.units` for physical constants  
- `astropy.cosmology` for cosmological model (Planck18)

## 🔬 Main Steps

### 1. Data Preparation
- Read the CSV data.
- Aggregate by `objid` to remove duplicates and get a single entry per galaxy.

### 2. Redshift Analysis
- Compute the **mean** and **standard deviation** of spectroscopic redshift (`specz`).
- Apply a 3-sigma filter to isolate galaxies likely belonging to a cluster.

### 3. Visualizations
- **Boxplot** and **histogram** of redshift distribution.
- Histogram of velocity and angular separation.

### 4. Velocity Calculations
- Convert redshift to expansion velocity using the relativistic formula.
- Calculate **velocity dispersion**.

### 5. Spatial Analysis
- Compute average angular separation from the cluster center.

### 6. Mass Estimation
- Estimate physical diameter from redshift and separation.
- Apply virial theorem to estimate **dynamical mass** of the cluster.

## 📊 Results
- **Cluster Redshift**: Computed as the mean of filtered redshift values.
- **Velocity Dispersion**: Derived from relativistic redshift velocities.
- **Dynamical Mass**: Estimated using velocity dispersion and physical diameter.
