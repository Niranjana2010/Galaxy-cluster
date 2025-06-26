# Importing Necessary Libraries
# We begin by importing Python libraries commonly used in data analysis and visualization:
# - `numpy` for numerical operations
# - `matplotlib.pyplot` for plotting graphs
# - `pandas` for handling CSV data, which is especially useful for tabular data such as redshift catalogs

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.constants import G, c
from astropy.cosmology import Planck18 as cosmo
import astropy.units as u

# ### step 1:  Before we begin calculations, we define key physical constants used throughout.
# Constants
#H_0 = Hubble constant in SI
#c = Speed of light in m/s
#G = Gravitational constant 
q0 = -0.534  # Deceleration parameter (assumed from Planck fit KEEP it as it is)


# Read the csv data into the python using the method below
df = pd.read_csv('Skyserver_SQL6_22_2025 6_29_58 AM.csv') 


# Calculating the Average Spectroscopic Redshift (`specz`) for Each Object
 
# When working with astronomical catalogs, an object (identified by a unique `objid`) might have multiple entries — for example, due to repeated observations. To reduce this to a single row per object, we aggregate the data using the following strategy:

averaged_df = df.groupby('objid').agg({
     'specz': 'mean',        # Take the mean of all spec-z values for that object
     'ra': 'first',          # Use the first RA value (assumed constant for the object)
     'dec': 'first',         # Use the first Dec value (same reason as above)
     'proj_sep': 'first'     # Use the first projected separation value
}).reset_index()

averaged_df.describe()['specz']


# To create a cut in the redshift so that a cluster can be identified. We must use some logic. Most astronomers prefer anything beyond 3*sigma away from the mean to be not part of the same group. 
#  Find the mean, standard deviation and limits of the redshift from the data

mean_val = averaged_df['specz'].mean()
std_dev = averaged_df['specz'].std()
lower_limit = mean_val - 3 * std_dev
upper_limit = mean_val + 3 * std_dev

print(f"Mean redshift: {mean_val:.6f}")
print(f"Standard deviation: {std_dev:.6f}")
print(f"Redshift range for cluster members: {lower_limit:.6f} to {upper_limit:.6f}")


#boxplot visualization of overall data
plt.figure(figsize=(8, 4))
plt.boxplot(df['specz'], vert=False, patch_artist=True,
            boxprops=dict(facecolor='lightgreen', color='black'),
            medianprops=dict(color='red', linewidth=2))
plt.title("Boxplot of Spectroscopic Redshift", fontsize=14)
plt.xlabel("Redshift (specz)", fontsize=12)
plt.grid(True, axis='x', linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()


# Plot the dsitribution of redshift as histogram and a boxplot 
# But the best plot would be a histogram to see where most of the objects downloaded lie in terms of redshift value

plt.figure(figsize=(10, 5))
plt.hist(df['specz'], bins=50, color='skyblue', edgecolor='black')
plt.title("Distribution of Spectroscopic Redshift", fontsize=14)
plt.xlabel("Redshift (specz)", fontsize=12)
plt.ylabel("Number of Galaxies", fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()


plt.figure(figsize=(10, 5))
plt.hist(averaged_df['specz'], bins=90, color='cornflowerblue', edgecolor='black', alpha=0.8)
plt.title("Distribution of Averaged Spectroscopic Redshift", fontsize=16)
plt.xlabel("Redshift (specz)", fontsize=13)
plt.ylabel("Number of Galaxies", fontsize=13)
plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
plt.tight_layout()
plt.show()


# Filter your data based on the 3-sigma limit of redshift. You should remove all data points which are 3-sigma away from mean of redshift
# Filtering the data based on specz values, used 3 sigma deviation from mean as upper limit.
filtered_df = averaged_df [( averaged_df['specz'] >= lower_limit) & ( averaged_df['specz'] <= upper_limit )].copy()
print(f"Number of galaxies in the cluster: {len(filtered_df)}")


# Use the relation between redshift and velocity to add a column named velocity in the data. This would tell the expansion velocity at that redshift 
filtered_df['velocity'] = c.value * filtered_df['specz'] 


#plot the velocity column created as hist
plt.figure(figsize=(10, 5))
plt.hist(filtered_df['velocity'], bins=50, color='tomato', edgecolor='black', alpha=0.85)
plt.title("Distribution of Galaxy Velocities", fontsize=16)
plt.xlabel("Velocity (km/s)", fontsize=13)
plt.ylabel("Number of Galaxies", fontsize=13)
plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
plt.tight_layout()
plt.show()


# ### Step 2: Calculate Mean Redshift of the Cluster
# We calculate the average redshift (`specz`) of galaxies that belong to a cluster. 
cluster_redshift = filtered_df['specz'].mean()

def relativistic_velocity(z, z_cluster):
    z1 = (1 + z)**2
    z2 = (1 + z_cluster)**2
    return c.value * (z1 - z2) / (z1 + z2)

filtered_df['v_dispersion'] = filtered_df['specz'].apply(lambda z: relativistic_velocity(z, cluster_redshift))
disp = filtered_df['v_dispersion'].std()

print(f"The value of the cluster redshift = {cluster_redshift:.4}")
print(f"The characteristic value of velocity dispersion of the cluster along the line of sight = {disp:.4} km/s.")
print(filtered_df['v_dispersion'].describe())


# ### Step 4: Visualizing Angular Separation of Galaxies
# We plot a histogram of the projected (angular) separation of galaxies from the cluster center. This helps us understand the spatial distribution of galaxies within the cluster field.
#Plot histogram for proj sep column
plt.hist(filtered_df['proj_sep'], bins=40, color='mediumseagreen', edgecolor='black')
plt.xlabel("Projected Angular Separation (arcminutes)")
plt.ylabel("Number of Galaxies")
plt.title("Angular Separation of Galaxies from Cluster Center")
plt.grid(True)
plt.show()


# ### Step 5: Estimating Physical Diameter of the Cluster
r = (c.value * cluster_redshift / cosmo.H0.value) * (1 - 0.5 * cluster_redshift * (1 + q0))
ra = r / (1 + cluster_redshift)

theta = filtered_df['proj_sep'].mean() * (np.pi / 180) / 60
diameter =  ra * theta 

# ### Step 6: Calculating the Dynamical Mass of the Cluster
# This mass estimate assumes the cluster is in dynamical equilibrium and bound by gravity.
### Calculating the dynamical mass in solar masses:
M_dyn =3*((disp*1000)**2)*(diameter*0.5*10**6*3*10**16)/(G*2*10**30)
print(f"Dynamical Mass of the cluster is {M_dyn:.2e} solar mass")
