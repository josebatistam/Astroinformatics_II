'''
Data file: abellN1989.dat

Description: This file contains a catalog of galaxy clusters from Abell, Corwin,
             and Olowin (1989).
Columns:
- GLON and GLAT: the angular coordinates of each galaxy cluster in galactic
                 coordinates (l, b). These coordinates are the observed angular
                 position of the cluster, as viewed from the Sun (i.e., our
                 viewing position). If, for example, a cluster has coordinates
                 (0°, 20°), it would be viewed 20 degrees above the Galactic
                 disk plane, in the direction of the Galactic center.
- Rich: richness, as from Abell, ApJS, 3, 211-288 (1958). It is a number related
        to the galaxy counts, between 0 (30-40 galaxies) and 5 (over 300
        galaxies).
- Dclass: distance class, as from Abell, ApJS, 3, 211-288 (1958). It is a number
          related to the distance, between 1 (near) and 7 (extremely distant).
- m10: red magnitude of the tenth brightest cluster member, as from Abell, ApJS,
       3, 211-288 (1958).
- LUMDIST: luminosity distance D_L (Mpc), calculated with the ΛCDM standard
           cosmology (H_0=70 km/s/Mpc, Ω_M=0.3, Ω_Λ=0.7). Note: if LUMDIST is
           -99.0, that means there is no estimate.
'''

import numpy as np, matplotlib.pyplot as plt, pandas as pd, timeit

# Set the path to the data file
file = 'Practices/Practice_1/abellN1989.dat'

def main():
    # Create a DataFrame by reading the data file, specifying column positions (fix)
    df = pd.read_csv(file, sep='\s+', usecols = range(0, 6),
                     names = ["GLON", "GLAT", "Rich", "Dclass", "m10", "LUMDIST"],
                     header = 0)
    # Task 1: Access Data and Plot
    # Wrap longitudes to [-180, 180] and convert cordinates to radians
    lon = np.radians((df['GLON'].values + 180) % 360 - 180)
    lat = np.radians(df['GLAT'].values)
    # Plot the data using Hammer projection
    plt.figure(figsize=(10, 5))
    ax = plt.subplot(111, projection = "hammer")
    ax.scatter(lon, lat, s = 10, color = "#0080ff", edgecolors = "none",
        alpha = 0.5)
    ax.set_title("AbellN1989 Galaxy Clusters Positons on the Sky", pad = 25)
    ax.set_xlabel("Galactic Longitude [degrees]", labelpad = 25)
    ax.set_ylabel("Galactic Latitude [degrees]")
    # Force ticks at 30 and 15 degree intervals for x and y axes, respectively
    x_ticks = np.arange(-180, 181, 30)
    y_ticks = np.arange(-90, 91, 15)
    ax.set_xticks(np.radians(x_ticks))
    ax.set_xticklabels([f"{x_tick}°" for x_tick in x_ticks])
    ax.set_yticks(np.radians(y_ticks))
    ax.set_yticklabels([f"{y_tick}°" for y_tick in y_ticks])
    # Add grid lines
    ax.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig('Practices/Practice_1/clusters_all_sky.pdf', format='pdf')
    plt.close()

    # Task 2: Filter Data and Plot

    # Transfrom galactic coordinates to Cartesian and add them to the DataFrame
    df['x'] = df['LUMDIST'].values * np.cos(lon) * np.cos(lat)
    df['y'] = df['LUMDIST'].values * np.sin(lon) * np.cos(lat)
    df['z'] = df['LUMDIST'].values * np.sin(lat)
    # Filter out clusters without measured luminosity distance (LUMDIST != -99)
    df_mld = df[df['LUMDIST'] != -99.0]
    # Item A: Plot the X-Y and X-Z projection of the filtered clusters
    fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(10,5))
    fig.suptitle("AbellN1989 Galaxy Clusters Positons on the Sky")
    ax1.scatter(df_mld['x'], df_mld['y'], s = 10, color = "#ff8000",
        edgecolors = "none", alpha = 0.5)
    ax1.set_title("X-Y Projection")
    ax1.set_xlabel("X [Mpc]")
    ax1.set_ylabel("Y [Mpc]")
    ax1.grid(True, linestyle='--', alpha=0.5)
    ax2.scatter(df_mld['x'], df_mld['z'], s = 10, color = "#0080ff",
        edgecolors = "none", alpha = 0.5)
    ax2.set_title("X-Z Projection")
    ax2.set_xlabel("X [Mpc]")
    ax2.set_ylabel("Z [Mpc]")
    ax2.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig('Practices/Practice_1/clusters_xy_xz_all.pdf', format='pdf')
    plt.close()
    # Filter out faint clusters (m10 < 17.0)
    df_bright = df_mld[df_mld['m10'] < 17.0]
    # Filter out bright clusters (m10 >= 17.0)
    df_faint = df_mld[df_mld['m10'] >= 17.0]
    # Item B: Plot the X-Y and X-Z projection of the further-filtered clusters
    # Bright clusters (m10 < 17.0)
    fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(10,5))
    fig.suptitle("AbellN1989 Galaxy Clusters Positons on the Sky (m10 < 17.0)")
    ax1 = plt.subplot(121)
    ax1.scatter(df_bright['x'], df_bright['y'], s = 10, color = "#ff8000",
        edgecolors = "none", alpha = 0.5)
    ax1.set_title("X-Y Projection")
    ax1.set_xlabel("X [Mpc]")
    ax1.set_ylabel("Y [Mpc]")
    ax1.grid(True, linestyle='--', alpha=0.5)
    ax2.scatter(df_bright['x'], df_bright['z'], s = 10, color = "#0080ff",
        edgecolors = "none", alpha = 0.5)
    ax2.set_title("X-Z Projection")
    ax2.set_xlabel("X [Mpc]")
    ax2.set_ylabel("Z [Mpc]")
    ax2.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig('Practices/Practice_1/clusters_xy_xz_bright.pdf', format='pdf')
    plt.close()
    # Faint clusters (m10 < 17.0)
    fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(10,5))
    fig.suptitle("AbellN1989 Galaxy Clusters Positons on the Sky (m10 >= 17.0)")
    ax1.scatter(df_faint['x'], df_faint['y'], s = 10, color = "#ff8000",
        edgecolors = "none", alpha = 0.5)
    ax1.set_title("X-Y Projection")
    ax1.set_xlabel("X [Mpc]")
    ax1.set_ylabel("Y [Mpc]")
    ax1.grid(True, linestyle='--', alpha=0.5)
    ax2.scatter(df_faint['x'], df_faint['z'], s = 10, color = "#0080ff",
        edgecolors = "none", alpha = 0.5)
    ax2.set_title("X-Z Projection")
    ax2.set_xlabel("X [Mpc]")
    ax2.set_ylabel("Z [Mpc]")
    ax2.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig('Practices/Practice_1/clusters_xy_xz_faint.pdf', format='pdf')
    plt.close()
    # Item C: Compare the distributions of m10 for the two subsets
    fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(10,5))
    fig.suptitle("AbellN1989 Bright and Faint Galaxy Clusters")
    ax1.scatter(df_bright['x'], df_bright['y'], s = 10, color = "#ff8000",
        edgecolors = "none", alpha = 0.5, label = "Bright (m10 < 17.0)")
    ax1.scatter(df_faint['x'], df_faint['y'], s = 10, color = "#0080ff",
        edgecolors = "none", alpha = 0.5, label = "Faint (m10 < 17.0)")
    ax1.set_title("X-Y Projection")
    ax1.set_xlabel("X [Mpc]")
    ax1.set_ylabel("Y [Mpc]")
    ax1.grid(True, linestyle='--', alpha=0.5)
    ax1.legend()
    ax2.scatter(df_bright['x'], df_bright['z'], s = 10, color = "#ff8000",
        edgecolors = "none", alpha = 0.5, label = "Bright (m10 < 17.0)")
    ax2.scatter(df_faint['x'], df_faint['z'], s = 10, color = "#0080ff",
        edgecolors = "none", alpha = 0.5, label = "Faint (m10 < 17.0)")
    ax2.set_title("X-Z Projection")
    ax2.set_xlabel("X [Mpc]")
    ax2.set_ylabel("Z [Mpc]")
    ax2.grid(True, linestyle='--', alpha=0.5)
    ax2.legend()
    plt.tight_layout()
    plt.savefig('Practices/Practice_1/clusters_xy_xz_both.pdf', format='pdf')
    plt.close()

    # Task 3: Plot More Information

    # Plot richness against against luminonsity distance
    plt.figure(figsize=(10, 6))
    plt.scatter(df_faint['LUMDIST'], df_faint['Rich'], s = 20, color = "#0080ff",
        edgecolors = "none", alpha = 0.5)
    plt.title("AbellN1989 Galaxy Clusters Richness vs. Luminosity Distance")
    plt.xlabel("Luminosity Distance [Mpc]")
    plt.ylabel("Richness")
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig('Practices/Practice_1/clusters_richness_lumdist.pdf', format='pdf')
    plt.close()

if __name__ == '__main__':
    # Measure execution: 5 repetitions, 10 loops per repetition
    times = timeit.repeat("main()", globals=globals(), repeat=5, number=10)
    print(f"Times for 5 repetitions of 10 loops each: {times}")
    print(f"Best average per run: {min(times)/10:.3f} s")

# Task 4: Code Optimization

'''
With your code written so far, try to find parts of the code that could be sped
up with the techniques learned in lecture 2 and 3. Measure code execution time.
You might change your code to the usage of more Pythonic code such as generators.
Describe which part can be speed up, and try to implementing it and again measure
the code execution time.
'''

# 1 - Separating computations and plotting
# 2 - Serializing the data (pickle)
# 3 - Being more careful with specifying your datatypes
# 4 - Use list comprehensions (with Conditional Filtering when applicable)
# 5 - Use file handling (Context Managers)
# 6 - Use Generator Functions (lazy loading)
# 7 - Use Caching (saving results to files) (functools module when applicable)
# 8 - Using NumPy arrays instead of lists where possible
# 9 - Using NumPy arrays instead of Pandas DataFrame where possible
# 10 - Using profilers to monitor the memory usage of your program
# 11 - Using built-in functions and libraries when more efficient than 3rd party
# 12 - Using local instead of global variables
# 13 - Avoiding the use of recursion. Using iteration instead
# 14 - Using vectorized operations and broadcasting when performing calculations
# 15 - Using multi-processing, multi-threading, async I/O and GPU acceleration