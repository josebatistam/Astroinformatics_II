'''
Optimized version workflow:
Check for pickle:
- If exists:
  - Load data (already filtered, Cartesian, correct dtypes, etc.)
  - Call a general plotting function (with tick labels generated efficiently).
- If doesn't exist:
  - Load .dat file
  - Convert to NumPy arrays (since all pandas features are not needed)
  - Apply dtype optimization (int8, float32) 
  - Compute radians, Cartesian coords, NumPy masks for filtering (bright/faint) 
  - Bundle results into a dict 
  - Serialize with pickle
  - Call the plotting function.
'''

import os, pickle, numpy as np, matplotlib.pyplot as plt, pandas as pd, timeit

# Set the path to the data and pickle file
data_file = 'Practices/Practice_1/abellN1989.dat'
pickle_file = 'Practices/Practice_1/abellN1989_ready.pkl'

def plot_clusters(x, y, subtitle, filename, rect_proj=True):
    # Set number of subplots and their properties acording to len(x)
    if len(x) == 1:
        if rect_proj is True:
            fig, ax = plt.subplots(figsize=(10,5)) 
            title = 'AbellN1989 Galaxy Clusters Positons on the Sky'
            xlabel, ylabel = 'Galactic Longitude [degrees]', 'Galactic Latitude [degrees]'
        else:
            fig, ax = plt.subplots(figsize=(10,5), subplot_kw={'projection':'hammer'})
            title = 'AbellN1989 Galaxy Clusters Richness vs. Luminosity Distance'
            xlabel, ylabel = 'Luminosity Distance [Mpc]', 'Richness'
        ax.scatter(x, y, s=15, color='#0080ff', edgecolors='none', alpha=0.5)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(True, linestyle='--', alpha=0.5)
    else:
        fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(10,5))
        fig.suptitle(subtitle)
        if len(x) == 2:
            x1, x2, y1, y2 = x[0], x[1], y[0], y[1]
            ax1.scatter(x1, y1, s=10, color='#ff8000', edgecolors='none', alpha=0.5)
            ax2.scatter(x2, y2, s=10, color='#0080ff', edgecolors='none', alpha=0.5)
        else:
            x1, x2, x3, x4, y1, y2, y3, y4 = x[0], x[1], x[2], x[3], y[0], y[1], y[2], y[3]
            llabel, rlabel = 'Bright (m10 < 17.0)', 'Faint (m10 < 17.0)'
            ax1.scatter(x1, y1, s=10, color = '#ff8000', edgecolors='none', alpha=0.5, label=llabel)
            ax1.scatter(x3, y3, s=10, color = '#0080ff', edgecolors='none', alpha=0.5, label=rlabel)
            ax2.scatter(x2, y2, s=10, color = '#ff8000', edgecolors='none', alpha=0.5, label=llabel)
            ax2.scatter(x4, y4, s=10, color = '#0080ff', edgecolors='none', alpha=0.5, label=rlabel)
        ax1.set_title('X-Y Projection')
        ax1.set_xlabel('X [Mpc]')
        ax1.set_ylabel('Y [Mpc]')
        ax1.grid(True, linestyle='--', alpha=0.5)
        ax2.set_title('X-Z Projection')
        ax2.set_xlabel('X [Mpc]')
        ax2.set_ylabel('Z [Mpc]')
        ax2.grid(True, linestyle='--', alpha=0.5)
        ax1.legend()
        ax2.legend()
    plt.tight_layout()
    plt.savefig(filename, format='pdf')
    plt.close(fig)

def main():
    if os.path.exists(pickle_file):
        # Load preprocessed data
        with open(pickle_file, 'rb') as f:
             data = pickle.load(f)
        print('Loaded from pickle')
    else: # Create pickle file with the processed data
        print('Pickle not found. Processing data...')
        # Load raw data with explicit data types
        arr = np.genfromtxt(data_file, skip_header = 1,
                            dtype = [('GLON', np.float32), ('GLAT', np.float32),
                                    ('Rich', np.int8),('Dclass', np.int8),
                                    ('m10', np.float32),('LUMDIST', np.float32)])
        # Convert to radians
        lon = np.radians((arr['GLON'] + 180) % 360 - 180)
        lat = np.radians(arr['GLAT'])
        # Calculate cartesian coords
        mask = arr['LUMDIST'] != -99.0
        x = arr['LUMDIST'][mask] * np.cos(lon[mask]) * np.cos(lat[mask])
        y = arr['LUMDIST'][mask] * np.sin(lon[mask]) * np.cos(lat[mask])
        z = arr['LUMDIST'][mask] * np.sin(lat[mask])
        # Filter by lumdist and m10 values
        bright_clusters = arr['m10'][mask] < 17.0
        faint_clusters  = arr['m10'][mask] >= 17.0
        # Serialize
        data = {'lon': lon, 'lat': lat, 'x': x, 'y': y, 'z': z, 'mask': mask,
                'bright': bright_clusters, 'faint': faint_clusters,
                'rich': arr['Rich'], 'lumdist': arr['LUMDIST'],
                'm10': arr['m10']}
        with open(pickle_file, 'wb') as f:
             pickle.dump(data, f)
        print('Data processed and saved to pickle')
    # Create plots
    # 1. All-sky Hammer projection (lon vs lat)
    plot_clusters([data['lon']], [data['lat']], None,
                  'Practices/Practice_1/clusters_all_sky.pdf', False)
    # 2. Valid clusters: X-Y vs X-Z
    plot_clusters((data['x'], data['x']), (data['y'], data['z']),
                  'AbellN1989 Galaxy Clusters Positons on the Sky',
                  'Practices/Practice_1/clusters_xy_xz_all.pdf', True)
    # 3. Bright clusters: X-Y vs X-Z
    plot_clusters((data['x'][data['bright']], data['x'][data['bright']]),
                  (data['y'][data['bright']], data['z'][data['bright']]),
                  'AbellN1989 Galaxy Clusters Positons on the Sky (m10 < 17.0)',
                  'Practices/Practice_1/clusters_xy_xz_bright.pdf', True)
    # 4. Faint clusters: X-Y vs X-Z
    plot_clusters((data['x'][data['faint']], data['x'][data['faint']]),
                  (data['y'][data['faint']], data['z'][data['faint']]),
                  'AbellN1989 Galaxy Clusters Positons on the Sky (m10 >= 17.0)',
                  'Practices/Practice_1/clusters_xy_xz_faint.pdf', True)
    # 5. Bright vs Faint: both in two panels (with legends)
    plot_clusters((data['x'][data['bright']], data['x'][data['faint']],
                   data['x'][data['bright']], data['x'][data['faint']]),
                  (data['y'][data['bright']], data['y'][data['faint']],
                   data['z'][data['bright']], data['z'][data['faint']]),
                   'AbellN1989 Bright and Faint Galaxy Clusters',
                   'Practices/Practice_1/clusters_xy_xz_both.pdf', True)
    # 6. Richness vs Luminosity Distance (Cartesian)
    plot_clusters([data['lumdist'][data['mask']]], [data['rich'][data['mask']]],
                  None, 'Practices/Practice_1/clusters_richness_lumdist.pdf',
                  True)

if __name__ == '__main__':
    # Measure execution: 5 repetitions, 10 loops per repetition
    times = timeit.repeat("main()", globals=globals(), repeat=5, number=10)
    print(f"Times for 5 repetitions of 10 loops each: {times}")
    print(f"Best average per run: {min(times)/10:.3f} s")
