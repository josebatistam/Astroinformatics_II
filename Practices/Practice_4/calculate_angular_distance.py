from calculate_angular_distance_helper import calculate_angular_distance_caller

def read_file(filename):
    ras_deg = []
    decs_deg = []
    zs = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#'): # Skip header
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            ra_deg = float(parts[0])
            dec_deg = float(parts[1])
            z_val = float(parts[2])
            ras_deg.append(ra_deg)
            decs_deg.append(dec_deg)
            zs.append(z_val)
    return ras_deg, decs_deg, zs

def find_closest_pair(ras_deg, decs_deg):
    n = len(ras_deg)
    if n < 2:
        raise ValueError('At least two galaxies are needed to calculate the angular distance')
    min_dist_deg = 180.0
    min_i = 0
    min_j = 1
    # Calculate for all possible galaxy pairs
    for i in range(n):
        ra_i = ras_deg[i]
        dec_i = decs_deg[i]
        for j in range(i + 1, n):
            ra_j = ras_deg[j]
            dec_j = decs_deg[j]
            dist_deg = calculate_angular_distance_caller(ra_i, dec_i, ra_j, dec_j)
            if dist_deg < min_dist_deg:
                min_dist_deg = dist_deg
                min_i = i
                min_j = j
    return min_i, min_j, min_dist_deg

filename = 'Practices/Practice_4/coma_members_Jim2025.dat'
ras_deg, decs_deg, zs = read_file(filename = filename)

print(f'Read {len(ras_deg)} galaxies from {filename}')

min_i, min_j, min_dist_deg = find_closest_pair(ras_deg, decs_deg)

print('\nClosest pair of galaxies:')
print(f'- Galaxy #{min_i + 1}: RA = {ras_deg[min_i]:.6f} deg, DEC = {decs_deg[min_i]:.6f} deg, z = {zs[min_i]:.6f}')
print(f'- Galaxy #{min_j + 1}: RA = {ras_deg[min_j]:.6f} deg, DEC = {decs_deg[min_j]:.6f} deg, z = {zs[min_j]:.6f}')
print(f'Minimum angular distance = {min_dist_deg:.10f} deg ({min_dist_deg*3600:.5f} arcsec)')