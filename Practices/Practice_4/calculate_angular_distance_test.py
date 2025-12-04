from calculate_angular_distance_helper import calculate_angular_distance_caller

# Test with first two 'coma_members_Jim2025.dat' values
ra1_deg, dec1_deg = 195.48264888888886, 29.323511111111106
ra2_deg, dec2_deg = 193.7912769444444, 26.00520111111111

dist = calculate_angular_distance_caller(ra1_deg, dec1_deg, ra2_deg, dec2_deg)
print(f'Angular distance between ({ra1_deg:.5f}, {dec1_deg:.5f}) and '
      f'({ra2_deg:.5f}, {dec2_deg:.5f}) = {dist:.5f} deg')

# Sanity check: distance of a point with itself should be 0
dist_zero = calculate_angular_distance_caller(ra1_deg, dec1_deg, ra1_deg, dec1_deg)
print(f'Self-distance = {dist_zero:.2f} deg')