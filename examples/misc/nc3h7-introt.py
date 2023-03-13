#! /usr/bin/env python3

"""
Sample calculation of reduced rotational constants of internal rotors

output:
CH3CH2CH2: A = 1.10189 [cm^-1]
CH3CH2CH2: B = 0.300552 [cm^-1]
CH3CH2CH2: C = 0.259693 [cm^-1]
[CH3CH2]-[CH2] rotor: reduced B = 10.5666 [cm^-1]
[CH3]-[CH2CH2] rotor: reduced B = 6.40172 [cm^-1]
"""

import me2d

xyz_geom = """10
n-C3H7, wb97xd/6-311++g(d,p) int=ultrafine
 C   -1.1744721654  -0.4088328261  -0.1298116247
 C    0.0009611900   0.5597098446  -0.0037821686
 C    1.3288170792  -0.1097868627  -0.0224493804
 H   -2.1299994802   0.1200445188  -0.0927586721
 H   -1.1641727016  -1.1404181191   0.6835112384
 H   -1.1277496421  -0.9580928940  -1.0740179565
 H   -0.1084158009   1.1367619043   0.9296230626
 H   -0.0455173692   1.3107287716  -0.8031276804
 H    2.2382175335   0.4631304819  -0.1606794360
 H    1.4223267768  -1.1639367593   0.2152427076
"""

atoms, coord = me2d.read_geom(xyz_geom)
A, B, C = me2d.external_rotc(atoms, coord)
#print(me2d.suggest_rotatable_bonds(atoms, coord))
B23 = me2d.internal_rotc(atoms, coord, 2, 3) # rotor around C2-C3 axis
B12 = me2d.internal_rotc(atoms, coord, 1, 2) # rotor around C1-C2 axis
print("CH3CH2CH2: A = %g [cm^-1]" % (A))
print("CH3CH2CH2: B = %g [cm^-1]" % (B))
print("CH3CH2CH2: C = %g [cm^-1]" % (C))
print("[CH3CH2]-[CH2] rotor: reduced B = %g [cm^-1]" % (B23))
print("[CH3]-[CH2CH2] rotor: reduced B = %g [cm^-1]" % (B12))
