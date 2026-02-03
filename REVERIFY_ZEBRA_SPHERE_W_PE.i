c MCNP Input of VERUS Sheilding problem
c Point source simulation with multiple shielding layers
c Based on working Co60_gammas_electrons.i format
c
c === Cell Cards ===
c Universe Cell (containing everything and excluding some areas)
100    0    50:999                        IMP:n,p,e=0    $ Outside World - Zero Importance
200    0   -50 20 -2                      IMP:n,p,e=1    $ Air gap (6 ft) with tracking
c Source Region                                  
900    0   -50 -20                        IMP:n,p,e=1    $ Source region - air around point source
c Concrete wall bins (6 ft = 182.88 cm thick)   
1000   2  -0.93   -50 2 -3               IMP:n,p,e=1    $ Bin 1
1100   1  -19.3   -50 3 -4               IMP:n,p,e=1    $ Bin 2
1200   2  -0.93   -50 4 -5               IMP:n,p,e=1    $ Bin 3
1300   1  -19.3   -50 5 -6               IMP:n,p,e=1    $ Bin 4
1400   2  -0.93   -50 6 -7               IMP:n,p,e=1    $ Bin 5
1500   1  -19.3   -50 7 -8               IMP:n,p,e=1    $ Bin 6
1600   2  -0.93   -50 8 -9               IMP:n,p,e=1    $ Bin 7
1700   1  -19.3   -50 9 -10              IMP:n,p,e=1    $ Bin 8
1800   2  -0.93   -50 10 -11             IMP:n,p,e=1    $ Bin 9
1900   1  -19.3   -50 11 -12             IMP:n,p,e=1    $ Bin 10
c Past the wall Air region                
2000   0          -50 12                 IMP:n,p,e=1    $ Finish Line Bin
c

c === Surface Cards ===
c Outer universe boundary
999  so   400                   $ Universe boundary sphere
c Bounding Box Surface
50   SPH 0 0 0  300             $ Rectangular box for geometry containment
c Material boundaries (converted from mm to cm)
2    SPH 0 0 0  182.88          $ 6 ft from origin (start of wall)
3    SPH 0 0 0  190.88          $ Bin 1 edge  (wall)
4    SPH 0 0 0  198.88          $ Bin 2 edge  (wall)
5    SPH 0 0 0  206.88          $ Bin 3 edge  (wall)
6    SPH 0 0 0  214.88          $ Bin 4 edge  (wall)
7    SPH 0 0 0  222.88          $ Bin 5 edge  (wall)
8    SPH 0 0 0  230.88          $ Bin 6 edge  (wall)
9    SPH 0 0 0  238.88          $ Bin 7 edge  (wall)
10   SPH 0 0 0  246.88          $ Bin 8 edge  (wall)
11   SPH 0 0 0  254.88          $ Bin 9 edge  (wall)
12   SPH 0 0 0  262.88          $ End of wall (wall)
20   SPH 0 0 0  1               $ Source card Dont Need

c === Material Cards ===
c W
m1  74000 -1.000000     $19.3g/cm3
c PE
m2  1001  -0.143686     $0.93g/cm3
    1002  -0.000033
	6000  -0.856276
c Tissue Equivalent MS20, (PNNL-15870 Rev 1)
m3  1001   -0.08119   $1g/cm3 (To 3 SIG-Figs)
    6000   -0.583442
    7014   -0.017728
    7015   -0.000069
    8000   -0.186381
    12000  -0.130287
    17035  -0.000673
    17037  -0.000227
c Ordinary concrete composition (PNNL-15870 Rev 1) LANL-MIX
c m1 1001   -0.004529
c    1002   -0.000001
c    8016   -0.511211
c    8017   -0.000207
c    8018   -0.001182
c    14000  -0.360360
c    13027  -0.035550
c    11023  -0.015270
c    20000  -0.057910
c    26000  -0.013780
c Water for quick doese calculation
c m2  1000   -0.111902    $0.997g/cm3
c    1002   -0.000026
c    8016   -0.885692
c    8017   -0.000359
c    8018   -0.002048
c Tissue Equivalent MS20, (PNNL-15870 Rev 1)
c m5  1001   -0.81171     $1g/cm3 (To 3 SIG-Figs)
c    1002   -0.000019
c    6000   -0.583442
c    7014   -0.017728
c    7015   -0.000069
c    8016   -0.185876
c    8017   -0.000075
c    8018   -0.000430
c    12000  -0.130287
c    17035  -0.000673
c    17037  -0.000227
c Air dry at sea level
c m2  6000   -0.000124
c     7014   -0.752316
c     7015   -0.002944
c     8016   -0.231153
c     8017   -0.000094
c     8018   -0.000535
c     18000  -0.012827
c === Mode ===
phys:p
mphys on
mode n p e
c mode n p e $ This is the full implementation of importances
c === Source Definition ===
c 14 MeV point source at origin
sdef pos=0 0 0 erg=14.0
c === Tallies ===
f1:n        2
c1  0 1
*f12:n      2 3 4 5 6 7 8 9 10 11 12    $ Energy flux in each bin (Neutrons)
e12 2.5e-5 0.001 0.511 100i 20.0 $ 50 Energy bins
c12 0 1
*f22:p      2 3 4 5 6 7 8 9 10 11 12  $ Energy flux in each bin (Photons)
e22 2.5e-5 0.001 0.511 100i 20.0 $ 50 Energy bins
c22 0 1
*f32:e      2 3 4 5 6 7 8 9 10 11 12  $ Energy flux in each bin (Electrons)
*f42:n,p,e  2 3 4 5 6 7 8 9 10 11 12  $ Flux in each bin
e42 2.5e-5 0.001 0.511 100i 20.0 $ 50 Energy bins
c --- F4 tally with ICRP 116 AP neutron dose (mrem) ---
c f14:n,p,e   1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000
c
c de14 0.01 0.015 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.1 0.15 0.2 &
c    0.3 0.4 0.5 0.511 0.6 0.662 0.8 1.0 1.117 1.33 1.5 2.0 &
c    3.0 4.0 5.0 6.0 6.129 8.0 10.0 15.0 20.0
c
c df14 6.85E-09 1.56E-08 2.25E-08 3.13E-08 3.51E-08 3.70E-08 3.90E-08 4.13E-08 4.44E-08 5.19E-08 7.48E-08 1.00E-07 &
c     1.51E-07 2.00E-07 2.47E-07 2.52E-07 2.91E-07 3.17E-07 3.73E-07 4.49E-07 4.90E-07 5.59E-07 6.12E-07 7.48E-07 &
c	 9.75E-07 1.17E-06 1.34E-06 1.50E-06 1.51E-06 1.78E-06 2.05E-06 2.61E-06 3.08E-06
c
c === Run Control ===
nps  5e6
prdmp  0 0 1 0     $ Report and dump every 1e9 particles
print  10 110 140 170