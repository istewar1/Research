Title: Simple geometry test (shield, reflector, soil, air-detector)
c CELLS
1 2 -0.001   -1 $ inner cylinder
2 1 -3.667 1 -2 $ Cut-out NaI detector
3 3 -11.25 1 2 -3 $ Lead shield surrounding NaI
4 4 -2.30 -4
5 0 1 2 3 4 -5
99 0 5
c

c SURFACES
1 rcc 0 0 2.54  0 0 27.94 7.62  $ Inner cylinder (open)
2 rcc 0 0 0     0 0 30.48 12.7  $ Cut-out NaI detector
3 rcc 0 0 -2.54 0 0 35.56 15.24 $ Lead cylinder outside detector
4 rpp -150 -115 -50 50 -70 130 $ soil
5 rpp -160 20 -60 60 -80 140     $ Surrounding atmosphere

c
imp:p 1 4r 0 
c NaI(Tl) Detector
m1 81000.01p 0.001 11000.01p 0.499 53000.01p 0.5 
c Atmosphere (air)
m2 7014 -0.01197
     7015 -4.68E-05
     8016 -0.98531
     8017 -0.00039889
     8018 -0.0022785
c Lead (shielding)
m3 82204 -0.013781 82206 -0.23956 82207 -0.22074 82208 -0.52592
c soil
m4 8016 -0.87477
     8017 -0.00035414
     8018 -0.0020229
     11023 -0.0013145
     12024 -0.0021001
     12025 -0.00027696
     12026 -0.0003171
     13027 -0.012504
     14028 -0.043661
     14029 -0.0022973
     14030 -0.0015683
     19039 -0.0016762
     19040 -2.16E-07
     19041 -0.00012717
     20040 -0.0060741
     20042 -4.26E-05
     20043 -9.09E-06
     20044 -0.00014376
     20046 -2.88E-07
     20048 -1.41E-05
     22046 -3.75E-05
     22047 -3.46E-05
     22048 -0.00034962
     22049 -2.62E-05
     22050 -2.56E-05
     25055 -0.00064102
     26054 -0.0028007
     26056 -0.045591
     26057 -0.0010717
     26058 -0.00014513
c
C ----------------------------------------
C      Source Cards:  Soil
C ----------------------------------------
 sdef  par=2 x=d1 y=d2 z=d3  erg=d4 wgt=1.0
c
  sc1  x coordinate
  si1  -150 -115
  sp1  0.0   1.0
c
  sc2  y coordinate
  si2  -50   50
  sp2  0.0   1.0
c
  sc3  z coordinate
  si3  -70   130
  sp3  0.0   1.0
c
c
 read file=soilDist.sdef  noEcho
c
C      -------------------------------------------------------------------------
C      Tallies
C          multiplier is product of
C          soil multiplier is product of 
C               309.57 gammas/kg
C               0.00152 kg/cm3
C               700,000 cm^3
c               = 329,382
C      -------------------------------------------------------------------------
f4:p  2
f4m 329382
e4  0.010 1023i 4.000
ft4 geb  -0.014 0.11 0
c
mode p
nps 1E9
