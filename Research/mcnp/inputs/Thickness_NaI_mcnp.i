Title: Simple geometry test (shield, reflector, soil, air-detector)
c CELLS
1 2 -0.001   -1 $ inner cylinder
2 1 -3.667 1 -2 $ Cut-out NaI detector
99 0 2
c

c SURFACES
1 rcc 0 0 2.54  0 0 27.94 7.62  $ Inner cylinder (open)
2 rcc 0 0 0     0 0 30.48 17.78 $ Cut-out NaI detector

c
imp:p 1 1r 0 
c NaI(Tl) Detector
m1 81000.01p 0.001 11000.01p 0.499 53000.01p 0.5 
c Atmosphere (air)
m2 7014 -0.01197
     7015 -4.68E-05
     8016 -0.98531
     8017 -0.00039889
     8018 -0.0022785
c
C ----------------------------------------
C      Source Cards:  Soil
C ----------------------------------------
 sdef  par=2 x=0 y=0 z=7.62  erg=d14 wgt=1.0
#         si14         #sp14
          L            D
          0.400		   1.0
fmesh4:p geom=xyz origin=0 -30 -30
           imesh=30 iints=100
           jmesh=30 jints=100
           kmesh=30 kints=100
f8:p  2
e8  0.010 1023i 3.000
c fm8 370
ft8 geb  -0.014 0.11 0
c
mode p
nps 1E7

