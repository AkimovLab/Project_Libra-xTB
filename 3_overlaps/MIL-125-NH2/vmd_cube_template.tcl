mol load cube job1/Diag_992-WFN_04608_1-1_0.cube
display projection Orthographic
mol modstyle 0 0 CPK 1.200000 0.300000 15.000000 12.000000
mol representation CPK 1.200000 0.300000 15.000000 12.000000

mol showperiodic 0 0 x
mol numperiodic 0 0 1
mol showperiodic 0 0 xX
mol numperiodic 0 0 1
mol showperiodic 0 0 xyX
mol numperiodic 0 0 1
mol showperiodic 0 0 xyXY
mol numperiodic 0 0 1
mol showperiodic 0 0 xyXYZ
mol numperiodic 0 0 1
mol showperiodic 0 0 xyzXYZ
mol numperiodic 0 0 1

mol addrep 0
mol modstyle 1 0 Isosurface 0.01 0 0 0 1 1
mol modmaterial 1 0 Opaque
mol modcolor 1 0 ColorID 0

mol showperiodic 0 1 x
mol numperiodic 0 1 1
mol showperiodic 0 1 xY
mol numperiodic 0 1 1
mol showperiodic 0 1 xXY
mol numperiodic 0 1 1
mol showperiodic 0 1 xyXY
mol numperiodic 0 1 1
mol showperiodic 0 1 xyzXY
mol numperiodic 0 1 1
mol showperiodic 0 1 xyzXYZ
mol numperiodic 0 1 1

mol addrep 0
mol modcolor 2 0 ColorID 1
mol modstyle 2 0 Isosurface -0.01 0 0 0 1 1

mol showperiodic 0 2 x
mol numperiodic 0 2 1
mol showperiodic 0 2 xX
mol numperiodic 0 2 1
mol showperiodic 0 2 xXY
mol numperiodic 0 2 1
mol showperiodic 0 2 xXYZ
mol numperiodic 0 2 1
mol showperiodic 0 2 xyXYZ
mol numperiodic 0 2 1
mol showperiodic 0 2 xyzXYZ
mol numperiodic 0 2 1

color Display Background white

#rotate x by 225.0
#scale by 1.6
axes location off

scale by 0.833000
scale by 0.833000
scale by 0.833000
scale by 0.833000
#scale by 1.200000

render Tachyon Diag_992-WFN_04608_1-1_0 "/util/academic/vmd/1.9.2/lib/vmd/tachyon_LINUXAMD64 -aasamples 12 %s -format TGA -res 2048 2048 -o %s.tga"
