#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import numpy as np

# First find the coordinate of the file_2 indices in file_1
# The second argument must contain a subset of the coordinates from the first file
def find_coord_indices(file_1, file_2):
    # File_1 is larger than File_2
    # Coord 1
    file = open(file_1, 'r')
    lines1 = file.readlines()
    file.close()
    # Coord 2
    file = open(file_2, 'r')
    lines2 = file.readlines()
    file.close()
    # Difference
    coord1 = []
    for i in range(2,len(lines1)):
        tmp1 = lines1[i].split()
        coord1.append([float(tmp1[1]),float(tmp1[2]),float(tmp1[3])])
    coord1 = np.array(coord1)
    coord2 = []
    for i in range(2,len(lines2)):
        tmp1 = lines2[i].split()
        coord2.append([float(tmp1[1]),float(tmp1[2]),float(tmp1[3])])
    coord2 = np.array(coord2)
    indices = []

    for i in range(coord1.shape[0]):
        #flag = False
        for j in range(coord2.shape[0]):
            if coord1[i,0]==coord2[j,0] and coord1[i,1]==coord2[j,1] and coord1[i,2]==coord2[j,2]:
                indices.append(i)
        #if flag:
        #    indices.append(lines1[i])
    return indices


# In[64]:

# Compute the standard deviation of the position in the MD trajectory
def std_pos_v2(trajectory_filename, coords_main, sample_coords_filename, istep, fstep):
    
    file = open(trajectory_filename,'r')
    lines = file.readlines()
    file.close()
    
    natoms = int(lines[0].split()[0])
    md_coord = []
    for i in range(istep,fstep):
        tmp1 = []
        start_line = 2+i*(natoms+2)
        end_line = start_line+natoms
        for j in range(start_line, end_line):
            tmp2 = lines[j].split()
            #print(tmp1)
            tmp1.append([float(tmp2[1]),float(tmp2[2]),float(tmp2[3])])
        md_coord.append(tmp1)
    md_coord = np.array(md_coord)
    
    file = open('sample.xyz', 'w')
    for i in range(natoms+2):
        file.write(lines[i])
    file.close()
    
    #_, diff, diff_indices, left_indices = find_coord_diff('sample.xyz', sample_coords_filename)
    indices = find_coord_indices(coords_main, sample_coords_filename)
    
    #print(indices, len(indices))
    stds = []
    for j in indices:
        ave_1 = []
        for i in range(len(md_coord)):
            ave_1.append(md_coord[i][j,:])
        ave_1 = np.array(ave_1)
        ave_1_val = np.average(ave_1,axis=0)
        std_1_val = np.std(ave_1,axis=0)
        stds.append(std_1_val)
    stds = np.average(np.array(stds),axis=0)
    average_all = np.average(stds)
    return stds, average_all


print('------------------------------------------------------ Si Core atoms')

print(std_pos_v2('../../3_overlaps/Si1009H412/Si1009H412-pos-1.xyz','Coordinates/Si1009H412.xyz',
                 'Coordinates/Si1009H412-core.xyz',0, 4000))
print(std_pos_v2('../../3_overlaps/Si501H228/Si501H228-pos-1.xyz','Coordinates/Si501H228.xyz',
                 'Coordinates/Si501H228-core.xyz', 2000, 6000))
print(std_pos_v2('../../3_overlaps/Si329H172/Si329H172-pos-1.xyz','Coordinates/Si329H172.xyz',
                      'Coordinates/Si329H172-core.xyz',2000, 6000))
print(std_pos_v2('../../3_overlaps/Si265H140/Si265H140-pos-1.xyz',
                 'Coordinates/Si265H140.xyz','Coordinates/Si265H140-core.xyz',2000, 6000))
print(std_pos_v2('../../3_overlaps/Si123H100/Si123H100-pos-1.xyz','Coordinates/Si123H100.xyz',
                      'Coordinates/Si123H100-core.xyz',2000, 6000))
print(std_pos_v2('../../3_overlaps/Si59H60/Si59H60-pos-1.xyz','Coordinates/Si59H60.xyz','Coordinates/Si59H60-core.xyz',2000, 6000))


# In[65]:

print('------------------------------------------------------ Si surface atoms')

print(std_pos_v2('../../3_overlaps/Si1009H412/Si1009H412-pos-1.xyz','Coordinates/Si1009H412.xyz',
                 'Coordinates/Si1009H412-surface-Si-only.xyz', 0, 4000))
print(std_pos_v2('../../3_overlaps/Si501H228/Si501H228-pos-1.xyz','Coordinates/Si501H228.xyz',
                 'Coordinates/Si501H228-surface-Si-only.xyz',2000, 6000))
print(std_pos_v2('../../3_overlaps/Si329H172/Si329H172-pos-1.xyz',
                      'Coordinates/Si329H172.xyz','Coordinates/Si329H172-surface-Si-only.xyz',2000, 6000))
print(std_pos_v2('../../3_overlaps/Si265H140/Si265H140-pos-1.xyz',
                      'Coordinates/Si265H140.xyz','Coordinates/Si265H140-surface-Si-only.xyz',2000, 6000))
print(std_pos_v2('../../3_overlaps/Si123H100/Si123H100-pos-1.xyz',
                      'Coordinates/Si123H100.xyz','Coordinates/Si123H100-surface-Si-only.xyz',2000, 6000))
print(std_pos_v2('../../3_overlaps/Si59H60/Si59H60-pos-1.xyz','Coordinates/Si59H60.xyz',
                 'Coordinates/Si59H60-surface-Si-only.xyz',2000, 6000))


# In[66]:

print('------------------------------------------------------ H atoms')

print(std_pos_v2('../../3_overlaps/Si1009H412/Si1009H412-pos-1.xyz','Coordinates/Si1009H412.xyz',
                 'Coordinates/Si1009H412-surface-H-only.xyz', 0, 4000))
print(std_pos_v2('../../3_overlaps/Si501H228/Si501H228-pos-1.xyz','Coordinates/Si501H228.xyz',
                 'Coordinates/Si501H228-surface-H-only.xyz',2000, 6000))
print(std_pos_v2('../../3_overlaps/Si329H172/Si329H172-pos-1.xyz',
                      'Coordinates/Si329H172.xyz','Coordinates/Si329H172-surface-H-only.xyz',2000, 6000))
print(std_pos_v2('../../3_overlaps/Si265H140/Si265H140-pos-1.xyz',
                      'Coordinates/Si265H140.xyz','Coordinates/Si265H140-surface-H-only.xyz',2000, 6000))
print(std_pos_v2('../../3_overlaps/Si123H100/Si123H100-pos-1.xyz',
                      'Coordinates/Si123H100.xyz','Coordinates/Si123H100-surface-H-only.xyz',2000, 6000))
print(std_pos_v2('../../3_overlaps/Si59H60/Si59H60-pos-1.xyz','Coordinates/Si59H60.xyz',
                 'Coordinates/Si59H60-surface-H-only.xyz',2000, 6000))


