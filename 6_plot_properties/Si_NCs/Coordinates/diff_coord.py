import numpy as np

def find_coord_diff(file_1, file_2):
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
    diff = []
    diff_indices = []
    left_indices = []
    for i in range(coord1.shape[0]):
        flag = True
        for j in range(coord2.shape[0]):
            if coord1[i,0]==coord2[j,0] and coord1[i,1]==coord2[j,1] and coord1[i,2]==coord2[j,2]:
                flag = False
                break
        if flag:
            diff.append(lines1[2+i])
            diff_indices.append(i+1)
        else:
            left_indices.append(i+1)
    return len(diff), diff, diff_indices, left_indices
    

names = ['Si59H60', 'Si123H100', 'Si265H140', 'Si329H172', 'Si501H228', 'Si1009H412']
for name in names:
    natoms, diff, diff_indices, left_indices = find_coord_diff(F'{name}.xyz',F'{name}-core.xyz')
    file = open(F'{name}-surface.xyz','w')
    file.write(F'{natoms}\n\n')
    for i in range(len(diff)):
        file.write(diff[i])
    file.close()
    print(name, '\n', diff_indices, '\n--------\n', left_indices, '\n--------')

