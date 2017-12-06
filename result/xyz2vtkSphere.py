
import math
import numpy as np
import os
import glob

runConfigFile = '../runConfig.txt'
with open(runConfigFile) as fp:
    for i, line in enumerate(fp):
        if i == 5:
            temp = np.fromstring(line, dtype='float64', sep=' ')
            box = np.zeros((2, 3))
            box[0, 0] = temp[0]
            box[0, 1] = temp[1]
            box[0, 2] = temp[2]
            box[1, 0] = temp[3]
            box[1, 1] = temp[4]
            box[1, 2] = temp[5]
            print(box)

f = open('box.vtk', 'w')
f.write('# vtk DataFile Version 3.0\n')
f.write('vtk file\n')
f.write('ASCII\n')
f.write('DATASET RECTILINEAR_GRID\n')
f.write('DIMENSIONS 2 2 2\n')
f.write('X_COORDINATES 2 float\n')
f.write(str(box[0,0])+' '+str(box[1,0])+'\n')
f.write('Y_COORDINATES 2 float\n')
f.write(str(box[0,1])+' '+str(box[1,1])+'\n')
f.write('Z_COORDINATES 2 float\n')
f.write(str(box[0,2])+' '+str(box[1,2])+'\n')
f.write('CELL_DATA 1\n')
f.write('POINT_DATA 8\n')




class Sphere:
    def __init__(self, linestring):
        data = linestring.split()
        self.gid = int(data[1])
        self.radius = float(data[2])
        self.pos = np.array([float(data[3]), float(data[4]), float(data[5])])


class Frame:
    def __init__(self, filename):
        self.SList = []
        file = open(filename, 'r')
        for line in file:
            if (line.startswith('S\t')):  # parse tubule data
                self.SList.append(Sphere(linestring=line))


def writeVTK(infile, filenamebase):

    vtkFileName = "SP"+filenamebase+'.vtk'
    if(os.path.isfile(vtkFileName) ):
    	print('File Exist')
   # 	return None

    frame=Frame(infile)
    # sort Tubule/Protein with gid
    frame.SList.sort(key=lambda S:S.gid)

    # write file header
    f = open(vtkFileName, 'w')
    f.write('# vtk DataFile Version 3.0\n')
    f.write('vtk file\n')
    f.write('ASCII\n')
    f.write('DATASET UNSTRUCTURED_GRID\n')

    # write points, three points for each Tubule, two points for each Protein
    f.write("POINTS "+str(len(frame.SList))+" float\n")
    for S in frame.SList:
        f.write(str(S.pos[0])+' '+str(S.pos[1]) +' '+ str(S.pos[2])+'\n')

    # write color look up data
    f.write('POINT_DATA ' +  str(len(frame.SList)) + ' \nSCALARS radius float 1 \n')
    f.write('LOOKUP_TABLE DEFAULT\n')
    for S in frame.SList:
        f.write(str(S.radius)+'\n')    # radius

# for infile in sorted(glob.glob('/mnt/xfs1/home/wyan/run_cluster/BD_Mono/L2B30T1000Pm1VP0_M/result/*.xyz')):
for infile in sorted(glob.glob(r'./*.xyz')):
    filename=os.path.split(infile)[-1]
    filenamebase=filename.split('.')[0]
    writeVTK(infile, filenamebase)
    print(filenamebase)
    print("Current File Being Processed is: " + infile)

