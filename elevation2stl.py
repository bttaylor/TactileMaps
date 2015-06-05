import csv
import numpy as np
import MapMaker3d

output_x = 1.0
output_y = 1.0
output_z = .2
#zScale = output_z/(385 - 215)
#zFloor = 215
min_thick = 0.05

filename = ''

def elevation2stl(fname):
	global filename
	filename = fname;
	startSTL(filename)
	genSTL(filename)
	endSTL(filename)


def startSTL(filename):
	f = open(MapMaker3d.directory + 'StlFiles/' + filename + '.stl','w')
	print f
	f.write("solid \"test\"\n")
	f.flush()

def endSTL(filename):
	f = open(MapMaker3d.directory + 'StlFiles/' + filename + '.stl','a')
	f.write("endsolid \"test\"\n")


def genSTL(filename):
	global zScale
  	global zFloor

  	dataList = []
  	with open(MapMaker3d.directory + 'elevationFiles/' + filename + '.csv','rb') as csvfile:
  		csvreader = csv.reader(csvfile,delimiter=',')
  		for row in csvreader:
  			#print len(row)
  			dataList.append(row[5:]);

  	X = len(dataList)		#number of lines in .csv
  	Y = len(dataList[1])	#number of entries per line
  	
  	#force the file into a square shape
  	Xstep = output_x/(X-1)
  	Ystep = output_y/(Y-1)
  	height = np.zeros((X,Y))
  	cent = np.zeros((X-1,Y-1))
  	
  	#convert string lists to numpy array
  	for i in range(0,X):
  		for j in range(0,Y):
  			height[i,j] = float(dataList[i][j])

  	zScale = output_z/(height.max() - height.min())
  	zFloor = height.min()
  	print "height Min: " 
  	print height.min()
  	print "height Max: " 
  	print height.max()
	for i in range(1,X):
		for j in range(1,Y):
			cent[i-1,j-1] = (height[i-1,j-1] + height[i-1,j] + height[i,j-1] + height[i,j])/4
			tL = [(i-1)*Xstep, (j-1)*Ystep, (height[i-1,j-1] - zFloor)*zScale]
			tR = [(i)*Xstep,   (j-1)*Ystep, (height[i,j-1] - zFloor)*zScale]
			bL = [(i-1)*Xstep, (j)*Ystep,   (height[i-1,j] - zFloor)*zScale]
			bR = [i*Xstep,     (j)*Ystep,   (height[i,j] - zFloor)*zScale]
			c =  [(i - .5)*Xstep, (j - .5)*Ystep, (cent[i-1,j-1] - zFloor)*zScale]
			printTriangle(tL,c,tR)
			printTriangle(bL,c,tL)
			printTriangle(bR,c,bL)
			printTriangle(tR,c,bR)

	bXstep = output_x/(X-2)
	bYstep = output_y/(Y-2)

	#build a wall
	for i in range(1,X-1):
		tL = [(i-1)*Xstep, 0, (height[i-1,0] - zFloor)*zScale]
		tR = [i*Xstep,     0, (height[i,0] - zFloor)*zScale]
		bL = [(i-1)*bXstep,0, -1*min_thick]
		bR = [i*bXstep,    0, -1*min_thick]
		printTriangle(bL,tR,tL)
		printTriangle(bL,bR,tR)
	printTriangle([1,0,-1*min_thick],[1,0,(height[X-1,0] - zFloor)*zScale], [(X-2)*Xstep, 0, (height[X-2,0] - zFloor)*zScale])

	#build opposite wall
	for i in range(1,X-1):
		tL = [(i-1)*Xstep, output_y, (height[i-1,Y-1] - zFloor)*zScale]
		tR = [i*Xstep,     output_y, (height[i,Y-1] - zFloor)*zScale]
		bL = [(i-1)*bXstep,output_y, -1*min_thick]
		bR = [i*bXstep,    output_y, -1*min_thick]
		printTriangle(bL,tR,tL)
		printTriangle(bL,bR,tR)
	printTriangle([1,output_y,-1*min_thick],[1,output_y,(height[X-1,Y-1] - zFloor)*zScale], [(X-2)*Xstep, output_y, (height[X-2,Y-1] - zFloor)*zScale])

  	#build a wall
	for j in range(1,Y-1):
		tL = [0, (j-1)*Ystep, (height[0,j-1] - zFloor)*zScale]
		tR = [0, j*Ystep,     (height[0,j] - zFloor)*zScale]
		bL = [0, (j-1)*bYstep, -1*min_thick]
		bR = [0, j*bYstep, -1*min_thick]
		printTriangle(bL,tR,tL)
		printTriangle(bL,bR,tR)
	printTriangle([0,1,-1*min_thick],[0,1,(height[0,Y-1] - zFloor)*zScale], [0, (Y-2)*Ystep, (height[0,Y-2] - zFloor)*zScale])

  	#build a wall
	for j in range(1,Y-1):
		tL = [output_x, (j-1)*Ystep, (height[X-1,j-1] - zFloor)*zScale]
		tR = [output_x, j*Ystep,     (height[X-1,j] - zFloor)*zScale]
		bL = [output_x, (j-1)*bYstep, -1*min_thick]
		bR = [output_x, j*bYstep, -1*min_thick]
		printTriangle(bL,tR,tL)
		printTriangle(bL,bR,tR)
	printTriangle([output_x,1,-1*min_thick],[output_x,1,(height[X-1,Y-1] - zFloor)*zScale], [output_x, (Y-2)*Ystep, (height[X-1,Y-2] - zFloor)*zScale])

	printTriangle([0,0,-1*min_thick],[output_x,output_y,-1*min_thick],[output_x,0,-1*min_thick])
	printTriangle([0,0,-1*min_thick],[0,output_y,-1*min_thick],[output_x,output_y,-1*min_thick])
  

def printTriangleMirror(p1,p2,p3,name,width):  
	if name != 'default':
		filename = name;

	#need to reverse along longitude
	p1 = [width-p1[0], p1[1], p1[2]]
	p2 = [width-p2[0], p2[1], p2[2]]
	p3 = [width-p3[0], p3[1], p3[2]]
	#printTriangleBin(p1,p2,p3)

	f = open(MapMaker3d.directory + 'StlFiles/' + filename + '.stl','a')
	#f = open(filename + '.stl','a')
	#f = open('AAA' + '.stl','a')
	f.write("\tfacet normal 0 1 0\n")
	f.write("\t\touter loop\n")
	f.write("\t\t\tvertex " + str(p1[0]) + " " + str(p1[1]) + " " + str(p1[2]) + "\n")
	f.write("\t\t\tvertex " + str(p2[0]) + " " + str(p2[1]) + " " + str(p2[2]) + "\n")
	f.write("\t\t\tvertex " + str(p3[0]) + " " + str(p3[1]) + " " + str(p3[2]) + "\n")
	f.write("\t\tendloop\n")
	f.write("\tendfacet\n")
	f.flush()
	f.close()


def printTriangle(p1,p2,p3,name='default'):  
	if name != 'default':
		filename = name;

	#need to reverse along longitude
	p1 = [p1[0], p1[1], p1[2]]
	p2 = [p2[0], p2[1], p2[2]]
	p3 = [p3[0], p3[1], p3[2]]
	#printTriangleBin(p1,p2,p3)

	f = open(MapMaker3d.directory + 'StlFiles/' + filename + '.stl','a')
	#f = open(filename + '.stl','a')
	#f = open('AAA' + '.stl','a')
	f.write("\tfacet normal 0 1 0\n")
	f.write("\t\touter loop\n")
	f.write("\t\t\tvertex " + str(p1[0]) + " " + str(p1[1]) + " " + str(p1[2]) + "\n")
	f.write("\t\t\tvertex " + str(p2[0]) + " " + str(p2[1]) + " " + str(p2[2]) + "\n")
	f.write("\t\t\tvertex " + str(p3[0]) + " " + str(p3[1]) + " " + str(p3[2]) + "\n")
	f.write("\t\tendloop\n")
	f.write("\tendfacet\n")
	f.flush()
	f.close()


def printTriangleBin(p1,p2,p3,name='default'):
	if name != 'default':
		filename = name;

	from array import array
	f = open(MapMaker3d.directory + 'StlFiles/' + filename + '.stl','a')
	#f = open(filename + '_bin.stl', 'a')
	#f = open('whatthefuck.stl', 'a')
	float_array = array('f',[0.0,0.0,0.0,p1[0],p1[1],p1[2],p2[0],p2[1],p2[2],p3[0],p3[1],p3[2]]);
	#float_array = array('f',[0])
	float_array.tofile(f)
	follow = bytearray(2)
	f.write(follow)
	f.close();

def startSTLBin(filename):
	import struct

	f = open(MapMaker3d.directory + 'StlFiles/' + filename + '.stl','w')

	#Header is 80 bytes. Disregarded by modeling software
	header = '********************This is the header for a binary STL file********************'
	f.write(header)

	#4 byte unsigned long integer specifies the # of facets
	f.write(struct.pack('I',0));
	
	f.flush()

def endSTLBin(filename,count):
	import struct
	f = open(MapMaker3d.directory + 'StlFiles/' + filename + '.stl','r+')
	f.seek(80,0)
	f.write(struct.pack('I',count));