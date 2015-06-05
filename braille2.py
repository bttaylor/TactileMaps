from makeDome import makeCircle
from shapely.geometry import Point, box
import elevation2stl
import MapMaker3d

printSize = 80;
floor = 0;
base_height = 1;
dot_r = 0.75;


def workflow2(fname):
	streets = [	{'label':'Butler', 	'x':37, 	'y':65,	'angle':-60},
				{'label':'42nd',	'x':76,		'y':23,	'angle':45},
				{'label':'44th',	'x':20,		'y':17,	'angle':45},
				{'label':'45th',	'x':53,		'y':75,	'angle':45},
				{'label':'46th',	'x':14,		'y':63,	'angle':45}]
	streetLabels(streets,fname)

def workflow(fname):
	streets = [	{'label':'Butler', 'x': 49, 'y': 22, 'angle':60},
				{'label':'42nd','x':4,'y':34,'angle':-45},
				{'label':'44th','x':60,'y':27,'angle':-45},
				{'label':'45th','x':30,'y':86,'angle':-45},
				{'label':'46th','x':65,'y':75,'angle':-45}]
	streetLabels(streets,fname)

def streetLabels(streets,fname):
	elevation2stl.startSTL(fname)
	label_cnt = 0;
	frames = []
	for street in streets:
		name = street['label'];
		pos = (street['x'],street['y'])
		angle = street['angle']
		print 'adding: ' + name
		label_cnt = label_cnt + 1
		frame = buildBasePolyFile(name,pos,fname,3,angle,label_cnt)
		frames.append(frame)
		#stl(name,fname,3,pos,angle)
		string2braille(name,fname,pos,3,angle)
	buildFramePolyFile(fname,frames)
	stl(name,fname,3,pos,angle,label_cnt)
	elevation2stl.endSTL


def total(string,fname,coarseness):
	elevation2stl.startSTL(fname)
	#buildBasePolyFile(string,(1,1),fname,coarseness);
	angle = 45;
	stl(string,fname,coarseness,angle)
	string2braille(string,fname,(20,20),coarseness,angle)
	elevation2stl.endSTL(fname)


def stl(string,fname,coarseness,pos,angle,label_cnt):
	global printSize
	import triangle
	import triangle.plot
	import matplotlib.pyplot as plt
	import numpy as np
	from shapely.geometry import box, LineString
	#buildBasePolyFile(string,pos,fname,coarseness,angle)


	print '/Users/bttaylor/Documents/TactileMap/PolyFiles/b_' + fname + '_frame'
	frames = triangle.get_data('/Users/bttaylor/Documents/TactileMap/PolyFiles/b_' + fname + '_frame');
	framesT = triangle.triangulate(frames,'p');
	for tri in framesT['triangles']:
		pt1 = Point(framesT['vertices'][tri[0]][0],framesT['vertices'][tri[0]][1])
		pt2 = Point(framesT['vertices'][tri[1]][0],framesT['vertices'][tri[1]][1])
		pt3 = Point(framesT['vertices'][tri[2]][0],framesT['vertices'][tri[2]][1])
		tri_v1 = [pt1.x,pt1.y,base_height]
		tri_v2 = [pt2.x,pt2.y,base_height]
		tri_v3 = [pt3.x,pt3.y,base_height]
		elevation2stl.printTriangle(tri_v1,tri_v2,tri_v3,fname)

	for i in range(1,label_cnt + 1):
		print '/Users/bttaylor/Documents/TactileMap/PolyFiles/b_' + fname + '_' + str(i)
		nondot = triangle.get_data('/Users/bttaylor/Documents/TactileMap/PolyFiles/b_' + fname + '_' + str(i));
		non = triangle.triangulate(nondot,'p');
		for tri in non['triangles']:
			pt1 = Point(non['vertices'][tri[0]][0],non['vertices'][tri[0]][1])
			pt2 = Point(non['vertices'][tri[1]][0],non['vertices'][tri[1]][1])
			pt3 = Point(non['vertices'][tri[2]][0],non['vertices'][tri[2]][1])
			tri_v1 = [pt1.x,pt1.y,base_height]
			tri_v2 = [pt2.x,pt2.y,base_height]
			tri_v3 = [pt3.x,pt3.y,base_height]
			elevation2stl.printTriangle(tri_v1,tri_v2,tri_v3,fname)



	# nondot = triangle.get_data('/Users/bttaylor/Documents/TactileMap/PolyFiles/b_' + fname);
	# #print nondot
	# edges = box(0,0,printSize,printSize);
	# edge_pts = []

	# non = triangle.triangulate(nondot,'p');
	# for tri in non['triangles']:
	# 	pt1 = Point(non['vertices'][tri[0]][0],non['vertices'][tri[0]][1])
	# 	pt2 = Point(non['vertices'][tri[1]][0],non['vertices'][tri[1]][1])
	# 	pt3 = Point(non['vertices'][tri[2]][0],non['vertices'][tri[2]][1])
	# 	tri_v1 = [pt1.x,pt1.y,base_height]
	# 	tri_v2 = [pt2.x,pt2.y,base_height]
	# 	tri_v3 = [pt3.x,pt3.y,base_height]
	# 	elevation2stl.printTriangle(tri_v1,tri_v2,tri_v3,fname)
	# for seg in non['segments']:
	# 	line = LineString(non['vertices'][seg])
	# 	if line.touches(edges):
	# 		XYpt1 = Point(non['vertices'][seg[0]])
	# 		XYpt2 = Point(non['vertices'][seg[1]])
	# 		edge_pts.append(XYpt1.coords[0]);
	# 		edge_pts.append(XYpt2.coords[0]);
	# #print edge_pts
	# edge_pts.append([0,0])
	# edge_pts.append([0,printSize])
	# edge_pts.append([printSize,0])
	# edge_pts.append([printSize,printSize])
	# edge_pt_array = np.array(edge_pts)
	# base = dict(vertices=edge_pt_array);
	# baseT = triangle.triangulate(base);
	# #triangle.plot.compare(plt,base,baseT);


	# for tri in baseT['triangles']:
	# 	pt1 = Point(baseT['vertices'][tri[0]][0],baseT['vertices'][tri[0]][1])
	# 	pt2 = Point(baseT['vertices'][tri[1]][0],baseT['vertices'][tri[1]][1])
	# 	pt3 = Point(baseT['vertices'][tri[2]][0],baseT['vertices'][tri[2]][1])
	# 	tri_v1 = [pt1.x,pt1.y,floor]
	# 	tri_v2 = [pt2.x,pt2.y,floor]
	# 	tri_v3 = [pt3.x,pt3.y,floor]

	# 	elevation2stl.printTriangle(tri_v1,tri_v2,tri_v3,fname)

	#print edges
	pt1 = [0, 0, floor]
	pt2 = [0, 0, base_height]
	pt3 = [0, printSize, floor]
	pt4 = [0, printSize, base_height]
	elevation2stl.printTriangle(pt1,pt2,pt4,fname)
	elevation2stl.printTriangle(pt1,pt3,pt4,fname)
	pt3 = [printSize, 0, floor]
	pt4 = [printSize, 0, base_height]
	elevation2stl.printTriangle(pt1,pt2,pt4,fname)
	elevation2stl.printTriangle(pt1,pt3,pt4,fname)
	pt1 = [printSize, printSize, floor]
	pt2 = [printSize, printSize, base_height]
	elevation2stl.printTriangle(pt1,pt2,pt4,fname)
	elevation2stl.printTriangle(pt1,pt3,pt4,fname)
	pt3 = [0, printSize, floor]
	pt4 = [0, printSize, base_height]
	elevation2stl.printTriangle(pt1,pt2,pt4,fname)
	elevation2stl.printTriangle(pt1,pt3,pt4,fname)
	#print bottom
	pt1 = [0, 0, floor]
	pt2 = [0, printSize, floor]
	pt3 = [printSize, printSize, floor]
	pt4 = [printSize, 0, floor]
	elevation2stl.printTriangle(pt1,pt3,pt2,fname)
	elevation2stl.printTriangle(pt1,pt4,pt3,fname)


	#plt.show()

def buildFramePolyFile(fname,frames):
	f = open(MapMaker3d.directory +  'PolyFiles/b_' + fname + '_frame.poly','w')
	f.write("# Braille " + fname + "_frame.poly\n")
	vertices = 4 + 4*len(frames)
	f.write(str(vertices) + ' 2 0 1\n')

	f.write("#the exteriors of the frame\n")
	f.write('1 0 0 2\n')
	f.write('2 0 ' + str(printSize) + ' 2\n')
	f.write('3 ' + str(printSize) + ' ' + str(printSize) + ' 2\n')
	f.write('4 ' + str(printSize) + ' 0 2\n')
	f.write("#the end of the exteriors\n")

	f_cnt = 1;
	v_cnt = 5;
	boundary = 3;
	f.write("#the interior sections\n")
	for frame in frames:
		f.write("#Hole " + str(f_cnt) + '\n')
		f.write(str(v_cnt) + ' ' + str(frame.bounds[0]) + ' ' +str(frame.bounds[1]) + ' ' + str(boundary) + '\n')
		v_cnt = v_cnt + 1
		f.write(str(v_cnt) + ' ' + str(frame.bounds[0]) + ' ' +str(frame.bounds[3]) + ' ' + str(boundary) + '\n')
		v_cnt = v_cnt + 1
		f.write(str(v_cnt) + ' ' + str(frame.bounds[2]) + ' ' +str(frame.bounds[3]) + ' ' + str(boundary) + '\n')
		v_cnt = v_cnt + 1
		f.write(str(v_cnt) + ' ' + str(frame.bounds[2]) + ' ' +str(frame.bounds[1]) + ' ' + str(boundary) + '\n')
		v_cnt = v_cnt + 1
		boundary = boundary + 1
		f_cnt = f_cnt + 1

	f.write("#the segments\n")
	f.write(str(vertices) + ' 1\n')
	f.write("#exterior segments\n")
	f.write("1 1 2 2\n")
	f.write("2 2 3 2\n")
	f.write("3 3 4 2\n")
	f.write("4 4 1 2\n")

	f_cnt = 1;
	v_cnt = 5;
	boundary = 3;
	f.write("#the interior segments\n")
	for frame in frames:
		f.write("#Interior Segment " + str(f_cnt) + '\n')
		f.write(str(v_cnt) + ' ' + str(v_cnt) + ' ' +str(v_cnt + 1) + ' ' + str(boundary) + '\n')
		v_cnt = v_cnt + 1
		f.write(str(v_cnt) + ' ' + str(v_cnt) + ' ' +str(v_cnt + 1) + ' ' + str(boundary) + '\n')
		v_cnt = v_cnt + 1
		f.write(str(v_cnt) + ' ' + str(v_cnt) + ' ' +str(v_cnt + 1) + ' ' + str(boundary) + '\n')
		v_cnt = v_cnt + 1
		f.write(str(v_cnt) + ' ' + str(v_cnt) + ' ' +str(v_cnt - 3) + ' ' + str(boundary) + '\n')
		v_cnt = v_cnt + 1
		boundary = boundary + 1
		f_cnt = f_cnt + 1
	f.write("#Holes\n")
	f.write(str(len(frames)) + '\n')
	f_cnt = 1
	for frame in frames:
		f.write(str(f_cnt) + ' ' + str((frame.bounds[0] + frame.bounds[2])/2) + ' ' + str((frame.bounds[1] + frame.bounds[3])/2) + '\n')
		f_cnt = f_cnt + 1


def buildBasePolyFile(string,pos,fname,coarseness,angle,cnt):
	global printSize

	symbol = getSymbolicString(string,pos,coarseness,angle);
	#symbol is a collection of points making up the outermost part of each braille dot
	vertices = len(symbol)*4*coarseness + 4
	#coarseness indicates the number of points used per quandrant (4*coarseness)
	#There are the 4 corner verticies  ( + 4 )
	print 'characters: ',len(string) 
	width = (len(string)-1)*6.5 + 3.9
	#this isn't correct if there is a number. that adds another character
	height = 14.8
	if (pos[1] - 7.8) < 0:
		frame = box(pos[0],pos[1],pos[0]+width,pos[1]+6.3)
	else:
		frame = box(pos[0],pos[1]-7.8,pos[0]+width,pos[1]+6.3)
	
	#f = open(dir + 'PolyFiles/b_' + filename + '.poly','w')
	f = open(MapMaker3d.directory +  'PolyFiles/b_' + fname + '_' + str(cnt) + '.poly','w')
	f.write("# Braille " + fname + ".poly\n")
	f.write(str(vertices) + ' 2 0 1\n')

	f.write("#the exteriors of the frame\n")
	#f.write('1 0 0 2\n')
	#f.write('2 0 ' + str(printSize) + ' 2\n')
	#f.write('3 ' + str(printSize) + ' ' + str(printSize) + ' 2\n')
	#f.write('4 ' + str(printSize) + ' 0 2\n')
	f.write('1 ' + str(pos[0]) + ' ' + str(frame.bounds[1]) + ' 2\n')
	#f.write('1 ' + str(pos[0]) + ' ' + str(pos[1] - 8.5) + ' 2\n')
	f.write('2 ' + str(pos[0]) + ' ' + str(pos[1] + 6.3) + ' 2\n')
	f.write('3 ' + str(pos[0] + width) + ' ' + str(pos[1] + 6.3) + ' 2\n')
	f.write('4 ' + str(pos[0] + width) + ' ' + str(frame.bounds[1]) + ' 2\n')
	#f.write('4 ' + str(pos[0] + width) + ' ' + str(pos[1] - 8.5) + ' 2\n')
	f.write("#the end of the exteriors\n")

	v_cnt = 5;
	boundary = 3;
	dot_cnt = 1;

	hole_centers = []
	for dot in symbol:
		f.write("#Hole " + str(dot_cnt) + "\n")
		dot_cnt = dot_cnt + 1
		hole_x = 0;
		hole_y = 0;
		for point in dot:
			hole_x = hole_x + point.x;
			hole_y = hole_y + point.y;
			f.write(str(v_cnt) + ' ' + str(point.x) + ' ' + str(point.y) + ' ' + str(boundary) + '\n')
			v_cnt = v_cnt + 1
		hole_x = hole_x / len(dot)
		hole_y = hole_y / len(dot)
		hole_centers.append((hole_x,hole_y));
		boundary = boundary + 1
	f.write("#the end of holes\n")
	f.write("#the segments\n")
	f.write(str(vertices) + ' 1\n')
	f.write("#exterior segments\n")

	f.write('1 1 2 2\n')
	f.write('2 2 3 2\n')
	f.write('3 3 4 2\n')
	f.write('4 4 1 2\n')

	boundary = 3;
	f.write("#interior segments\n")
	cnt = 5;
	for dot in symbol:
		f.write('#interior segment ' + str(boundary - 2) + '\n')
		for i in range(0,len(dot)):
			if(i < len(dot) - 1):
				f.write(str(i + cnt) + ' ' + str(i + cnt) + ' ' + str(i+ cnt + 1) + ' ' + str(boundary) + '\n')
			else:
				f.write(str(i + cnt) + ' ' + str(i + cnt) + ' ' + str(cnt) + ' ' + str(boundary) + '\n')
		cnt = cnt + len(dot) 
		boundary = boundary + 1

	f.write("# X holes\n")
	f.write(str(len(hole_centers)) + '\n')
	cnt = 1
	for hole in hole_centers:
		f.write(str(cnt) + ' ' + str(hole[0]) + ' ' + str(hole[1]) + '\n')
		cnt = cnt + 1
	return frame

def string2braille(string,fname,pos,coarseness,angle):
	#elevation2stl.startSTL(fname)
	#pos = (0,0)
	prev_num = False;
	for char in string:
		if ord(char) >= 48 and ord(char) <= 57:
			#is a number
			if not prev_num:
				pos = char2braille('#',pos,fname,coarseness,angle)

			if ord(char) == 48:
				#'0' -> 'j'
				char = 'j'
			else:
				char = chr(ord(char) + 48) 
			pos = char2braille(char,pos,fname,coarseness,angle)
			prev_num = True;
		else:
			pos = char2braille(char,pos,fname,coarseness,angle)
			prev_num = False
	#elevation2stl.endSTL(fname)

def getSymbolicString(string,pos,coarseness,angle):
	import math
	import numpy as np
	rad = math.radians(angle);
	rotation_matrix = np.array([[math.cos(rad), math.sin(-rad)],[math.sin(rad), math.cos(rad)]]);
	symbols = [];
	prev_num = False;
	for char in string:
		#symbols = symbols + getSymbolicChar(char,pos,coarseness,angle)
		if ord(char) >= 48 and ord(char) <= 57:
			if not prev_num:
				symbols = symbols + getSymbolicChar('#',pos,coarseness,angle)
				offset = np.array([6.5,0])
				offset = np.dot(rotation_matrix,offset);
				pos = (pos[0] + offset[0], pos[1] + offset[1])

			if ord(char) == 48:
				#'0' -> 'j'
				char = 'j'
			else:
				char = chr(ord(char) + 48) 

			symbols = symbols + getSymbolicChar(char,pos,coarseness,angle)
			prev_num = True
			offset = np.array([6.5,0])
		else:
			symbols = symbols + getSymbolicChar(char,pos,coarseness,angle)
			prev_num = False
			offset = np.array([6.5, 0]);
		offset = np.dot(rotation_matrix,offset);
		pos = (pos[0] + offset[0], pos[1] + offset[1])
	return symbols

def getSymbolicChar(char,pos,coarseness,angle):
	global dot_r
	import blar
	import math
	import numpy as np
	rad = math.radians(angle);
	rotation_matrix = np.array([[math.cos(rad), math.sin(-rad)],[math.sin(rad), math.cos(rad)]]);
	symbol = [];
	pattern = getBraillePattern(char);
	if(pattern[0]):
		center = np.array([dot_r, dot_r + 4.8]);
		center = np.dot(rotation_matrix,center);
		(tri,outer) = blar.makeCircle(dot_r,Point(pos[0] + center[0],pos[1] + center[1]),coarseness);
		#(tri,outer) = blar.makeCircle(dot_r,Point(pos[0] + dot_r * math.cos(rad),pos[1] + 4.8 * math.sin(rad)),coarseness);
		symbol.append(outer);
	if(pattern[1]):
		center = np.array([dot_r, dot_r + 2.4]);
		center = np.dot(rotation_matrix,center);
		(tri,outer) = blar.makeCircle(dot_r,Point(pos[0] + center[0],pos[1] + center[1]),coarseness);
		#(tri,outer) = blar.makeCircle(dot_r,Point(pos[0] + 0,pos[1] + 2.4),coarseness);
		symbol.append(outer);
	if(pattern[2]):
		center = np.array([dot_r, dot_r]);
		center = np.dot(rotation_matrix,center);
		(tri,outer) = blar.makeCircle(dot_r,Point(pos[0] + center[0],pos[1] + center[1]),coarseness);
		#(tri,outer) = blar.makeCircle(dot_r,Point(pos[0] + 0,pos[1] + 0),coarseness);
		symbol.append(outer);
	if(pattern[3]):
		center = np.array([dot_r + 2.4, dot_r + 4.8]);
		center = np.dot(rotation_matrix,center);
		(tri,outer) = blar.makeCircle(dot_r,Point(pos[0] + center[0],pos[1] + center[1]),coarseness);
		#(tri,outer) = blar.makeCircle(dot_r,Point(pos[0] + 2.4,pos[1] + 4.8),coarseness);
		symbol.append(outer);
	if(pattern[4]):
		center = np.array([dot_r + 2.4, dot_r + 2.4]);
		center = np.dot(rotation_matrix,center);
		(tri,outer) = blar.makeCircle(dot_r,Point(pos[0] + center[0],pos[1] + center[1]),coarseness);
		#(tri,outer) = blar.makeCircle(dot_r,Point(pos[0] + 2.4,pos[1] + 2.4),coarseness);
		symbol.append(outer);
	if(pattern[5]):
		center = np.array([dot_r + 2.4, dot_r]);
		center = np.dot(rotation_matrix,center);
		(tri,outer) = blar.makeCircle(dot_r,Point(pos[0] + center[0],pos[1] + center[1]),coarseness);
		#(tri,outer) = blar.makeCircle(dot_r,Point(pos[0] + 2.4,pos[1] + 0),coarseness);
		symbol.append(outer);

	if ord(char) >= 48 and ord(char) <= 57:
		#numbers are represented with an additional symbol

		sym2 = []
		for dot in symbol:
			dot2 = [];
			for point in dot:
				offset = np.array([6.5, 0]);
				offset = np.dot(rotation_matrix,offset);
				dot2.append(Point(point.x + offset[0],point.y + offset[1],point.z))
			sym2.append(dot2)
		symbol = getSymbolicChar('#',pos,coarseness,angle) + sym2

	return symbol

def char2braille(char,pos,fname,coarseness,angle):
	
	import math
	import numpy as np
	import blar
	rad = math.radians(angle);
	rotation_matrix = np.array([[math.cos(rad), math.sin(-rad)],[math.sin(rad), math.cos(rad)]]);

	if ord(char) >= 48 and ord(char) <= 57:
		#numbers get additional symbol
		pos = char2braille('#',pos,fname,coarseness,angle)
	pattern = getBraillePattern(char);
	if(pattern[0]):
		center = np.array([dot_r, dot_r + 4.8]);
		center = np.dot(rotation_matrix,center);
		(tri,outer) = blar.makeCircle(dot_r,Point(pos[0] + center[0],pos[1] + center[1]),coarseness);
		# makeDome.makeDome(.3,Point(pos[0] + 0,pos[1] + 4.8),coarseness);
		#(tri,outer) = makeCircle(.75,Point(pos[0] + 0,pos[1] + 4.8),coarseness);
		for i in range(0,len(tri)):
			    p1 = [tri[i][0].x,tri[i][0].y,tri[i][0].z + base_height]
			    p2 = [tri[i][1].x,tri[i][1].y,tri[i][1].z + base_height]
			    p3 = [tri[i][2].x,tri[i][2].y,tri[i][2].z + base_height]
			    elevation2stl.printTriangle(p1,p2,p3,fname)
	if(pattern[1]):
		center = np.array([dot_r, dot_r + 2.4]);
		center = np.dot(rotation_matrix,center);
		(tri,outer) = blar.makeCircle(dot_r,Point(pos[0] + center[0],pos[1] + center[1]),coarseness);
		# makeDome.makeDome(.3,Point(pos[0] + 0,pol[1] + 2.4),coarseness);
		#(tri,outer) = makeCircle(.75,Point(pos[0] + 0,pos[1] + 2.4),coarseness);
		for i in range(0,len(tri)):
			    p1 = [tri[i][0].x,tri[i][0].y,tri[i][0].z + base_height]
			    p2 = [tri[i][1].x,tri[i][1].y,tri[i][1].z + base_height]
			    p3 = [tri[i][2].x,tri[i][2].y,tri[i][2].z + base_height]
			    elevation2stl.printTriangle(p1,p2,p3,fname)
	if(pattern[2]):
		center = np.array([dot_r, dot_r ]);
		center = np.dot(rotation_matrix,center);
		(tri,outer) = blar.makeCircle(dot_r,Point(pos[0] + center[0],pos[1] + center[1]),coarseness);
		# makeDome.makeDome(.3,Point(pos[0] + 0,pos[1] + 0),coarseness);
		#(tri,outer) = makeCircle(.75,Point(pos[0] + 0,pos[1] + 0),coarseness);
		for i in range(0,len(tri)):
			    p1 = [tri[i][0].x,tri[i][0].y,tri[i][0].z + base_height]
			    p2 = [tri[i][1].x,tri[i][1].y,tri[i][1].z + base_height]
			    p3 = [tri[i][2].x,tri[i][2].y,tri[i][2].z + base_height]
			    elevation2stl.printTriangle(p1,p2,p3,fname)
	if(pattern[3]):
		center = np.array([dot_r + 2.4, dot_r + 4.8]);
		center = np.dot(rotation_matrix,center);
		(tri,outer) = blar.makeCircle(dot_r,Point(pos[0] + center[0],pos[1] + center[1]),coarseness);
		# makeDome.makeDome(.3,Point(pos[0] + 2.4,pos[1] + 4.8),coarseness);
		#(tri,outer) = makeCircle(.75,Point(pos[0] + 2.4,pos[1] + 4.8),coarseness);
		for i in range(0,len(tri)):
			    p1 = [tri[i][0].x,tri[i][0].y,tri[i][0].z + base_height]
			    p2 = [tri[i][1].x,tri[i][1].y,tri[i][1].z + base_height]
			    p3 = [tri[i][2].x,tri[i][2].y,tri[i][2].z + base_height]
			    elevation2stl.printTriangle(p1,p2,p3,fname)
	if(pattern[4]):
		center = np.array([dot_r + 2.4, dot_r + 2.4]);
		center = np.dot(rotation_matrix,center);
		(tri,outer) = blar.makeCircle(dot_r,Point(pos[0] + center[0],pos[1] + center[1]),coarseness);
		# makeDome.makeDome(.3,Point(pos[0] + 2.4,pos[1] + 2.4),coarseness);
		#(tri,outer) = makeCircle(.75,Point(pos[0] + 2.4,pos[1] + 2.4),coarseness);
		for i in range(0,len(tri)):
			    p1 = [tri[i][0].x,tri[i][0].y,tri[i][0].z + base_height]
			    p2 = [tri[i][1].x,tri[i][1].y,tri[i][1].z + base_height]
			    p3 = [tri[i][2].x,tri[i][2].y,tri[i][2].z + base_height]
			    elevation2stl.printTriangle(p1,p2,p3,fname)
	if(pattern[5]):
		center = np.array([dot_r + 2.4, dot_r]);
		center = np.dot(rotation_matrix,center);
		(tri,outer) = blar.makeCircle(dot_r,Point(pos[0] + center[0],pos[1] + center[1]),coarseness);
		#makeDome.makeDome(.3,Point(pos[0] + 2.4,pos[1] + 0),coarseness)
		#(tri,outer) = makeCircle(.75,Point(pos[0] + 2.4,pos[1] + 0),coarseness);
		for i in range(0,len(tri)):
			    p1 = [tri[i][0].x,tri[i][0].y,tri[i][0].z + base_height]
			    p2 = [tri[i][1].x,tri[i][1].y,tri[i][1].z + base_height]
			    p3 = [tri[i][2].x,tri[i][2].y,tri[i][2].z + base_height]
			    elevation2stl.printTriangle(p1,p2,p3,fname)
	offset = np.array([6.5, 0]);
	offset = np.dot(rotation_matrix,offset);
	return (pos[0] + offset[0], pos[1] + offset[1])

def getBraillePattern(char):
	#currently ignoring capitalization
	if ord(char) >= 65 and ord(char) <=90:
		char = chr(ord(char) + 32) 
	return{
		'a':[True ,False,False,False,False,False],
		'b':[True ,True ,False,False,False,False],
		'c':[True ,False,False,True ,False,False],
		'd':[True ,False,False,True ,True ,False],
		'e':[True ,False,False,False,True ,False],
		'f':[True ,True ,False,True ,False,False],
		'g':[True ,True ,False,True ,True ,False],
		'h':[True ,True ,False,False,True ,False],
		'i':[False,True ,False,True ,False,False],
		'j':[False,True ,False,True ,True ,False],
		'k':[True ,False,True ,False,False,False],
		'l':[True ,True ,True ,False,False,False],
		'm':[True ,False,True ,True ,False,False],
		'n':[True ,False,True ,True ,True ,False],
		'o':[True ,False,True ,False,True ,False],
		'p':[True ,True ,True ,True ,False,False],
		'q':[True ,True ,True ,True ,True ,False],
		'r':[True ,True ,True ,False,True ,False],
		's':[False,True ,True ,True ,False,False],
		't':[False,True ,True ,True ,True ,False],
		'u':[True ,False,True ,False,False,True ],
		'v':[True ,True ,True ,False,False,True ],
		'w':[False,True ,False,True ,True ,True ],
		'x':[True ,False,True ,True ,False,True ],
		'y':[True ,False,True ,True ,True ,True ],
		'z':[True ,False,True ,False,True ,True ],
		'1':[True ,False,False,False,False,False],
		'2':[True ,True ,False,False,False,False],
		'3':[True ,False,False,True ,False,False],
		'4':[True ,False,False,True ,True ,False],
		'5':[True ,False,False,False,True ,False],
		'6':[True ,True ,False,True ,False,False],
		'7':[True ,True ,False,True ,True ,False],
		'8':[True ,True ,False,False,True ,False],
		'9':[False,True ,False,True ,False,False],
		'0':[False,True ,False,True ,True ,False],
		'.':[False,True ,False,False,True ,True ],
		',':[False,True ,False,False,False,False],
		'?':[False,True ,True ,False,False,True ],
		';':[False,True ,True ,False,False,False],
		'!':[False,True ,True ,False,True ,False],
		'"':[False,True ,True ,False,False,True ],
		'#':[False,False,True ,True ,True ,True ],
	}.get(char,[False,False,False,False,False,False])