from bs4 import BeautifulSoup
from shapely.geometry import box, Polygon, LineString, Point, MultiLineString, MultiPolygon
from shapely.ops import cascaded_union
import Image, ImageDraw
from scipy.interpolate import interp1d, interp2d
import btpoly_tri
import json

#imgSize = 1024
#minlat = 0
#minlon = 0
#maxlat = 0
#maxlon = 0

directory = '/Users/bttaylor/Documents/TactileMap/'
#directory= '/usr0/home/bttaylor/TactileMap/'

def workflow(name,lon,lat,*locs,**attr):
	import downloadMap
	import getElevation
	import shutil

	filename = setParams(name,lon,lat,*locs,**attr);
	print 'in workflow'
	params = json.load(open(directory + 'ParamFiles/' + filename + '.txt'))
	if params['elevation']:
		getElevation.getElevation(params['minlon'],params['maxlat'],params['maxlon'],params['minlat'],filename)

	downloadMap.downloadMap(params['minlon'],params['maxlat'],params['maxlon'],params['minlat'],filename)
	roads = getRoads(filename);
	road_poly = roads2poly(filename,roads)
	
	poly2polyFile(filename,road_poly,False) #Think this is always false.
	grid = box(params['minlon'],params['minlat'],params['maxlon'],params['maxlat'])
	ipoly = grid.difference(road_poly)
	poly2polyFile(filename,ipoly,True)

	#polyFile2stl(filename,poly)
	polyFile2stlBin(filename,road_poly)
	shutil.copy2(directory + 'StlFiles/' + filename + '.stl','./' + params['name'] + '.stl')


def workflow2(name,lon,lat,*locs,**attr):
	import downloadMap
	import getElevation
	import shutil
	import featurePrinter

	filename = setParams(name,lon,lat,*locs,**attr);
	print 'in workflow 2'
	params = json.load(open(directory + 'ParamFiles/' + filename + '.txt'))
	if params['elevation']:
		getElevation.getElevation(params['minlon'],params['maxlat'],params['maxlon'],params['minlat'],filename)

	downloadMap.downloadMap(params['minlon'],params['maxlat'],params['maxlon'],params['minlat'],filename)
	roads = getRoads(filename);
	road_poly = roads2poly(filename,roads)
	water = getWater(filename);
	if len(water) > 1:
		water = MultiPolygon(water)
	#park = getPark(filename);
	#if len(park) > 1:
	#	park = MultiPolygon(park)
	water = water.difference(road_poly);
	#park = park.difference(road_poly);
	#park = park.difference(water);

	poly2polyFile2(filename,'w',water,False)
	wcnt = 1;
	feat_tris = [];
	edge_pts = [];
	for ww in water:
		(x, edge_points) = featurePrinter.featureAreaPoly(filename,str(wcnt),ww,'Dome')
		feat_tris.append(x)
		edge_pts.append(edge_points)
		#featurePrinter.featureAreaPoly(filename,water)
		wcnt = wcnt + 1
	#poly2polyFile2(filename,'p',park,False)

	print 'feat_tris', (len(feat_tris))

	poly2polyFile(filename,road_poly,False) #Think this is always false.
	grid = box(params['minlon'],params['minlat'],params['maxlon'],params['maxlat'])
	ipoly = grid.difference(road_poly)
	ipoly = ipoly.difference(water)
	#ipoly = ipoly.difference(park)
	poly2polyFile(filename,ipoly,True)

	#polyFile2stl(filename,poly)
	polyFile2stlBin2(filename,road_poly,feat_tris,edge_pts)
	shutil.copy2(directory + 'StlFiles/' + filename + '.stl','./' + params['name'] + '.stl')


#creates a json file with relevant parameters
def setParams(name,lon,lat,*locs,**attr):
	import time
	import math
	filename = name + str(int(math.floor(time.time())))
	print "filename: " + filename
	print "lat,lon: " + "(" + str(lon) + ',' + str(lat) + ')'
	print len(locs)
	if len(locs) == 0:
		ratio = getLonLatRatio(lat)
		print 'ratio: ' + str(ratio)
		latDelta = dist2LatDegree(.5,lat);
		minlon = lon - latDelta * ratio;
		maxlon = lon + latDelta * ratio;
		minlat = lat - latDelta;
		maxlat = lat + latDelta;
	elif len(locs) == 2:
		if lon < locs[0]:
			minlon = lon;
			maxlon = locs[0];
		else:
			minlon = locs[0];
			maxlon = lon;
		if lat < locs[1]:
			minlat = lat
			maxlat = locs[1]
		else:
			minlat = locs[1]
			maxlat = lat
	if 'elev' in attr:
		print "set elevation to: " + str(attr.get('elev'))
		elevation = attr.get('elev')
	else:
		print "set elevation to default: False"
		elevation = False
	if 'zmin' in attr:
		zmin = True;
		zminVal = attr.get('zmin')
	else:
		zmin = False;
		zminVal = 0

	print minlon,maxlat,maxlon,minlat
	d = {"name":name, "minlon":minlon, 'minlat':minlat, 'maxlon':maxlon, 'maxlat':maxlat, 
			'elevation':elevation, 'samples':20, 'ratio':1, 'road_width':40, 'road_height':2,
			'imgSize':1024, 'base_size':100, 'base_height':1, 'print_size':100, 'dir':directory,
			'z_scale':0.2, 'zmin':zmin, 'zminVal':zminVal}
	print d
	json.dump(d, open(directory+ 'ParamFiles/' + filename + '.txt','w'))
	return filename

def dist2LatDegree(miles,lat):
	return .0036;
	# (2 * .0036) = .0072 degree latitude = ~1/2 mile at ~40.47 latitude (Pittsburgh)
	# 1/2 mile seems to be a decent distance for city block scale maps
	# .009 degree lat ~= 1km size might be good too.
	# could make this a function to get the degrees equivilant to ~1/2 mile at given lat.

def getLonLatRatio(lat):
	import math
	lat = math.radians(lat);
	# equation from http://www.csgnetwork.com/degreelenllavcalc.html
	m1 = 111132.92;		# latitude calculation term 1
	m2 = -559.82;		# latitude calculation term 2
	m3 = 1.175;			# latitude calculation term 3
	m4 = -0.0023;		# latitude calculation term 4
	p1 = 111412.84;		# longitude calculation term 1
	p2 = -93.5;			# longitude calculation term 2
	p3 = 0.118;			# longitude calculation term 3

	latlen = m1 + (m2 * math.cos(2*lat)) + (m3 * math.cos(4*lat)) + (m4 * math.cos(6*lat));
	lonlen = (p1 * math.cos(lat)) + (p2*math.cos(3*lat)) + (p3 * math.cos(5*lat));

	return latlen / lonlen;

def getWater(filename):
	params = json.load(open(directory + 'ParamFiles/' + filename + '.txt'))
	#global minlat, minlon, maxlat, maxlon
	minlat = params['minlat']
	minlon = params['minlon']
	maxlat = params['maxlat']
	maxlon = params['maxlon']
	#print minlon, maxlat, maxlon, minlat
	imgSize = params['imgSize']

	boundary = box(minlon,minlat,maxlon,maxlat);
	map_water = Image.new("RGB",(imgSize,imgSize),"white")
	draw_water = ImageDraw.Draw(map_water)

	f=open(directory + 'mapFiles/' + filename + '.osm')

	soup = BeautifulSoup(f)

	#ways are the features we are interested in
	ways = soup.find_all('way')
	#nodes contain the coordinates of points
	nodes = soup.find_all('node')

	maplat_img = interp1d([minlat,maxlat],[0,imgSize])
	maplon_img = interp1d([minlon,maxlon],[0,imgSize])
	#interp1d([minlon,maxlon],[0,imgSize])

	water_pts=[]
	water_objs = []
	#All features (road, water, parks anyway) are specified as ways
	for way in ways:
		name = ''
		size = 0
		Water = False
		water_pts = [];
		OutOfBounds = False
		tags = way.find_all('tag')
		#tags have attributes of the ways
		for tag in tags:
			if str(tag.get('k')) == 'name':
				#store the name of the particular way
				name = tag.get('v')
			if str(tag.get('k')) == 'waterway':
				#water has a 'waterway' tag
				if tag.get('v') == 'riverbank':
					print '*riverbank'
					#Water = True;
				Water = True;
					#print tag.get('v')
			if str(tag.get('k')) == 'natural':
				if tag.get('v') == 'water':
					Water = True;

		if(Water == True):
			print name + ' - Water'
			#nd are the nodes (lon,lat) of this particular way
			nds = way.find_all('nd')
			for nd in nds:
				loc = getLocfromNode(nodes,int(nd.get('ref')))
				water_pts.append((loc[0],loc[1]))

			#if len(water_pts) > 2:
			if water_pts[0] == water_pts[-1]:
				print 'loopeder'
				water_edges = Polygon(water_pts)
				if water_edges.is_valid:
					water_edge = water_edges.intersection(boundary);

					#for coord in water_edge.exterior.coords:

					water_objs.append(water_edge);
					loc = water_edge.exterior.coords[0]
					for coord in water_edge.exterior.coords[1:]:
						if isOutOfBounds(filename,coord[0],coord[1]):
							print 'OUT',coord
						draw_water.line((maplon_img(loc[0]),imgSize - maplat_img(loc[1]),maplon_img(coord[0]),imgSize - maplat_img(coord[1])), fill=(0,0,255))
						loc = coord

	map_water.save(directory + 'mapImages/' + filename + '_water.png')
	return water_objs


def getPark(filename):
	params = json.load(open(directory + 'ParamFiles/' + filename + '.txt'))
	#global minlat, minlon, maxlat, maxlon
	minlat = params['minlat']
	minlon = params['minlon']
	maxlat = params['maxlat']
	maxlon = params['maxlon']
	#print minlon, maxlat, maxlon, minlat
	imgSize = params['imgSize']

	boundary = box(minlon,minlat,maxlon,maxlat);
	map_park = Image.new("RGB",(imgSize,imgSize),"white")
	draw_park = ImageDraw.Draw(map_park)

	f=open(directory + 'mapFiles/' + filename + '.osm')

	soup = BeautifulSoup(f)

	#ways are the features we are interested in
	ways = soup.find_all('way')
	#nodes contain the coordinates of points
	nodes = soup.find_all('node')

	maplat_img = interp1d([minlat,maxlat],[0,imgSize])
	maplon_img = interp1d([minlon,maxlon],[0,imgSize])
	#interp1d([minlon,maxlon],[0,imgSize])

	park_pts=[]
	park_objs = []
	#All features (road, water, parks anyway) are specified as ways
	for way in ways:
		name = ''
		size = 0
		Park = False
		park_pts = [];
		OutOfBounds = False
		tags = way.find_all('tag')
		#tags have attributes of the ways
		for tag in tags:
			if str(tag.get('k')) == 'name':
				#store the name of the particular way
				name = tag.get('v')
			if str(tag.get('k')) == 'leisure':
				print tag.get('v')
				if tag.get('v') == 'park':
					Park = True;
		
		if(Park == True):
			print name + ' - Park'
			#nd are the nodes (lon,lat) of this particular way
			nds = way.find_all('nd')
			for nd in nds:
				loc = getLocfromNode(nodes,int(nd.get('ref')))
				park_pts.append((loc[0],loc[1]))

			park_edges = Polygon(park_pts)
			if park_edges.is_valid:
				park_edges = park_edges.intersection(boundary);
				park_objs.append(park_edges);
				loc = park_edges.exterior.coords[0]
				for coord in park_edges.exterior.coords[1:]:
					if isOutOfBounds(filename,coord[0],coord[1]):
						print 'OUT',coord
					draw_park.line((maplon_img(loc[0]),imgSize - maplat_img(loc[1]),maplon_img(coord[0]),imgSize - maplat_img(coord[1])), fill=(0,255,0))
					loc = coord

	map_park.save(directory + 'mapImages/' + filename + '_park.png')
	return park_objs


#Opens an existing (i.e. already downloaded) .osm file
#Returns a shapely MultiLineString containing a LineString for each road
def getRoads(filename):
	params = json.load(open(directory + 'ParamFiles/' + filename + '.txt'))
	#global minlat, minlon, maxlat, maxlon
	minlat = params['minlat']
	minlon = params['minlon']
	maxlat = params['maxlat']
	maxlon = params['maxlon']
	print minlon, maxlat, maxlon, minlat
	imgSize = params['imgSize']

	map_road = Image.new("RGB",(imgSize,imgSize),"white")
	draw_road = ImageDraw.Draw(map_road)
	
	f=open(directory + 'mapFiles/' + filename + '.osm')

	soup = BeautifulSoup(f)

	#ways are the features we are interested in
	ways = soup.find_all('way')
	#nodes contain the coordinates of points
	nodes = soup.find_all('node')

	maplat_img = interp1d([minlat,maxlat],[0,imgSize])
	maplon_img = interp1d([minlon,maxlon],[0,imgSize])
	#interp1d([minlon,maxlon],[0,imgSize])
	print maplat_img(maxlat)

	road_lines=[]
	road_objs = []
	road_pts =[]
	#All features (road, water, parks anyway) are specified as ways
	for way in ways:
		name = ''
		size = 0
		Road = False
		#road_pts = []
		OutOfBounds = False
		tags = way.find_all('tag')
		#tags have attributes of the ways
		for tag in tags:
			if str(tag.get('k')) == 'name':
				#store the name of the particular way
				name = tag.get('v')
			if str(tag.get('k')) == 'highway':
				#roads have 'highway' tags
				size = roadSize(str(tag.get('v')))
				if( size > 0):
					#filter out small roads
					Road = True

		if(Road == True):
			#print name + ' - Road'
			#nd are the nodes (lon,lat) of this particular way
			nds = way.find_all('nd')
			pt1ind = 0
			#we assume the first point (nd) is not located within the bounds of our map
			startInbounds = False
			while startInbounds == False:
				loc = getLocfromNode(nodes,int(nds[pt1ind].get('ref')))
				lon1 = loc[0]
				lat1 = loc[1]
				if(isOutOfBounds(filename,lon1,lat1) == False):
					#first point in bounds, start plotting
					startInbounds = True
				else:
					loc = getLocfromNode(nodes,int(nds[pt1ind+1].get('ref')))
					lon2 = loc[0]
					lat2 = loc[1]
					if(isOutOfBounds(filename,lon2,lat2) == False):
						#point 2 is in bounds, plot point 1 on border
						loc = adjustPoints(filename,lon1,lat1,lon2,lat2)
						lon1 = loc[0]
						lat1 = loc[1]
						startInbounds = True
					else:
						#ignore point 1 and move on
						pt1ind = pt1ind + 1
			#Point 1 is now in bounds (or on border)
			road_pts = [Point(lon1,lat1)]
			#Plot remaining points
			for nd in nds[pt1ind+1:len(nds)]:
				loc = getLocfromNode(nodes,int(nd.get('ref')))
				lon2 = loc[0]
				lat2 = loc[1]
				if( (isOutOfBounds(filename,lon1,lat1) == True) & (isOutOfBounds(filename,lon2,lat2) == True) ):
					print "BOTH OUT OF BOUNDS. Some error?"
					OutOfBounds = True
				elif isOutOfBounds(filename,lon2,lat2) == True:
					#2nd Point is OOB, needs to be plotted to border
					loc = adjustPoints(filename,lon2,lat2,lon1,lat1)
					lon2 = loc[0]
					lat2 = loc[1]
					if((lon1 == lon2) & (lat1 == lat2)):
						#if the road goes out of bounds, we stop adding points
						#otherwise it will just interpolate it to the edge and we'll have a bunch of identical points
						#might be an issue if the road exits and re-enters the boundaries
						OutOfBounds = True
				elif isOutOfBounds(filename,lon1,lat1) == True:
					print "Point 1 OOB. error?"
					loc = adjustPoints(filename,lon1,lat1,lon2,lat2)
					lon1 = loc[0]
					lat1 = loc[1]
				if (OutOfBounds == False):
					road_pts.append(Point(lon2,lat2))
					draw_road.line((maplon_img(lon1),imgSize - maplat_img(lat1),maplon_img(lon2),imgSize - maplat_img(lat2)), fill=(0,0,0))
					#road_dxf.add(dxf.line( (maplon_img(lon1),imgSize - maplat_img(lat1)), (maplon_img(lon2),imgSize - maplat_img(lat2)), color=7))
				lon1 = lon2
				lat1 = lat2
			road_objs.append({'name':name,'size':size,'path':LineString(road_pts)})
			road_lines.append(LineString(road_pts));

	allRoads = MultiLineString(road_lines);

	map_road.save(directory + 'mapImages/' + filename + '_road.png')
	return road_objs
	#return allRoads



def roads2poly(filename,roads=None):
	params = json.load(open(directory + 'ParamFiles/' + filename + '.txt'))
	#global minlat, minlon, maxlat, maxlon
	imgSize = params['imgSize']
	map_road = Image.new("RGB",(imgSize,imgSize),"white")
	draw_road = ImageDraw.Draw(map_road)

	if roads is None:
		roads = getRoads(filename)

	minlat = params['minlat']
	minlon = params['minlon']
	maxlat = params['maxlat']
	maxlon = params['maxlon']
	width = params['road_width']

	grid = box(minlon,minlat,maxlon,maxlat)
	#width = 10
	maplat_img = interp1d([minlat,maxlat],[imgSize,0])
	maplon_img = interp1d([minlon,maxlon],[0,imgSize])
	
	print "number of roads: " + str(len(roads))

	road_outlines = []
	for road in roads:
		#print road
		#this metric is in the road points (lat, lon), needs to convert to road size?
		#buff_width = width * .000002
		#this converts the roads to a fixed latitude degree unit which is later converted to an image based
		#on maplat_img.  There might be an issue w/ lat and lon being stretched differently?
		buff_width = width * .00025 * (maxlat - minlat);
		#this should fix the road to a set size relative to the scaled image so that different zooms will 
		#produce the same road widths

		#I didn't see a difference between cap styles, but the round (=1) JointStyle better avoids narrow cracks
		fat = road['path'].buffer(buff_width*road['size']/2,16,2,1)
		#fat = road.buffer(buff_width,16,2,2)
		fat = fat.intersection(grid)
		road_outlines.append(fat);


	road_poly = cascaded_union(road_outlines)

	print road_poly.geom_type
	if road_poly.geom_type == 'MultiPolygon':
		print "number of polygons: " + str(len(road_poly.geoms))
		color = (255,0,0)
		for single in road_poly:
			bound = single.exterior.coords
			for j in range(1, len(list(bound))):
				draw_road.line((maplon_img(bound[j-1][0]),maplat_img(bound[j-1][1]),maplon_img(bound[j][0]),maplat_img(bound[j][1])), fill=color,width=2)
				#draw_road.line((bound[j-1][0],bound[j-1][1],bound[j][0],bound[j][1]), fill=(255,0,0),width=2)
			color = (0,0,255)
			drawHoles(filename,single.interiors,draw_road)
	elif road_poly.geom_type == 'Polygon':
		print road_poly.is_valid
		#print list(road_poly.interiors)
		bound = road_poly.exterior.coords
		#print list(bound)
		color = (255,0,0)
		for j in range(1, len(list(bound))):
			#print bound[j]
			draw_road.line((maplon_img(bound[j-1][0]),maplat_img(bound[j-1][1]),maplon_img(bound[j][0]),maplat_img(bound[j][1])), fill=color,width=2)
			#draw_road.line((bound[j-1][0],bound[j-1][1],bound[j][0],bound[j][1]), fill=(255,0,0),width=2)
		drawHoles(filename,road_poly.interiors,draw_road)

	map_road.save(directory + 'PngFiles/' + filename + '.png')

	return road_poly

def drawHoles(filename,holes,draw):
	params = json.load(open(directory + 'ParamFiles/' + filename + '.txt'))
	minlat = params['minlat']
	minlon = params['minlon']
	maxlat = params['maxlat']
	maxlon = params['maxlon']
	imgSize = params['imgSize']
	#global minlat, minlon, maxlat, maxlon
	maplat_img = interp1d([minlat,maxlat],[imgSize,0])
	maplon_img = interp1d([minlon,maxlon],[0,imgSize])
	#Draw the holes w/ points located inside
	print "_____________"
	print "# of holes: " + str(len(holes))
	for hole in holes:
		p = btpoly_tri.draw_poly(hole.coords[1:]);
		x = p[0][0] + .5*(p[1][0]-p[0][0]) + .25*(p[2][0]-p[0][0])
		y = p[0][1] + .5*(p[1][1]-p[0][1]) + .25*(p[2][1]-p[0][1])
		pp = Point(x,y)
		#print pp
		draw.ellipse((maplon_img(pp.x)-5,maplat_img(pp.y)-5,maplon_img(pp.x)+5,maplat_img(pp.y)+5), fill=(0,0,255))
		
		bound = hole.coords
		for j in range(1,len(list(hole.coords))):
			draw.line((maplon_img(bound[j-1][0]),maplat_img(bound[j-1][1]),maplon_img(bound[j][0]),maplat_img(bound[j][1])), fill=(0,255,0),width=2)


def poly2polyFile(filename,poly=None,inverseRoads=False):
	Hole_Area_Threshold = 10e-20
	#hole = 0; #Need to fix this to deal w/ multiple holes
	params = json.load(open(directory + 'ParamFiles/' + filename + '.txt'))
	minlat = params['minlat']
	minlon = params['minlon']
	maxlat = params['maxlat']
	maxlon = params['maxlon']
	imgSize = params['imgSize']

	if poly is None:
		poly = roads2poly(filename);

	print minlat, minlon, maxlat, maxlon
	maplat_img = interp1d([minlat,maxlat],[imgSize,0])
	maplon_img = interp1d([minlon,maxlon],[0,imgSize])

	#Count the vertices
	vertices = 0
	if poly.geom_type == 'MultiPolygon':
		for p in poly:
			#polygon exteriors have repeated coordinates @ first,last so subtract 1
			vertices = vertices + len(p.exterior.coords) - 1
			for hole in p.interiors:
				vertices = vertices + len(hole.coords) - 1
	else:
		vertices = len(poly.exterior.coords) - 1
		for hole in poly.interiors:
			vertices = vertices + len(hole.coords) - 1

	segments = vertices

	if params['elevation']:
		if inverseRoads:
			#import testipoly
			vs = getElevationVertices(filename,poly)
			print len(vs)
			vertices = vertices + len(vs)
	else:
		vs = []

	#Create .poly file
	if inverseRoads:
		f = open(directory + 'PolyFiles/i_' + filename + '.poly','w')
	else:
		f = open(directory + 'PolyFiles/' + filename + '.poly','w')
	#Write Header
	f.write("#" + filename + ".poly\n")
	f.write(str(vertices) + ' 2 0 1\n')

	#Write the vertices
	if poly.geom_type == 'MultiPolygon':
		f.write("#the exteriors of the roads polygon\n")
		v_cnt = 1
		boundary = 2;
		ext_cnt = 0;
		#First write the exterior vertices
		for p in poly:


			#test for markers
			skip = False;
			

			if not skip:
				ext_cnt = ext_cnt + 1
				f.write("#exterior " + str(boundary - 1) + "\n")
				for i in range(1,len(p.exterior.coords)):
					f.write(str(v_cnt) + ' ' + str(maplon_img(p.exterior.coords[i][0])) + ' ' + str(maplat_img(p.exterior.coords[i][1])) + ' ' + str(boundary) + '\n')
					v_cnt = v_cnt + 1
				boundary = boundary + 1
		f.write("#the end of the exterior polygons\n")

		#Next write the hole vertices
		for p in poly:
			for hole in p.interiors:
				f.write("#Hole " + str(boundary - ext_cnt - 1) + "\n")
				#f.write("#Hole " + str(boundary - len(poly.geoms) - 1) + "\n")
				for i in range(1,len(hole.coords)):
					f.write(str(v_cnt) + ' ' + str(maplon_img(hole.coords[i][0])) + ' ' + str(maplat_img(hole.coords[i][1])) + ' ' + str(boundary) + '\n')
					v_cnt = v_cnt + 1
				boundary = boundary + 1

	else:
		f.write("#the exterior of the roads polygon\n")

		v_cnt = 1
		boundary = 2
		for i in range(1,len(poly.exterior.coords)):
			f.write(str(v_cnt) + ' ' + str(maplon_img(poly.exterior.coords[i][0])) + ' ' + str(maplat_img(poly.exterior.coords[i][1])) + ' ' + str(boundary) + '\n')
			v_cnt = v_cnt + 1
		boundary = boundary + 1

		f.write("#the first hole\n");
		for hole in poly.interiors:
			f.write("#Hole 1\n")
			for i in range(1,len(hole.coords)):
				f.write(str(v_cnt) + ' ' + str(maplon_img(hole.coords[i][0])) + ' ' + str(maplat_img(hole.coords[i][1])) + ' ' + str(boundary) + '\n')
				v_cnt = v_cnt + 1
			boundary = boundary + 1

	if inverseRoads:
		f.write("#Sample Points\n")
		for i in range(0,len(vs)):
			f.write(str(v_cnt) + ' ' + str(maplon_img(vs[i][0])) + ' ' + str(maplat_img(vs[i][1])) + ' 0\n')
			v_cnt = v_cnt + 1

	#write segment header
	f.write("#the segments\n")
	f.write(str(segments) + ' 1\n')
	f.write("#exterior segment\n")

	if poly.geom_type == 'MultiPolygon':
		boundary = 2;
		seg_cnt = 1;
		first_vert = 1;
		for p in poly:


			#test for markers
			skip = False;
	

			if not skip:
				f.write("#exterior seg " + str(boundary - 1) + "\n")
				for i in range(1,len(p.exterior.coords)):
					if(i < len(p.exterior.coords) - 1):
						#segments are incremental 
						f.write(str(seg_cnt) + ' ' + str(seg_cnt) + ' ' + str(seg_cnt + 1) + ' ' + str(boundary) + '\n')
					else:
						#last segment loops back to first
						f.write(str(seg_cnt) + ' ' + str(seg_cnt) + ' ' + str(first_vert) + ' ' + str(boundary) + '\n')
					seg_cnt = seg_cnt + 1
				first_vert = seg_cnt
				boundary = boundary + 1

		for p in poly:
			for hole in p.interiors:
				f.write("#interior segment " + str(boundary - ext_cnt - 1) + "\n")
				#f.write("#interior segment " + str(boundary - len(poly.geoms) - 1) + "\n")
				for i in range(1,len(hole.coords)):
					if(i < len(hole.coords) -1):
						f.write(str(seg_cnt) + ' ' + str(seg_cnt) + ' ' + str(seg_cnt + 1) + ' ' + str(boundary) + '\n')
					else:
						f.write(str(seg_cnt) + ' ' + str(seg_cnt) + ' ' + str(first_vert) + ' ' + str(boundary) + '\n')
					seg_cnt = seg_cnt + 1
				first_vert = seg_cnt
				boundary = boundary + 1
	else:
		cnt = 0
		for i in range(1,len(poly.exterior.coords)):
			if(i < len(poly.exterior.coords) - 1):
				f.write(str(i + cnt) + ' ' + str(i+cnt) + ' ' + str(i+cnt+1) + ' 2\n')
			else:
				f.write(str(i + cnt) + ' ' + str(i+cnt) + ' ' + str(1+cnt) + ' 2\n')

		cnt = cnt + len(poly.exterior.coords) - 1;

		boundary = 3;
		for hole in poly.interiors:
			f.write("#interior segment " + str(boundary - 2) + '\n')
			for i in range(1,len(hole.coords)):
				if(i < len(hole.coords) - 1):
					f.write(str(i + cnt) + ' ' + str(i+cnt) + ' ' + str(i+cnt+1) +' ' + str(boundary) + '\n')
				else:
					f.write(str(i + cnt) + ' ' + str(i+cnt) + ' ' + str(1+cnt) + ' ' + str(boundary) + '\n')
			cnt = cnt + len(hole.coords) - 1
			boundary = boundary + 1

	f.write("# X holes\n")

	hole_cnt = 0
	if poly.geom_type == 'MultiPolygon':
		for p in poly:
			hole_cnt = hole_cnt + len(p.interiors)
	else:
		hole_cnt = len(poly.interiors)

	f.write(str(hole_cnt) + '\n')
	if hole_cnt > 0:
	
		hole_num = 1
		if poly.geom_type == 'MultiPolygon':
			for p in poly:
				for hole in p.interiors:
					test_poly = Polygon(hole)

					#test for markers
					if test_poly.intersects(Point(11.24494,43.780094)):
						print hole_num, ' has the Hotel'

					if test_poly.intersects(Point(11.248853,43.781198)):
						print hole_num, ' has the Fortress'

					if test_poly.intersects(Point(11.247373,43.776893)):
						print hole_num, ' has the Station'


					if test_poly.area > Hole_Area_Threshold:
						pt = btpoly_tri.draw_poly(hole.coords[1:]);
						x = pt[0][0] + .5*(pt[1][0]-pt[0][0]) + .25*(pt[2][0]-pt[0][0])
						y = pt[0][1] + .5*(pt[1][1]-pt[0][1]) + .25*(pt[2][1]-pt[0][1])
						f.write(str(hole_num) + ' ' + str(maplon_img(x)) + ' ' + str(maplat_img(y)) +'\n')
						if not (test_poly.contains(Point(x,y))):
							print "Point " + str(hole_num) + " is not interior?"
							print pt
							print x,y
						hole_num = hole_num + 1
					else:
						print "Skipped hole " + str(hole_num) + " is too small: " + str(test_poly.area)


		else:
			for hole in poly.interiors:
				test_poly = Polygon(hole)
				if test_poly.area > Hole_Area_Threshold:
					pt = btpoly_tri.draw_poly(hole.coords[1:]);
					x = pt[0][0] + .5*(pt[1][0]-pt[0][0]) + .25*(pt[2][0]-pt[0][0])
					y = pt[0][1] + .5*(pt[1][1]-pt[0][1]) + .25*(pt[2][1]-pt[0][1])
					f.write(str(hole_num) + ' ' + str(maplon_img(x)) + ' ' + str(maplat_img(y)) +'\n')
					if not (test_poly.contains(Point(x,y))):
						print "Point " + str(hole_num) + " is not interior?"
					hole_num = hole_num + 1
				else:
					print "Skipped hole " + str(hole_num) + " is too small: " + str(test_poly.area)


def poly2polyFileItaly(filename,poly=None,inverseRoads=False):
	Hole_Area_Threshold = 10e-20
	#hole = 0; #Need to fix this to deal w/ multiple holes
	params = json.load(open(directory + 'ParamFiles/' + filename + '.txt'))
	minlat = params['minlat']
	minlon = params['minlon']
	maxlat = params['maxlat']
	maxlon = params['maxlon']
	imgSize = params['imgSize']

	if poly is None:
		poly = roads2poly(filename);

	print minlat, minlon, maxlat, maxlon
	maplat_img = interp1d([minlat,maxlat],[imgSize,0])
	maplon_img = interp1d([minlon,maxlon],[0,imgSize])

	#Count the vertices
	vertices = 0
	if poly.geom_type == 'MultiPolygon':
		for p in poly:
			#polygon exteriors have repeated coordinates @ first,last so subtract 1
			vertices = vertices + len(p.exterior.coords) - 1
			for hole in p.interiors:
				vertices = vertices + len(hole.coords) - 1
	else:
		vertices = len(poly.exterior.coords) - 1
		for hole in poly.interiors:
			vertices = vertices + len(hole.coords) - 1

	segments = vertices

	if params['elevation']:
		if inverseRoads:
			#import testipoly
			vs = getElevationVertices(filename,poly)
			print len(vs)
			vertices = vertices + len(vs)
	else:
		vs = []

	#Create .poly file
	if inverseRoads:
		f = open(directory + 'PolyFiles/i_' + filename + '.poly','w')
	else:
		f = open(directory + 'PolyFiles/' + filename + '.poly','w')
	#Write Header
	f.write("#" + filename + ".poly\n")
	f.write(str(vertices) + ' 2 0 1\n')

	#Write the vertices
	if poly.geom_type == 'MultiPolygon':
		f.write("#the exteriors of the roads polygon\n")
		v_cnt = 1
		boundary = 2;
		ext_cnt = 0;
		#First write the exterior vertices
		for p in poly:


			#test for markers
			skip = False;
			if p.intersects(Point(11.24494,43.780094)):
				print boundary, ' has the Hotel', len(p.exterior.coords)
				skip = True;

			if p.intersects(Point(11.248853,43.781198)):
				print boundary, ' has the Fortress', len(p.exterior.coords)
				skip = True;

			if p.intersects(Point(11.247373,43.776893)):
				print boundary, ' has the Station', len(p.exterior.coords)
				skip = True;

			if not skip:
				ext_cnt = ext_cnt + 1
				f.write("#exterior " + str(boundary - 1) + "\n")
				for i in range(1,len(p.exterior.coords)):
					f.write(str(v_cnt) + ' ' + str(maplon_img(p.exterior.coords[i][0])) + ' ' + str(maplat_img(p.exterior.coords[i][1])) + ' ' + str(boundary) + '\n')
					v_cnt = v_cnt + 1
				boundary = boundary + 1
		f.write("#the end of the exterior polygons\n")

		#Next write the hole vertices
		for p in poly:
			for hole in p.interiors:
				f.write("#Hole " + str(boundary - ext_cnt - 1) + "\n")
				#f.write("#Hole " + str(boundary - len(poly.geoms) - 1) + "\n")
				for i in range(1,len(hole.coords)):
					f.write(str(v_cnt) + ' ' + str(maplon_img(hole.coords[i][0])) + ' ' + str(maplat_img(hole.coords[i][1])) + ' ' + str(boundary) + '\n')
					v_cnt = v_cnt + 1
				boundary = boundary + 1

	else:
		f.write("#the exterior of the roads polygon\n")

		v_cnt = 1
		boundary = 2
		for i in range(1,len(poly.exterior.coords)):
			f.write(str(v_cnt) + ' ' + str(maplon_img(poly.exterior.coords[i][0])) + ' ' + str(maplat_img(poly.exterior.coords[i][1])) + ' ' + str(boundary) + '\n')
			v_cnt = v_cnt + 1
		boundary = boundary + 1

		f.write("#the first hole\n");
		for hole in poly.interiors:
			f.write("#Hole 1\n")
			for i in range(1,len(hole.coords)):
				f.write(str(v_cnt) + ' ' + str(maplon_img(hole.coords[i][0])) + ' ' + str(maplat_img(hole.coords[i][1])) + ' ' + str(boundary) + '\n')
				v_cnt = v_cnt + 1
			boundary = boundary + 1

	if inverseRoads:
		f.write("#Sample Points\n")
		for i in range(0,len(vs)):
			f.write(str(v_cnt) + ' ' + str(maplon_img(vs[i][0])) + ' ' + str(maplat_img(vs[i][1])) + ' 0\n')
			v_cnt = v_cnt + 1

	#write segment header
	f.write("#the segments\n")
	f.write(str(segments) + ' 1\n')
	f.write("#exterior segment\n")

	if poly.geom_type == 'MultiPolygon':
		boundary = 2;
		seg_cnt = 1;
		first_vert = 1;
		for p in poly:


			#test for markers
			skip = False;
			if p.intersects(Point(11.24494,43.780094)):
				print boundary, ' has the Hotel'
				skip = True;

			if p.intersects(Point(11.248853,43.781198)):
				print boundary, ' has the Fortress'
				skip = True;

			if p.intersects(Point(11.247373,43.776893)):
				print boundary, ' has the Station'
				skip = True;

			if not skip:
				f.write("#exterior seg " + str(boundary - 1) + "\n")
				for i in range(1,len(p.exterior.coords)):
					if(i < len(p.exterior.coords) - 1):
						#segments are incremental 
						f.write(str(seg_cnt) + ' ' + str(seg_cnt) + ' ' + str(seg_cnt + 1) + ' ' + str(boundary) + '\n')
					else:
						#last segment loops back to first
						f.write(str(seg_cnt) + ' ' + str(seg_cnt) + ' ' + str(first_vert) + ' ' + str(boundary) + '\n')
					seg_cnt = seg_cnt + 1
				first_vert = seg_cnt
				boundary = boundary + 1

		for p in poly:
			for hole in p.interiors:
				f.write("#interior segment " + str(boundary - ext_cnt - 1) + "\n")
				#f.write("#interior segment " + str(boundary - len(poly.geoms) - 1) + "\n")
				for i in range(1,len(hole.coords)):
					if(i < len(hole.coords) -1):
						f.write(str(seg_cnt) + ' ' + str(seg_cnt) + ' ' + str(seg_cnt + 1) + ' ' + str(boundary) + '\n')
					else:
						f.write(str(seg_cnt) + ' ' + str(seg_cnt) + ' ' + str(first_vert) + ' ' + str(boundary) + '\n')
					seg_cnt = seg_cnt + 1
				first_vert = seg_cnt
				boundary = boundary + 1
	else:
		cnt = 0
		for i in range(1,len(poly.exterior.coords)):
			if(i < len(poly.exterior.coords) - 1):
				f.write(str(i + cnt) + ' ' + str(i+cnt) + ' ' + str(i+cnt+1) + ' 2\n')
			else:
				f.write(str(i + cnt) + ' ' + str(i+cnt) + ' ' + str(1+cnt) + ' 2\n')

		cnt = cnt + len(poly.exterior.coords) - 1;

		boundary = 3;
		for hole in poly.interiors:
			f.write("#interior segment " + str(boundary - 2) + '\n')
			for i in range(1,len(hole.coords)):
				if(i < len(hole.coords) - 1):
					f.write(str(i + cnt) + ' ' + str(i+cnt) + ' ' + str(i+cnt+1) +' ' + str(boundary) + '\n')
				else:
					f.write(str(i + cnt) + ' ' + str(i+cnt) + ' ' + str(1+cnt) + ' ' + str(boundary) + '\n')
			cnt = cnt + len(hole.coords) - 1
			boundary = boundary + 1

	f.write("# X holes\n")

	hole_cnt = 0
	if poly.geom_type == 'MultiPolygon':
		for p in poly:
			hole_cnt = hole_cnt + len(p.interiors)
	else:
		hole_cnt = len(poly.interiors)

	f.write(str(hole_cnt) + '\n')
	if hole_cnt > 0:
	
		hole_num = 1
		if poly.geom_type == 'MultiPolygon':
			for p in poly:
				for hole in p.interiors:
					test_poly = Polygon(hole)

					#test for markers
					if test_poly.intersects(Point(11.24494,43.780094)):
						print hole_num, ' has the Hotel'

					if test_poly.intersects(Point(11.248853,43.781198)):
						print hole_num, ' has the Fortress'

					if test_poly.intersects(Point(11.247373,43.776893)):
						print hole_num, ' has the Station'


					if test_poly.area > Hole_Area_Threshold:
						pt = btpoly_tri.draw_poly(hole.coords[1:]);
						x = pt[0][0] + .5*(pt[1][0]-pt[0][0]) + .25*(pt[2][0]-pt[0][0])
						y = pt[0][1] + .5*(pt[1][1]-pt[0][1]) + .25*(pt[2][1]-pt[0][1])
						f.write(str(hole_num) + ' ' + str(maplon_img(x)) + ' ' + str(maplat_img(y)) +'\n')
						if not (test_poly.contains(Point(x,y))):
							print "Point " + str(hole_num) + " is not interior?"
							print pt
							print x,y
						hole_num = hole_num + 1
					else:
						print "Skipped hole " + str(hole_num) + " is too small: " + str(test_poly.area)


		else:
			for hole in poly.interiors:
				test_poly = Polygon(hole)
				if test_poly.area > Hole_Area_Threshold:
					pt = btpoly_tri.draw_poly(hole.coords[1:]);
					x = pt[0][0] + .5*(pt[1][0]-pt[0][0]) + .25*(pt[2][0]-pt[0][0])
					y = pt[0][1] + .5*(pt[1][1]-pt[0][1]) + .25*(pt[2][1]-pt[0][1])
					f.write(str(hole_num) + ' ' + str(maplon_img(x)) + ' ' + str(maplat_img(y)) +'\n')
					if not (test_poly.contains(Point(x,y))):
						print "Point " + str(hole_num) + " is not interior?"
					hole_num = hole_num + 1
				else:
					print "Skipped hole " + str(hole_num) + " is too small: " + str(test_poly.area)


def poly2polyFile2(filename,addchar,poly=None,inverseRoads=False):

	Hole_Area_Threshold = 10e-20
	#hole = 0; #Need to fix this to deal w/ multiple holes
	params = json.load(open(directory + 'ParamFiles/' + filename + '.txt'))
	minlat = params['minlat']
	minlon = params['minlon']
	maxlat = params['maxlat']
	maxlon = params['maxlon']
	imgSize = params['imgSize']

	if poly is None:
		poly = roads2poly(filename);

	print minlat, minlon, maxlat, maxlon
	maplat_img = interp1d([minlat,maxlat],[imgSize,0])
	maplon_img = interp1d([minlon,maxlon],[0,imgSize])

	#Count the vertices
	vertices = 0
	if poly.geom_type == 'MultiPolygon':
		for p in poly:
			#polygon exteriors have repeated coordinates @ first,last so subtract 1
			vertices = vertices + len(p.exterior.coords) - 1
			for hole in p.interiors:
				vertices = vertices + len(hole.coords) - 1
	else:
		vertices = len(poly.exterior.coords) - 1
		for hole in poly.interiors:
			vertices = vertices + len(hole.coords) - 1

	segments = vertices

	if params['elevation']:
		if inverseRoads:
			#import testipoly
			vs = getElevationVertices(filename,poly)
			print len(vs)
			vertices = vertices + len(vs)
	else:
		vs = []

	#Create .poly file
	if inverseRoads:
		f = open(directory + 'PolyFiles/i_' + filename + '.poly','w')
	else:
		f = open(directory + 'PolyFiles/' + addchar + filename + '.poly','w')
	#Write Header
	f.write("#" + filename + ".poly\n")
	f.write(str(vertices) + ' 2 0 1\n')

	#Write the vertices
	if poly.geom_type == 'MultiPolygon':
		f.write("#the exteriors of the roads polygon\n")
		v_cnt = 1
		boundary = 2;
		#First write the exterior vertices
		for p in poly:
			f.write("#exterior " + str(boundary - 1) + "\n")
			for i in range(1,len(p.exterior.coords)):
				f.write(str(v_cnt) + ' ' + str(maplon_img(p.exterior.coords[i][0])) + ' ' + str(maplat_img(p.exterior.coords[i][1])) + ' ' + str(boundary) + '\n')
				v_cnt = v_cnt + 1
			boundary = boundary + 1
		f.write("#the end of the exterior polygons\n")

		#Next write the hole vertices
		for p in poly:
			for hole in p.interiors:
				f.write("#Hole " + str(boundary - len(poly.geoms) - 1) + "\n")
				for i in range(1,len(hole.coords)):
					f.write(str(v_cnt) + ' ' + str(maplon_img(hole.coords[i][0])) + ' ' + str(maplat_img(hole.coords[i][1])) + ' ' + str(boundary) + '\n')
					v_cnt = v_cnt + 1
				boundary = boundary + 1

	else:
		f.write("#the exterior of the roads polygon\n")

		v_cnt = 1
		boundary = 2
		for i in range(1,len(poly.exterior.coords)):
			f.write(str(v_cnt) + ' ' + str(maplon_img(poly.exterior.coords[i][0])) + ' ' + str(maplat_img(poly.exterior.coords[i][1])) + ' ' + str(boundary) + '\n')
			v_cnt = v_cnt + 1
		boundary = boundary + 1

		f.write("#the first hole\n");
		for hole in poly.interiors:
			f.write("#Hole 1\n")
			for i in range(1,len(hole.coords)):
				f.write(str(v_cnt) + ' ' + str(maplon_img(hole.coords[i][0])) + ' ' + str(maplat_img(hole.coords[i][1])) + ' ' + str(boundary) + '\n')
				v_cnt = v_cnt + 1
			boundary = boundary + 1

	if inverseRoads:
		f.write("#Sample Points\n")
		for i in range(0,len(vs)):
			f.write(str(v_cnt) + ' ' + str(maplon_img(vs[i][0])) + ' ' + str(maplat_img(vs[i][1])) + ' 0\n')
			v_cnt = v_cnt + 1

	#write segment header
	f.write("#the segments\n")
	f.write(str(segments) + ' 1\n')
	f.write("#exterior segment\n")

	if poly.geom_type == 'MultiPolygon':
		boundary = 2;
		seg_cnt = 1;
		first_vert = 1;
		for p in poly:
			f.write("#exterior seg " + str(boundary - 1) + "\n")
			for i in range(1,len(p.exterior.coords)):
				if(i < len(p.exterior.coords) - 1):
					#segments are incremental 
					f.write(str(seg_cnt) + ' ' + str(seg_cnt) + ' ' + str(seg_cnt + 1) + ' ' + str(boundary) + '\n')
				else:
					#last segment loops back to first
					f.write(str(seg_cnt) + ' ' + str(seg_cnt) + ' ' + str(first_vert) + ' ' + str(boundary) + '\n')
				seg_cnt = seg_cnt + 1
			first_vert = seg_cnt
			boundary = boundary + 1

		for p in poly:
			for hole in p.interiors:
				f.write("#interior segment " + str(boundary - len(poly.geoms) - 1) + "\n")
				for i in range(1,len(hole.coords)):
					if(i < len(hole.coords) -1):
						f.write(str(seg_cnt) + ' ' + str(seg_cnt) + ' ' + str(seg_cnt + 1) + ' ' + str(boundary) + '\n')
					else:
						f.write(str(seg_cnt) + ' ' + str(seg_cnt) + ' ' + str(first_vert) + ' ' + str(boundary) + '\n')
					seg_cnt = seg_cnt + 1
				first_vert = seg_cnt
				boundary = boundary + 1
	else:
		cnt = 0
		for i in range(1,len(poly.exterior.coords)):
			if(i < len(poly.exterior.coords) - 1):
				f.write(str(i + cnt) + ' ' + str(i+cnt) + ' ' + str(i+cnt+1) + ' 2\n')
			else:
				f.write(str(i + cnt) + ' ' + str(i+cnt) + ' ' + str(1+cnt) + ' 2\n')

		cnt = cnt + len(poly.exterior.coords) - 1;

		boundary = 3;
		for hole in poly.interiors:
			f.write("#interior segment " + str(boundary - 2) + '\n')
			for i in range(1,len(hole.coords)):
				if(i < len(hole.coords) - 1):
					f.write(str(i + cnt) + ' ' + str(i+cnt) + ' ' + str(i+cnt+1) +' ' + str(boundary) + '\n')
				else:
					f.write(str(i + cnt) + ' ' + str(i+cnt) + ' ' + str(1+cnt) + ' ' + str(boundary) + '\n')
			cnt = cnt + len(hole.coords) - 1
			boundary = boundary + 1

	f.write("# X holes\n")

	hole_cnt = 0
	if poly.geom_type == 'MultiPolygon':
		for p in poly:
			hole_cnt = hole_cnt + len(p.interiors)
	else:
		hole_cnt = len(poly.interiors)

	f.write(str(hole_cnt) + '\n')
	if hole_cnt > 0:
	
		hole_num = 1
		if poly.geom_type == 'MultiPolygon':
			for p in poly:
				for hole in p.interiors:
					test_poly = Polygon(hole)
					if test_poly.area > Hole_Area_Threshold:
						pt = btpoly_tri.draw_poly(hole.coords[1:]);
						x = pt[0][0] + .5*(pt[1][0]-pt[0][0]) + .25*(pt[2][0]-pt[0][0])
						y = pt[0][1] + .5*(pt[1][1]-pt[0][1]) + .25*(pt[2][1]-pt[0][1])
						f.write(str(hole_num) + ' ' + str(maplon_img(x)) + ' ' + str(maplat_img(y)) +'\n')
						if not (test_poly.contains(Point(x,y))):
							print "Point " + str(hole_num) + " is not interior?"
							print pt
							print x,y
						hole_num = hole_num + 1
					else:
						print "Skipped hole " + str(hole_num) + " is too small: " + str(test_poly.area)


		else:
			for hole in poly.interiors:
				test_poly = Polygon(hole)
				if test_poly.area > Hole_Area_Threshold:
					pt = btpoly_tri.draw_poly(hole.coords[1:]);
					x = pt[0][0] + .5*(pt[1][0]-pt[0][0]) + .25*(pt[2][0]-pt[0][0])
					y = pt[0][1] + .5*(pt[1][1]-pt[0][1]) + .25*(pt[2][1]-pt[0][1])
					f.write(str(hole_num) + ' ' + str(maplon_img(x)) + ' ' + str(maplat_img(y)) +'\n')
					if not (test_poly.contains(Point(x,y))):
						print "Point " + str(hole_num) + " is not interior?"
					hole_num = hole_num + 1
				else:
					print "Skipped hole " + str(hole_num) + " is too small: " + str(test_poly.area)



def polyFile2stl(filename,poly=None):
	global height, i2e, samples
	params = json.load(open(directory + 'ParamFiles/' + filename + '.txt'))
	minlat = params['minlat']
	minlon = params['minlon']
	maxlat = params['maxlat']
	maxlon = params['maxlon']
	imgSize = params['imgSize']
	#global minlon, minlat, maxlon, maxlat
	import triangle
	import triangle.plot
	#import matplotlib.pyplot as plt
	import elevation2stl
	#import testipoly
	import numpy as np
	samples = params['samples']
	road_offset = params['road_height']
	print_size = params['print_size']

	if poly is None:
		poly = roads2poly(filename)
		poly2polyFile(filename,poly,False);

		grid = box(minlon,minlat,maxlon,maxlat)
		ipoly = grid.difference(poly)
		poly2polyFile(filename,ipoly,True)

	if params['elevation']:
		print 'elevation is true?'
		#elevs = testipoly.getElevations(filename)
		elevs = getElevations(filename)
	else:
		print 'flat elevation?'
		elevs = np.zeros((2,2))

	z_scale = float(params['z_scale'])
	if params['zmin']:
		print 'zmin is true'
		height = interp2d(range(0,len(elevs)),range(0,len(elevs[0])),z_scale*(elevs - float(zminVal)));
	else:
		print 'zmin is false'
		print elevs.min(), elevs.max()
		print 'z_scale: ' + str(z_scale)
		print len(elevs), len(elevs[0])
		v = z_scale*(elevs - elevs.min());
		print elevs.shape
		height = interp2d(range(0,len(elevs)),range(0,len(elevs[0])),v);
	print elevs.min(), elevs.max()
	floor = 0
	base_height = params['base_height']
	#i2e = interp1d([0,imgSize],[0,samples-1])
	i2e = interp1d([0,imgSize],[0,samples-1])
	img2printx = interp1d([0,imgSize],[0,print_size])
	img2printy = interp1d([0,imgSize],[print_size,0])

	edges = box(0,0,imgSize,imgSize)
	edge_vs = []
	print 'gonna get the polyfiles'
	#poly2polyFile(filename);
	RoadPSLG = triangle.get_data(params['dir'] + 'PolyFiles/' + filename)
	RoadCDT = triangle.triangulate(RoadPSLG,'p')
	nonRoadPSLG = triangle.get_data(params['dir'] + 'PolyFiles/' + 'i_' + filename)
	nonRoadCDT = triangle.triangulate(nonRoadPSLG,'p')
	print 'gonna write the stl header'
	elevation2stl.startSTL(filename)
	for tri in RoadCDT['triangles']:
		#top layer of roads
		pt1 = Point(RoadCDT['vertices'][tri[0]][0],RoadCDT['vertices'][tri[0]][1])
		pt2 = Point(RoadCDT['vertices'][tri[1]][0],RoadCDT['vertices'][tri[1]][1])
		pt3 = Point(RoadCDT['vertices'][tri[2]][0],RoadCDT['vertices'][tri[2]][1])
		#points are scaled to imgSize
		tri_v1 = [img2printx(pt1.x),img2printy(pt1.y),get_height(pt1) + road_offset + base_height]
		tri_v2 = [img2printx(pt2.x),img2printy(pt2.y),get_height(pt2) + road_offset + base_height]
		tri_v3 = [img2printx(pt3.x),img2printy(pt3.y),get_height(pt3) + road_offset + base_height]
		elevation2stl.printTriangle(tri_v1,tri_v2,tri_v3,filename)
	for tri in nonRoadCDT['triangles']:
		#top layer of non-roads
		pt1 = Point(nonRoadCDT['vertices'][tri[0]][0],nonRoadCDT['vertices'][tri[0]][1])
		pt2 = Point(nonRoadCDT['vertices'][tri[1]][0],nonRoadCDT['vertices'][tri[1]][1])
		pt3 = Point(nonRoadCDT['vertices'][tri[2]][0],nonRoadCDT['vertices'][tri[2]][1])
		#points are scaled to imgSize
		tri_v1 = [img2printx(pt1.x),img2printy(pt1.y),get_height(pt1) + base_height]
		tri_v2 = [img2printx(pt2.x),img2printy(pt2.y),get_height(pt2) + base_height]
		tri_v3 = [img2printx(pt3.x),img2printy(pt3.y),get_height(pt3) + base_height]
		elevation2stl.printTriangle(tri_v1,tri_v2,tri_v3,filename)
	for seg in RoadCDT['segments']:
		#build the sides of the roads
		line = LineString(RoadCDT['vertices'][seg])
		if line.touches(edges):
			#along the edges of the print, fill from the base_height to the floor
			#fill road edge sides
			XYpt1 = Point(RoadCDT['vertices'][seg[0]])
			XYpt2 = Point(RoadCDT['vertices'][seg[1]])
			edge_vs.append(XYpt1.coords[0])
			edge_vs.append(XYpt2.coords[0])
			pt1 = [img2printx(XYpt1.x),img2printy(XYpt1.y),get_height(XYpt1) + base_height]
			pt2 = [img2printx(XYpt2.x),img2printy(XYpt2.y),get_height(XYpt2) + base_height]
			pt3 = [img2printx(XYpt1.x),img2printy(XYpt1.y),floor]
			pt4 = [img2printx(XYpt2.x),img2printy(XYpt2.y),floor]
			elevation2stl.printTriangle(pt1,pt2,pt3,filename)
			elevation2stl.printTriangle(pt3,pt2,pt4,filename)
		
		#fill all road sides from the road_offset to the base_height
		#this will cover the gap from base to road height along the edges
		XYpt1 = Point(RoadCDT['vertices'][seg[0]])
		XYpt2 = Point(RoadCDT['vertices'][seg[1]])
		pt1 = [img2printx(XYpt1.x),img2printy(XYpt1.y),get_height(XYpt1) + road_offset + base_height]
		pt2 = [img2printx(XYpt2.x),img2printy(XYpt2.y),get_height(XYpt2) + road_offset + base_height]
		pt3 = [img2printx(XYpt1.x),img2printy(XYpt1.y),get_height(XYpt1) + base_height]
		pt4 = [img2printx(XYpt2.x),img2printy(XYpt2.y),get_height(XYpt2) + base_height]
		elevation2stl.printTriangle(pt1,pt2,pt3,filename)
		elevation2stl.printTriangle(pt3,pt2,pt4,filename)

	for seg in nonRoadCDT['segments']:
		line = LineString(nonRoadCDT['vertices'][seg])
		#fill non-road edge sides
		if line.touches(edges):
			XYpt1 = Point(nonRoadCDT['vertices'][seg[0]])
			XYpt2 = Point(nonRoadCDT['vertices'][seg[1]])
			edge_vs.append(XYpt1.coords[0])
			edge_vs.append(XYpt2.coords[0])
			pt1 = [img2printx(XYpt1.x),img2printy(XYpt1.y),get_height(XYpt1) + base_height]
			pt2 = [img2printx(XYpt2.x),img2printy(XYpt2.y),get_height(XYpt2) + base_height]
			pt3 = [img2printx(XYpt1.x),img2printy(XYpt1.y),floor]
			pt4 = [img2printx(XYpt2.x),img2printy(XYpt2.y),floor]
			elevation2stl.printTriangle(pt1,pt2,pt3,filename)
			elevation2stl.printTriangle(pt3,pt2,pt4,filename)

	for i in range(1,13):
		print 'PolyFiles/f_' + str(i) + '_' + filename
		WaterPSLG = triangle.get_data(params['dir'] + 'PolyFiles/f_' + str(i) + '_' + filename)
		if(len(WaterPSLG['vertices']) > 3):
			WaterCDT = triangle.triangulate(RoadPSLG,'p')
			for tri in WaterCDT['triangles']:
				#top layer of roads
				pt1 = Point(WaterCDT['vertices'][tri[0]][0],WaterCDT['vertices'][tri[0]][1])
				pt2 = Point(WaterCDT['vertices'][tri[1]][0],WaterCDT['vertices'][tri[1]][1])
				pt3 = Point(WaterCDT['vertices'][tri[2]][0],WaterCDT['vertices'][tri[2]][1])
				#points are scaled to imgSize
				tri_v1 = [img2printx(pt1.x),img2printy(pt1.y),get_height(pt1) + 10 + base_height]
				tri_v2 = [img2printx(pt2.x),img2printy(pt2.y),get_height(pt2) + 10 + base_height]
				tri_v3 = [img2printx(pt3.x),img2printy(pt3.y),get_height(pt3) + 10 + base_height]
				elevation2stl.printTriangle(tri_v1,tri_v2,tri_v3,filename)

	#build the base
	edge_vs.append([0,0])
	edge_vs.append([0,imgSize])
	edge_vs.append([imgSize,0])
	edge_vs.append([imgSize,imgSize])
	edge_array = np.array(edge_vs)
	base = dict(vertices=edge_array)
	baseT = triangle.triangulate(base)

	for tri in baseT['triangles']:
		pt1 = Point(baseT['vertices'][tri[0]][0],baseT['vertices'][tri[0]][1])
		pt2 = Point(baseT['vertices'][tri[1]][0],baseT['vertices'][tri[1]][1])
		pt3 = Point(baseT['vertices'][tri[2]][0],baseT['vertices'][tri[2]][1])
		tri_v1 = [img2printx(pt1.x),img2printy(pt1.y),floor]
		tri_v2 = [img2printx(pt2.x),img2printy(pt2.y),floor]
		tri_v3 = [img2printx(pt3.x),img2printy(pt3.y),floor]

		elevation2stl.printTriangle(tri_v1,tri_v2,tri_v3,filename)



	elevation2stl.endSTL(filename)
	return edge_vs



def polyFile2stlBin(filename,poly=None):
	global height, i2e, samples
	facetCnt = 0
	params = json.load(open(directory + 'ParamFiles/' + filename + '.txt'))
	minlat = params['minlat']
	minlon = params['minlon']
	maxlat = params['maxlat']
	maxlon = params['maxlon']
	imgSize = params['imgSize']
	#global minlon, minlat, maxlon, maxlat
	import triangle
	import triangle.plot
	#import matplotlib.pyplot as plt
	import elevation2stl
	#import testipoly
	import numpy as np
	samples = params['samples']
	road_offset = params['road_height']
	print_size = params['print_size']

	if poly is None:
		poly = roads2poly(filename)
		poly2polyFile(filename,poly,False);

		grid = box(minlon,minlat,maxlon,maxlat)
		ipoly = grid.difference(poly)
		poly2polyFile(filename,ipoly,True)

	if params['elevation']:
		print 'elevation is true?'
		#elevs = testipoly.getElevations(filename)
		elevs = getElevations(filename)
	else:
		print 'flat elevation?'
		elevs = np.zeros((2,2))

	z_scale = float(params['z_scale'])
	if params['zmin']:
		print 'zmin is true'
		height = interp2d(range(0,len(elevs)),range(0,len(elevs[0])),z_scale*(elevs - float(zminVal)));
	else:
		print 'zmin is false'
		print elevs.min(), elevs.max()
		print 'z_scale: ' + str(z_scale)
		print len(elevs), len(elevs[0])
		v = z_scale*(elevs - elevs.min());
		print elevs.shape
		height = interp2d(range(0,len(elevs)),range(0,len(elevs[0])),v);
	print elevs.min(), elevs.max()
	floor = 0
	base_height = params['base_height']
	#i2e = interp1d([0,imgSize],[0,samples-1])
	i2e = interp1d([0,imgSize],[0,samples-1])
	img2printx = interp1d([0,imgSize],[0,print_size])
	img2printy = interp1d([0,imgSize],[print_size,0])

	edges = box(0,0,imgSize,imgSize)
	edge_vs = []
	print 'gonna get the polyfiles'
	#poly2polyFile(filename);
	RoadPSLG = triangle.get_data(params['dir'] + 'PolyFiles/' + filename)
	RoadCDT = triangle.triangulate(RoadPSLG,'p')
	nonRoadPSLG = triangle.get_data(params['dir'] + 'PolyFiles/' + 'i_' + filename)
	nonRoadCDT = triangle.triangulate(nonRoadPSLG,'p')
	print 'gonna write the stl header'
	elevation2stl.startSTLBin(filename)

	for tri in RoadCDT['triangles']:
		#top layer of roads
		pt1 = Point(RoadCDT['vertices'][tri[0]][0],RoadCDT['vertices'][tri[0]][1])
		pt2 = Point(RoadCDT['vertices'][tri[1]][0],RoadCDT['vertices'][tri[1]][1])
		pt3 = Point(RoadCDT['vertices'][tri[2]][0],RoadCDT['vertices'][tri[2]][1])
		#points are scaled to imgSize
		tri_v1 = [img2printx(pt1.x),img2printy(pt1.y),get_height(pt1) + road_offset + base_height]
		tri_v2 = [img2printx(pt2.x),img2printy(pt2.y),get_height(pt2) + road_offset + base_height]
		tri_v3 = [img2printx(pt3.x),img2printy(pt3.y),get_height(pt3) + road_offset + base_height]
		elevation2stl.printTriangleBin(tri_v1,tri_v2,tri_v3,filename)
		facetCnt = facetCnt + 1
	for tri in nonRoadCDT['triangles']:
		#top layer of non-roads
		pt1 = Point(nonRoadCDT['vertices'][tri[0]][0],nonRoadCDT['vertices'][tri[0]][1])
		pt2 = Point(nonRoadCDT['vertices'][tri[1]][0],nonRoadCDT['vertices'][tri[1]][1])
		pt3 = Point(nonRoadCDT['vertices'][tri[2]][0],nonRoadCDT['vertices'][tri[2]][1])
		#points are scaled to imgSize
		tri_v1 = [img2printx(pt1.x),img2printy(pt1.y),get_height(pt1) + base_height]
		tri_v2 = [img2printx(pt2.x),img2printy(pt2.y),get_height(pt2) + base_height]
		tri_v3 = [img2printx(pt3.x),img2printy(pt3.y),get_height(pt3) + base_height]
		elevation2stl.printTriangleBin(tri_v1,tri_v2,tri_v3,filename)
		facetCnt = facetCnt + 1
	for seg in RoadCDT['segments']:
		#build the sides of the roads
		line = LineString(RoadCDT['vertices'][seg])
		if line.touches(edges):
			#along the edges of the print, fill from the base_height to the floor
			#fill road edge sides
			XYpt1 = Point(RoadCDT['vertices'][seg[0]])
			XYpt2 = Point(RoadCDT['vertices'][seg[1]])
			edge_vs.append(XYpt1.coords[0])
			edge_vs.append(XYpt2.coords[0])
			pt1 = [img2printx(XYpt1.x),img2printy(XYpt1.y),get_height(XYpt1) + base_height]
			pt2 = [img2printx(XYpt2.x),img2printy(XYpt2.y),get_height(XYpt2) + base_height]
			pt3 = [img2printx(XYpt1.x),img2printy(XYpt1.y),floor]
			pt4 = [img2printx(XYpt2.x),img2printy(XYpt2.y),floor]
			elevation2stl.printTriangleBin(pt1,pt2,pt3,filename)
			facetCnt = facetCnt + 1
			elevation2stl.printTriangleBin(pt3,pt2,pt4,filename)
			facetCnt = facetCnt + 1
		
		#fill all road sides from the road_offset to the base_height
		#this will cover the gap from base to road height along the edges
		XYpt1 = Point(RoadCDT['vertices'][seg[0]])
		XYpt2 = Point(RoadCDT['vertices'][seg[1]])
		pt1 = [img2printx(XYpt1.x),img2printy(XYpt1.y),get_height(XYpt1) + road_offset + base_height]
		pt2 = [img2printx(XYpt2.x),img2printy(XYpt2.y),get_height(XYpt2) + road_offset + base_height]
		pt3 = [img2printx(XYpt1.x),img2printy(XYpt1.y),get_height(XYpt1) + base_height]
		pt4 = [img2printx(XYpt2.x),img2printy(XYpt2.y),get_height(XYpt2) + base_height]
		elevation2stl.printTriangleBin(pt1,pt2,pt3,filename)
		facetCnt = facetCnt + 1
		elevation2stl.printTriangleBin(pt3,pt2,pt4,filename)
		facetCnt = facetCnt + 1

	for seg in nonRoadCDT['segments']:
		line = LineString(nonRoadCDT['vertices'][seg])
		#fill non-road edge sides
		if line.touches(edges):
			XYpt1 = Point(nonRoadCDT['vertices'][seg[0]])
			XYpt2 = Point(nonRoadCDT['vertices'][seg[1]])
			edge_vs.append(XYpt1.coords[0])
			edge_vs.append(XYpt2.coords[0])
			pt1 = [img2printx(XYpt1.x),img2printy(XYpt1.y),get_height(XYpt1) + base_height]
			pt2 = [img2printx(XYpt2.x),img2printy(XYpt2.y),get_height(XYpt2) + base_height]
			pt3 = [img2printx(XYpt1.x),img2printy(XYpt1.y),floor]
			pt4 = [img2printx(XYpt2.x),img2printy(XYpt2.y),floor]
			elevation2stl.printTriangleBin(pt1,pt2,pt3,filename)
			facetCnt = facetCnt + 1
			elevation2stl.printTriangleBin(pt3,pt2,pt4,filename)
			facetCnt = facetCnt + 1


	# #WaterPSLG = triangle.get_data(params['dir'] + 'PolyFiles/f_' + filename)
	# for i in range(1,13):
	# 	print 'PolyFiles/f_' + str(i) + '_' + filename
	# 	WaterPSLG = triangle.get_data(params['dir'] + 'PolyFiles/f_' + str(i) + '_' + filename)
	# 	if(len(WaterPSLG['vertices']) > 3):
	# 		WaterCDT = triangle.triangulate(WaterPSLG,'p')
	# 		for tri in WaterCDT['triangles']:
	# 			#top layer of roads
	# 			pt1 = Point(WaterCDT['vertices'][tri[0]][0],WaterCDT['vertices'][tri[0]][1])
	# 			pt2 = Point(WaterCDT['vertices'][tri[1]][0],WaterCDT['vertices'][tri[1]][1])
	# 			pt3 = Point(WaterCDT['vertices'][tri[2]][0],WaterCDT['vertices'][tri[2]][1])
	# 			#points are scaled to imgSize
	# 			tri_v1 = [img2printx(pt1.x),img2printy(pt1.y),get_height(pt1) + 10 + base_height]
	# 			tri_v2 = [img2printx(pt2.x),img2printy(pt2.y),get_height(pt2) + 10 + base_height]
	# 			tri_v3 = [img2printx(pt3.x),img2printy(pt3.y),get_height(pt3) + 10 + base_height]
	# 			elevation2stl.printTriangleBin(tri_v1,tri_v2,tri_v3,filename)


	#build the base
	edge_vs.append([0,0])
	edge_vs.append([0,imgSize])
	edge_vs.append([imgSize,0])
	edge_vs.append([imgSize,imgSize])
	edge_array = np.array(edge_vs)
	base = dict(vertices=edge_array)
	baseT = triangle.triangulate(base)

	for tri in baseT['triangles']:
		pt1 = Point(baseT['vertices'][tri[0]][0],baseT['vertices'][tri[0]][1])
		pt2 = Point(baseT['vertices'][tri[1]][0],baseT['vertices'][tri[1]][1])
		pt3 = Point(baseT['vertices'][tri[2]][0],baseT['vertices'][tri[2]][1])
		tri_v1 = [img2printx(pt1.x),img2printy(pt1.y),floor]
		tri_v2 = [img2printx(pt2.x),img2printy(pt2.y),floor]
		tri_v3 = [img2printx(pt3.x),img2printy(pt3.y),floor]

		elevation2stl.printTriangleBin(tri_v1,tri_v2,tri_v3,filename)
		facetCnt = facetCnt + 1

	elevation2stl.endSTLBin(filename,facetCnt)
	return edge_vs



def polyFile2stlBin2(filename,poly,feat_tris,feat_edges):
	global height, i2e, samples
	facetCnt = 0
	params = json.load(open(directory + 'ParamFiles/' + filename + '.txt'))
	minlat = params['minlat']
	minlon = params['minlon']
	maxlat = params['maxlat']
	maxlon = params['maxlon']
	imgSize = params['imgSize']
	#global minlon, minlat, maxlon, maxlat
	import triangle
	import triangle.plot
	#import matplotlib.pyplot as plt
	import elevation2stl
	#import testipoly
	import numpy as np
	samples = params['samples']
	road_offset = params['road_height']
	print_size = params['print_size']

	if poly is None:
		poly = roads2poly(filename)
		poly2polyFile(filename,poly,False);

		grid = box(minlon,minlat,maxlon,maxlat)
		ipoly = grid.difference(poly)
		poly2polyFile(filename,ipoly,True)

	if params['elevation']:
		print 'elevation is true?'
		#elevs = testipoly.getElevations(filename)
		elevs = getElevations(filename)
	else:
		print 'flat elevation?'
		elevs = np.zeros((2,2))

	z_scale = float(params['z_scale'])
	if params['zmin']:
		print 'zmin is true'
		height = interp2d(range(0,len(elevs)),range(0,len(elevs[0])),z_scale*(elevs - float(zminVal)));
	else:
		print 'zmin is false'
		print elevs.min(), elevs.max()
		print 'z_scale: ' + str(z_scale)
		print len(elevs), len(elevs[0])
		v = z_scale*(elevs - elevs.min());
		print elevs.shape
		height = interp2d(range(0,len(elevs)),range(0,len(elevs[0])),v);
	print elevs.min(), elevs.max()
	floor = 0
	base_height = params['base_height']
	#i2e = interp1d([0,imgSize],[0,samples-1])
	i2e = interp1d([0,imgSize],[0,samples-1])
	img2printx = interp1d([0,imgSize],[0,print_size])
	img2printy = interp1d([0,imgSize],[print_size,0])

	edges = box(0,0,imgSize,imgSize)
	edge_vs = []
	print 'gonna get the polyfiles'
	#poly2polyFile(filename);
	RoadPSLG = triangle.get_data(params['dir'] + 'PolyFiles/' + filename)
	RoadCDT = triangle.triangulate(RoadPSLG,'p')
	nonRoadPSLG = triangle.get_data(params['dir'] + 'PolyFiles/' + 'i_' + filename)
	nonRoadCDT = triangle.triangulate(nonRoadPSLG,'p')
	print 'gonna write the stl header'
	elevation2stl.startSTLBin(filename)

	for feat_tri in feat_tris:
		for tri in feat_tri:
			elevation2stl.printTriangleBin(tri[0],tri[1],tri[2],filename)
			facetCnt = facetCnt + 1

	for tri in RoadCDT['triangles']:
		#top layer of roads
		pt1 = Point(RoadCDT['vertices'][tri[0]][0],RoadCDT['vertices'][tri[0]][1])
		pt2 = Point(RoadCDT['vertices'][tri[1]][0],RoadCDT['vertices'][tri[1]][1])
		pt3 = Point(RoadCDT['vertices'][tri[2]][0],RoadCDT['vertices'][tri[2]][1])
		#points are scaled to imgSize
		tri_v1 = [img2printx(pt1.x),img2printy(pt1.y),get_height(pt1) + road_offset + base_height]
		tri_v2 = [img2printx(pt2.x),img2printy(pt2.y),get_height(pt2) + road_offset + base_height]
		tri_v3 = [img2printx(pt3.x),img2printy(pt3.y),get_height(pt3) + road_offset + base_height]
		elevation2stl.printTriangleBin(tri_v1,tri_v2,tri_v3,filename)
		facetCnt = facetCnt + 1
	for tri in nonRoadCDT['triangles']:
		#top layer of non-roads
		pt1 = Point(nonRoadCDT['vertices'][tri[0]][0],nonRoadCDT['vertices'][tri[0]][1])
		pt2 = Point(nonRoadCDT['vertices'][tri[1]][0],nonRoadCDT['vertices'][tri[1]][1])
		pt3 = Point(nonRoadCDT['vertices'][tri[2]][0],nonRoadCDT['vertices'][tri[2]][1])
		#points are scaled to imgSize
		tri_v1 = [img2printx(pt1.x),img2printy(pt1.y),get_height(pt1) + base_height]
		tri_v2 = [img2printx(pt2.x),img2printy(pt2.y),get_height(pt2) + base_height]
		tri_v3 = [img2printx(pt3.x),img2printy(pt3.y),get_height(pt3) + base_height]
		elevation2stl.printTriangleBin(tri_v1,tri_v2,tri_v3,filename)
		facetCnt = facetCnt + 1
	for seg in RoadCDT['segments']:
		#build the sides of the roads
		line = LineString(RoadCDT['vertices'][seg])
		if line.touches(edges):
			#along the edges of the print, fill from the base_height to the floor
			#fill road edge sides
			XYpt1 = Point(RoadCDT['vertices'][seg[0]])
			XYpt2 = Point(RoadCDT['vertices'][seg[1]])
			edge_vs.append(XYpt1.coords[0])
			edge_vs.append(XYpt2.coords[0])
			pt1 = [img2printx(XYpt1.x),img2printy(XYpt1.y),get_height(XYpt1) + base_height]
			pt2 = [img2printx(XYpt2.x),img2printy(XYpt2.y),get_height(XYpt2) + base_height]
			pt3 = [img2printx(XYpt1.x),img2printy(XYpt1.y),floor]
			pt4 = [img2printx(XYpt2.x),img2printy(XYpt2.y),floor]
			elevation2stl.printTriangleBin(pt1,pt2,pt3,filename)
			facetCnt = facetCnt + 1
			elevation2stl.printTriangleBin(pt3,pt2,pt4,filename)
			facetCnt = facetCnt + 1
		
		#fill all road sides from the road_offset to the base_height
		#this will cover the gap from base to road height along the edges
		XYpt1 = Point(RoadCDT['vertices'][seg[0]])
		XYpt2 = Point(RoadCDT['vertices'][seg[1]])
		pt1 = [img2printx(XYpt1.x),img2printy(XYpt1.y),get_height(XYpt1) + road_offset + base_height]
		pt2 = [img2printx(XYpt2.x),img2printy(XYpt2.y),get_height(XYpt2) + road_offset + base_height]
		pt3 = [img2printx(XYpt1.x),img2printy(XYpt1.y),get_height(XYpt1) + base_height]
		pt4 = [img2printx(XYpt2.x),img2printy(XYpt2.y),get_height(XYpt2) + base_height]
		elevation2stl.printTriangleBin(pt1,pt2,pt3,filename)
		facetCnt = facetCnt + 1
		elevation2stl.printTriangleBin(pt3,pt2,pt4,filename)
		facetCnt = facetCnt + 1

	for seg in nonRoadCDT['segments']:
		line = LineString(nonRoadCDT['vertices'][seg])
		#fill non-road edge sides
		if line.touches(edges):
			XYpt1 = Point(nonRoadCDT['vertices'][seg[0]])
			XYpt2 = Point(nonRoadCDT['vertices'][seg[1]])
			edge_vs.append(XYpt1.coords[0])
			edge_vs.append(XYpt2.coords[0])
			pt1 = [img2printx(XYpt1.x),img2printy(XYpt1.y),get_height(XYpt1) + base_height]
			pt2 = [img2printx(XYpt2.x),img2printy(XYpt2.y),get_height(XYpt2) + base_height]
			pt3 = [img2printx(XYpt1.x),img2printy(XYpt1.y),floor]
			pt4 = [img2printx(XYpt2.x),img2printy(XYpt2.y),floor]
			elevation2stl.printTriangleBin(pt1,pt2,pt3,filename)
			facetCnt = facetCnt + 1
			elevation2stl.printTriangleBin(pt3,pt2,pt4,filename)
			facetCnt = facetCnt + 1


	#WaterPSLG = triangle.get_data(params['dir'] + 'PolyFiles/f_' + filename)
	
	for i in range(1,13):
		print 'PolyFiles/f_' + str(i) + '_' + filename
		WaterPSLG = triangle.get_data(params['dir'] + 'PolyFiles/f_' + str(i) + '_' + filename)
		if(len(WaterPSLG['vertices']) > 3):
			WaterCDT = triangle.triangulate(WaterPSLG,'p')
			if 'triangles' in WaterCDT.keys():
				for tri in WaterCDT['triangles']:
					#top layer of roads
					pt1 = Point(WaterCDT['vertices'][tri[0]][0],WaterCDT['vertices'][tri[0]][1])
					pt2 = Point(WaterCDT['vertices'][tri[1]][0],WaterCDT['vertices'][tri[1]][1])
					pt3 = Point(WaterCDT['vertices'][tri[2]][0],WaterCDT['vertices'][tri[2]][1])
					#points are scaled to imgSize
					tri_v1 = [img2printx(pt1.x),img2printy(pt1.y),get_height(pt1)  + base_height]
					tri_v2 = [img2printx(pt2.x),img2printy(pt2.y),get_height(pt2)  + base_height]
					tri_v3 = [img2printx(pt3.x),img2printy(pt3.y),get_height(pt3)  + base_height]
					elevation2stl.printTriangleBin(tri_v1,tri_v2,tri_v3,filename)
					facetCnt = facetCnt + 1
	

	#build the base
	edge_vs.append([0,0])
	edge_vs.append([0,imgSize])
	edge_vs.append([imgSize,0])
	edge_vs.append([imgSize,imgSize])
	edge_array = np.array(edge_vs)
	base = dict(vertices=edge_array)
	baseT = triangle.triangulate(base)

	for tri in baseT['triangles']:
		pt1 = Point(baseT['vertices'][tri[0]][0],baseT['vertices'][tri[0]][1])
		pt2 = Point(baseT['vertices'][tri[1]][0],baseT['vertices'][tri[1]][1])
		pt3 = Point(baseT['vertices'][tri[2]][0],baseT['vertices'][tri[2]][1])
		tri_v1 = [img2printx(pt1.x),img2printy(pt1.y),floor]
		tri_v2 = [img2printx(pt2.x),img2printy(pt2.y),floor]
		tri_v3 = [img2printx(pt3.x),img2printy(pt3.y),floor]

		elevation2stl.printTriangleBin(tri_v1,tri_v2,tri_v3,filename)
		facetCnt = facetCnt + 1

	# = facetCnt 
	elevation2stl.endSTLBin(filename,facetCnt)
	return edge_vs



def get_height(pt):#
	global height, i2e, samples
	return height((samples-1) - i2e(pt.x),i2e(pt.y))[0]

def getElevationVertices(filename,poly):
	import csv
	import numpy as np

	params = json.load(open(directory + 'ParamFiles/' + filename + '.txt'))
	dataList = []
	latList = []
	lonList = []

	samples = params['samples']
	#maplat_img = interp1d([minlat,maxlat],[imgSize,0])
	#maplon_img = interp1d([minlon,maxlon],[0,imgSize])

	#filename = 'thePoint'

  	with open(directory + 'elevationFiles/' + filename + '.csv','rb') as csvfile:
  		csvreader = csv.reader(csvfile,delimiter=',')
  		for row in csvreader:
  			#track lat of start and end of each row 
  			latList.append([float(row[1]),float(row[3])])
  			lonList.append([float(row[2]),float(row[4])])
  			dataList.append(row[5:]);

  	X = len(dataList)		#number of lines in .csv
  	Y = len(dataList[1])	#number of entries per line
  	print X,Y

  	lats = np.zeros((X,Y))
  	lons = np.zeros((X,Y))
  	height = np.zeros((X,Y))

  	vertices =[]
  	
  	inner_cnt = 0;
  	for i in range(0,X):
  		for j in range(0,Y):
  			lats[i,j] = latList[i][0] + (float(i)*(latList[i][0] - latList[i][1]))/(samples - 1)
  			lons[i,j] = lonList[j][0] - (float(j)*(lonList[j][0] - lonList[j][1]))/(samples - 1)
  			height[i,j] = float(dataList[i][j])
  			p = Point(lons[i,j],lats[i,j])
  			#p = Point(maplon_img(lons[i,j]),maplat_img(lats[i,j]))
  			#if p.within(poly):
  			if p.intersects(poly):
  				#vertices.append([p.x,p.y])
  				vertices.append([lons[i,j],lats[i,j]])
  				inner_cnt = inner_cnt + 1

  	print "inner_cnt: " + str(inner_cnt)
  	print "lats:"
  	#print lats[0,0], lats[0,10]
  	print lats[0,1], lats[10,0], lats[19,0]
  	print "lons: "
  	print lons[0,0], lons[0,10], lons[0,19]
  	#print lons[0,0], lons[10,0]
  	return vertices

def getElevations(filename):
	import csv
	import numpy as np
	dataList = []
  	with open(directory + 'elevationFiles/' + filename + '.csv','rb') as csvfile:
  		csvreader = csv.reader(csvfile,delimiter=',')
  		for row in csvreader:
  			dataList.append(row[5:]);

  	X = len(dataList)		#number of lines in .csv
  	Y = len(dataList[1])	#number of entries per line
  
  	height = np.zeros((X,Y))
  	
  	#convert string lists to numpy array
  	for i in range(0,X):
  		for j in range(0,Y):
  			height[i,j] = float(dataList[i][(Y-1)-j])
  	print "end of getElevations"
  	return height

def roadSize(x):
	return{
		'motorway' : 4,
		'trunk' : 4,
		'primary' : 4,
		'secondary' : 3,
		'tertiary' : 3,
		'unclassified' : 0,
		'residential' : 2,
		'service' : 1
		}.get(x,-1)

def isOutOfBounds(filename,lon,lat):
	#global minlat, minlon, maxlat, maxlon
	params = json.load(open(directory + 'ParamFiles/' + filename + '.txt'))
	#global minlat, minlon, maxlat, maxlon
	minlat = params['minlat']
	minlon = params['minlon']
	maxlat = params['maxlat']
	maxlon = params['maxlon']
	
	if lon > maxlon:
		#print "high lon: " + str(lon)
		return True
	if lon < minlon:
		#print "low lon: " + str(lon)
		return True
	if lat > maxlat:
		#print "high lat: " + str(lat)
		return True
	if lat < minlat:
		#print "low lat: " + str(lat)
		return True
	return False

def getLocfromNode(nodes,ref):
	for node in nodes:
		if ref == int(node.get('id')):
			lon = float(node.get('lon'))
			lat = float(node.get('lat'))
			return [lon,lat]

def adjustPoints(filename,lon_old,lat_old,lon1,lat1):
	#finds the point along a line between 2 points that intersects the map boundary
	#global minlat, minlon, maxlat, maxlon
	params = json.load(open(directory + 'ParamFiles/' + filename + '.txt'))
	minlat = params['minlat']
	minlon = params['minlon']
	maxlat = params['maxlat']
	maxlon = params['maxlon']
	#print "fixing lon: " + str(lon_old) + " lat: " + str(lat_old) + " to lon: " + str(lon1) + " lat: " + str(lat1)
	if lon_old < minlon:
		lat_old = lat_old + ((minlon - lon_old)/(lon1 - lon_old))*(lat1 - lat_old)
		lon_old = minlon
	if lon_old > maxlon:
		lat_old = lat_old + ((maxlon - lon_old)/(lon1 - lon_old))*(lat1 - lat_old)
		lon_old = maxlon
	if lat_old < minlat:
		lon_old = lon_old + ((minlat - lat_old)/(lat1 - lat_old))*(lon1 - lon_old)
		lat_old = minlat
	if lat_old > maxlat:
		lon_old = lon_old + ((maxlat - lat_old)/(lat1 - lat_old))*(lon1 - lon_old)
		lat_old = maxlat
	return [lon_old,lat_old]



def polyFile2stlBinItaly(filename):
	global height, i2e, samples
	facetCnt = 0
	params = json.load(open(directory + 'ParamFiles/' + filename + '.txt'))
	minlat = params['minlat']
	minlon = params['minlon']
	maxlat = params['maxlat']
	maxlon = params['maxlon']
	imgSize = params['imgSize']
	#global minlon, minlat, maxlon, maxlat
	import triangle
	import triangle.plot
	#import matplotlib.pyplot as plt
	import elevation2stl
	import makeDome
	#import testipoly
	import numpy as np
	samples = params['samples']
	road_offset = params['road_height']
	print_size = params['print_size']

	maplat_img = interp1d([minlat,maxlat],[imgSize,0])
	maplon_img = interp1d([minlon,maxlon],[0,imgSize])


	if params['elevation']:
		print 'elevation is true?'
		#elevs = testipoly.getElevations(filename)
		elevs = getElevations(filename)
	else:
		print 'flat elevation?'
		elevs = np.zeros((2,2))

	z_scale = float(params['z_scale'])
	if params['zmin']:
		print 'zmin is true'
		height = interp2d(range(0,len(elevs)),range(0,len(elevs[0])),z_scale*(elevs - float(zminVal)));
	else:
		print 'zmin is false'
		print elevs.min(), elevs.max()
		print 'z_scale: ' + str(z_scale)
		print len(elevs), len(elevs[0])
		v = z_scale*(elevs - elevs.min());
		print elevs.shape
		height = interp2d(range(0,len(elevs)),range(0,len(elevs[0])),v);
	print elevs.min(), elevs.max()
	floor = 0
	base_height = params['base_height']
	#i2e = interp1d([0,imgSize],[0,samples-1])
	i2e = interp1d([0,imgSize],[0,samples-1])
	img2printx = interp1d([0,imgSize],[0,print_size])
	img2printy = interp1d([0,imgSize],[print_size,0])

	edges = box(0,0,imgSize,imgSize)
	edge_vs = []
	print 'gonna get the polyfiles'
	#poly2polyFile(filename);
	RoadPSLG = triangle.get_data(params['dir'] + 'PolyFiles/' + filename)
	RoadCDT = triangle.triangulate(RoadPSLG,'p')
	nonRoadPSLG = triangle.get_data(params['dir'] + 'PolyFiles/' + 'i_' + filename)
	nonRoadCDT = triangle.triangulate(nonRoadPSLG,'p')
	print 'gonna write the stl header'
	elevation2stl.startSTLBin(filename)


	for tri in RoadCDT['triangles']:
		#top layer of roads
		pt1 = Point(RoadCDT['vertices'][tri[0]][0],RoadCDT['vertices'][tri[0]][1])
		pt2 = Point(RoadCDT['vertices'][tri[1]][0],RoadCDT['vertices'][tri[1]][1])
		pt3 = Point(RoadCDT['vertices'][tri[2]][0],RoadCDT['vertices'][tri[2]][1])
		#points are scaled to imgSize
		tri_v1 = [img2printx(pt1.x),img2printy(pt1.y),get_height(pt1) + road_offset + base_height]
		tri_v2 = [img2printx(pt2.x),img2printy(pt2.y),get_height(pt2) + road_offset + base_height]
		tri_v3 = [img2printx(pt3.x),img2printy(pt3.y),get_height(pt3) + road_offset + base_height]
		elevation2stl.printTriangleBin(tri_v1,tri_v2,tri_v3,filename)
		facetCnt = facetCnt + 1
	for tri in nonRoadCDT['triangles']:
		#top layer of non-roads
		pt1 = Point(nonRoadCDT['vertices'][tri[0]][0],nonRoadCDT['vertices'][tri[0]][1])
		pt2 = Point(nonRoadCDT['vertices'][tri[1]][0],nonRoadCDT['vertices'][tri[1]][1])
		pt3 = Point(nonRoadCDT['vertices'][tri[2]][0],nonRoadCDT['vertices'][tri[2]][1])
		#points are scaled to imgSize
		tri_v1 = [img2printx(pt1.x),img2printy(pt1.y),get_height(pt1) + base_height]
		tri_v2 = [img2printx(pt2.x),img2printy(pt2.y),get_height(pt2) + base_height]
		tri_v3 = [img2printx(pt3.x),img2printy(pt3.y),get_height(pt3) + base_height]
		elevation2stl.printTriangleBin(tri_v1,tri_v2,tri_v3,filename)
		facetCnt = facetCnt + 1
	for seg in RoadCDT['segments']:
		#build the sides of the roads
		line = LineString(RoadCDT['vertices'][seg])
		if line.touches(edges):
			#along the edges of the print, fill from the base_height to the floor
			#fill road edge sides
			XYpt1 = Point(RoadCDT['vertices'][seg[0]])
			XYpt2 = Point(RoadCDT['vertices'][seg[1]])
			edge_vs.append(XYpt1.coords[0])
			edge_vs.append(XYpt2.coords[0])
			pt1 = [img2printx(XYpt1.x),img2printy(XYpt1.y),get_height(XYpt1) + base_height]
			pt2 = [img2printx(XYpt2.x),img2printy(XYpt2.y),get_height(XYpt2) + base_height]
			pt3 = [img2printx(XYpt1.x),img2printy(XYpt1.y),floor]
			pt4 = [img2printx(XYpt2.x),img2printy(XYpt2.y),floor]
			elevation2stl.printTriangleBin(pt1,pt2,pt3,filename)
			facetCnt = facetCnt + 1
			elevation2stl.printTriangleBin(pt3,pt2,pt4,filename)
			facetCnt = facetCnt + 1
		
		#fill all road sides from the road_offset to the base_height
		#this will cover the gap from base to road height along the edges
		XYpt1 = Point(RoadCDT['vertices'][seg[0]])
		XYpt2 = Point(RoadCDT['vertices'][seg[1]])
		pt1 = [img2printx(XYpt1.x),img2printy(XYpt1.y),get_height(XYpt1) + road_offset + base_height]
		pt2 = [img2printx(XYpt2.x),img2printy(XYpt2.y),get_height(XYpt2) + road_offset + base_height]
		pt3 = [img2printx(XYpt1.x),img2printy(XYpt1.y),get_height(XYpt1) + base_height]
		pt4 = [img2printx(XYpt2.x),img2printy(XYpt2.y),get_height(XYpt2) + base_height]
		elevation2stl.printTriangleBin(pt1,pt2,pt3,filename)
		facetCnt = facetCnt + 1
		elevation2stl.printTriangleBin(pt3,pt2,pt4,filename)
		facetCnt = facetCnt + 1

	for seg in nonRoadCDT['segments']:
		line = LineString(nonRoadCDT['vertices'][seg])
		#fill non-road edge sides
		if line.touches(edges):
			XYpt1 = Point(nonRoadCDT['vertices'][seg[0]])
			XYpt2 = Point(nonRoadCDT['vertices'][seg[1]])
			edge_vs.append(XYpt1.coords[0])
			edge_vs.append(XYpt2.coords[0])
			pt1 = [img2printx(XYpt1.x),img2printy(XYpt1.y),get_height(XYpt1) + base_height]
			pt2 = [img2printx(XYpt2.x),img2printy(XYpt2.y),get_height(XYpt2) + base_height]
			pt3 = [img2printx(XYpt1.x),img2printy(XYpt1.y),floor]
			pt4 = [img2printx(XYpt2.x),img2printy(XYpt2.y),floor]
			elevation2stl.printTriangleBin(pt1,pt2,pt3,filename)
			facetCnt = facetCnt + 1
			elevation2stl.printTriangleBin(pt3,pt2,pt4,filename)
			facetCnt = facetCnt + 1


	HotelPSLG = triangle.get_data(params['dir'] + 'PolyFiles/h' + filename)
	HotelCDT = triangle.triangulate(HotelPSLG,'p')
	TrainPSLG = triangle.get_data(params['dir'] + 'PolyFiles/t' + filename)
	TrainCDT = triangle.triangulate(TrainPSLG,'p')
	FortPSLG = triangle.get_data(params['dir'] + 'PolyFiles/f' + filename)
	FortCDT = triangle.triangulate(FortPSLG,'p')
	
	for tri in HotelCDT['triangles']:
		pt1 = Point(HotelCDT['vertices'][tri[0]][0],HotelCDT['vertices'][tri[0]][1])
		pt2 = Point(HotelCDT['vertices'][tri[1]][0],HotelCDT['vertices'][tri[1]][1])
		pt3 = Point(HotelCDT['vertices'][tri[2]][0],HotelCDT['vertices'][tri[2]][1])
		#points are scaled to imgSize
		tri_v1 = [img2printx(pt1.x),img2printy(pt1.y),get_height(pt1)  + base_height]
		tri_v2 = [img2printx(pt2.x),img2printy(pt2.y),get_height(pt2)  + base_height]
		tri_v3 = [img2printx(pt3.x),img2printy(pt3.y),get_height(pt3)  + base_height]
		elevation2stl.printTriangleBin(tri_v1,tri_v2,tri_v3,filename)
		facetCnt = facetCnt + 1

	for tri in TrainCDT['triangles']:
		pt1 = Point(TrainCDT['vertices'][tri[0]][0],TrainCDT['vertices'][tri[0]][1])
		pt2 = Point(TrainCDT['vertices'][tri[1]][0],TrainCDT['vertices'][tri[1]][1])
		pt3 = Point(TrainCDT['vertices'][tri[2]][0],TrainCDT['vertices'][tri[2]][1])
		#points are scaled to imgSize
		tri_v1 = [img2printx(pt1.x),img2printy(pt1.y),get_height(pt1)  + base_height]
		tri_v2 = [img2printx(pt2.x),img2printy(pt2.y),get_height(pt2)  + base_height]
		tri_v3 = [img2printx(pt3.x),img2printy(pt3.y),get_height(pt3)  + base_height]
		elevation2stl.printTriangleBin(tri_v1,tri_v2,tri_v3,filename)
		facetCnt = facetCnt + 1

	for seg in TrainCDT['segments']:
		#build the sides of the roads
		line = LineString(TrainCDT['vertices'][seg])
		if line.touches(edges):
			#along the edges of the print, fill from the base_height to the floor
			#fill road edge sides
			XYpt1 = Point(TrainCDT['vertices'][seg[0]])
			XYpt2 = Point(TrainCDT['vertices'][seg[1]])
			edge_vs.append(XYpt1.coords[0])
			edge_vs.append(XYpt2.coords[0])
			pt1 = [img2printx(XYpt1.x),img2printy(XYpt1.y),get_height(XYpt1) + base_height]
			pt2 = [img2printx(XYpt2.x),img2printy(XYpt2.y),get_height(XYpt2) + base_height]
			pt3 = [img2printx(XYpt1.x),img2printy(XYpt1.y),floor]
			pt4 = [img2printx(XYpt2.x),img2printy(XYpt2.y),floor]
			elevation2stl.printTriangleBin(pt1,pt2,pt3,filename)
			facetCnt = facetCnt + 1
			elevation2stl.printTriangleBin(pt3,pt2,pt4,filename)
			facetCnt = facetCnt + 1

	for tri in FortCDT['triangles']:
		pt1 = Point(FortCDT['vertices'][tri[0]][0],FortCDT['vertices'][tri[0]][1])
		pt2 = Point(FortCDT['vertices'][tri[1]][0],FortCDT['vertices'][tri[1]][1])
		pt3 = Point(FortCDT['vertices'][tri[2]][0],FortCDT['vertices'][tri[2]][1])
		#points are scaled to imgSize
		tri_v1 = [img2printx(pt1.x),img2printy(pt1.y),get_height(pt1)  + base_height]
		tri_v2 = [img2printx(pt2.x),img2printy(pt2.y),get_height(pt2)  + base_height]
		tri_v3 = [img2printx(pt3.x),img2printy(pt3.y),get_height(pt3)  + base_height]
		elevation2stl.printTriangleBin(tri_v1,tri_v2,tri_v3,filename)
		facetCnt = facetCnt + 1


	markScale = 20000
	hotelLat = 43.780094
	hotelLon = 11.244573
	newlon = ((maplon_img(hotelLon) + 22)/1024.0) * (maxlon - minlon) + minlon
	(tri1,outer1) = makeDome.makeCircle((10.0/1024.0)*(maxlon-minlon),Point(newlon,hotelLat),5)
	print 'dome faces',len(tri1),tri1[0][0].x
	for tri in tri1:
		pt1 = Point(maplon_img(tri[0].x),maplat_img(tri[0].y),tri[0].z)
		pt2 = Point(maplon_img(tri[1].x),maplat_img(tri[1].y),tri[1].z)
		pt3 = Point(maplon_img(tri[2].x),maplat_img(tri[2].y),tri[2].z)
		tri_v1 = [img2printx(pt1.x),img2printy(pt1.y),get_height(pt1)  + base_height + pt1.z*markScale]
		tri_v2 = [img2printx(pt2.x),img2printy(pt2.y),get_height(pt2)  + base_height + pt2.z*markScale]
		tri_v3 = [img2printx(pt3.x),img2printy(pt3.y),get_height(pt3)  + base_height + pt3.z*markScale]
		elevation2stl.printTriangleBin(tri_v1,tri_v2,tri_v3,filename)
		facetCnt = facetCnt + 1

	fortLat = 43.781198
	fortLon = 11.24838
	newlon = ((maplon_img(fortLon) - 16 + 30)/1024.0) * (maxlon - minlon) + minlon
	(tri2,outer2) = makeDome.makeCircle((10.0/1024.0)*(maxlon-minlon),Point(newlon,fortLat + .00008),5)
	for tri in tri2:
		pt1 = Point(maplon_img(tri[0].x),maplat_img(tri[0].y),tri[0].z)
		pt2 = Point(maplon_img(tri[1].x),maplat_img(tri[1].y),tri[1].z)
		pt3 = Point(maplon_img(tri[2].x),maplat_img(tri[2].y),tri[2].z)
		tri_v1 = [img2printx(pt1.x),img2printy(pt1.y),get_height(pt1)  + base_height + pt1.z*markScale]
		tri_v2 = [img2printx(pt2.x),img2printy(pt2.y),get_height(pt2)  + base_height + pt2.z*markScale]
		tri_v3 = [img2printx(pt3.x),img2printy(pt3.y),get_height(pt3)  + base_height + pt3.z*markScale]
		elevation2stl.printTriangleBin(tri_v1,tri_v2,tri_v3,filename)
		facetCnt = facetCnt + 1
	newlon = ((maplon_img(fortLon) + 16 + 30)/1024.0) * (maxlon - minlon) + minlon
	(tri3,outer3) = makeDome.makeCircle((10.0/1024.0)*(maxlon-minlon),Point(newlon,fortLat + .00008),5)
	for tri in tri3:
		pt1 = Point(maplon_img(tri[0].x),maplat_img(tri[0].y),tri[0].z)
		pt2 = Point(maplon_img(tri[1].x),maplat_img(tri[1].y),tri[1].z)
		pt3 = Point(maplon_img(tri[2].x),maplat_img(tri[2].y),tri[2].z)
		tri_v1 = [img2printx(pt1.x),img2printy(pt1.y),get_height(pt1)  + base_height + pt1.z*markScale]
		tri_v2 = [img2printx(pt2.x),img2printy(pt2.y),get_height(pt2)  + base_height + pt2.z*markScale]
		tri_v3 = [img2printx(pt3.x),img2printy(pt3.y),get_height(pt3)  + base_height + pt3.z*markScale]
		elevation2stl.printTriangleBin(tri_v1,tri_v2,tri_v3,filename)
		facetCnt = facetCnt + 1
	

	trainLat = 43.776893
	trainLon = 11.247373
	newlat = ((maplat_img(trainLat) + 20)/1024.0) * (maxlat - minlat) + minlat
	(tri4,outer4) = makeDome.makeCircle((10.0/1024.0)*(maxlon-minlon),Point(trainLon,newlat),5)
	for tri in tri4:
		pt1 = Point(maplon_img(tri[0].x),1024-maplat_img(tri[0].y),tri[0].z)
		pt2 = Point(maplon_img(tri[1].x),1024-maplat_img(tri[1].y),tri[1].z)
		pt3 = Point(maplon_img(tri[2].x),1024-maplat_img(tri[2].y),tri[2].z)
		tri_v1 = [img2printx(pt1.x),img2printy(pt1.y),get_height(pt1)  + base_height + pt1.z*markScale]
		tri_v2 = [img2printx(pt2.x),img2printy(pt2.y),get_height(pt2)  + base_height + pt2.z*markScale]
		tri_v3 = [img2printx(pt3.x),img2printy(pt3.y),get_height(pt3)  + base_height + pt3.z*markScale]
		elevation2stl.printTriangleBin(tri_v1,tri_v2,tri_v3,filename)
		facetCnt = facetCnt + 1
	newlat = ((maplat_img(trainLat) - 20)/1024.0) * (maxlat - minlat) + minlat
	(tri5,outer5) = makeDome.makeCircle((10.0/1024.0)*(maxlon-minlon),Point(trainLon,newlat),5)
	for tri in tri5:
		pt1 = Point(maplon_img(tri[0].x),1024-maplat_img(tri[0].y),tri[0].z)
		pt2 = Point(maplon_img(tri[1].x),1024-maplat_img(tri[1].y),tri[1].z)
		pt3 = Point(maplon_img(tri[2].x),1024-maplat_img(tri[2].y),tri[2].z)
		tri_v1 = [img2printx(pt1.x),img2printy(pt1.y),get_height(pt1)  + base_height + pt1.z*markScale]
		tri_v2 = [img2printx(pt2.x),img2printy(pt2.y),get_height(pt2)  + base_height + pt2.z*markScale]
		tri_v3 = [img2printx(pt3.x),img2printy(pt3.y),get_height(pt3)  + base_height + pt3.z*markScale]
		elevation2stl.printTriangleBin(tri_v1,tri_v2,tri_v3,filename)
		facetCnt = facetCnt + 1





	#build the base
	edge_vs.append([0,0])
	edge_vs.append([0,imgSize])
	edge_vs.append([imgSize,0])
	edge_vs.append([imgSize,imgSize])
	edge_array = np.array(edge_vs)
	base = dict(vertices=edge_array)
	baseT = triangle.triangulate(base)

	for tri in baseT['triangles']:
		pt1 = Point(baseT['vertices'][tri[0]][0],baseT['vertices'][tri[0]][1])
		pt2 = Point(baseT['vertices'][tri[1]][0],baseT['vertices'][tri[1]][1])
		pt3 = Point(baseT['vertices'][tri[2]][0],baseT['vertices'][tri[2]][1])
		tri_v1 = [img2printx(pt1.x),img2printy(pt1.y),floor]
		tri_v2 = [img2printx(pt2.x),img2printy(pt2.y),floor]
		tri_v3 = [img2printx(pt3.x),img2printy(pt3.y),floor]

		elevation2stl.printTriangleBin(tri_v1,tri_v2,tri_v3,filename)
		facetCnt = facetCnt + 1

	# = facetCnt 
	elevation2stl.endSTLBin(filename,facetCnt)
	return edge_vs
