from bs4 import BeautifulSoup
from shapely.geometry import box, Polygon, LineString, Point, MultiLineString, MultiPolygon
from shapely.ops import cascaded_union
import Image, ImageDraw
from scipy.interpolate import interp1d, interp2d
import btpoly_tri
import json
import drawMap

MAP_REGIONS_ON = True;

directory = '/Users/bttaylor/Documents/TactileMap/'
minlon = 0;
minlat = 0;
maxlon = 0;
maxlat = 0;
base_height = 1;
print_size = 80; 	#mm
imgSize = 1000;		#pixel
road_height = 2;
zVal = [];
test_poly = []

def workflow(name,lon,lat,*locs,**attr):
	import downloadMap
	import downloadElevation
	import shutil
	import pprint

	global test_poly

	#set the paramaters to be used in various modules
	(filename,params) = setParams(name,lon,lat,*locs,**attr);
	print 'in workflow' + filename + '.txt'
	
	#filename = 'test'
	#(fime,params) = setParams('test',-80.05,40.4,-80.0,40.45)

	if params['elevation']:
		#Get Elevation data from Google
		downloadElevation.downloadElevation(params['minlon'],params['maxlat'],params['maxlon'],params['minlat'],filename)

	setzVal(filename)

	#Get the .osm Map data from OpenStreetMap
	downloadMap.downloadMap(params['minlon'],params['maxlat'],params['maxlon'],params['minlat'],filename)
	
	#Get the different Map sections as Polygons
	polygons = getRegions(filename)

	#waters = additionalWater(filename)

	#polygons[1]['poly'] = polygons[1]['poly'].union(waters)

	drawMap.drawMap('RoadWater',polygons[0]['poly'],False,(0,0,0))
	drawMap.drawMap('RoadWater',polygons[1]['poly'],False,(0,0,255))
	drawMap.drawMap('RoadWater',polygons[0]['poly'].intersection(polygons[1]['poly']),False)

	drawMap.drawMap('ParkCross',polygons[0]['poly'],False,(0,0,0))
	drawMap.drawMap('ParkCross',polygons[1]['poly'],False,(0,0,255))
	drawMap.drawMap('ParkCross',polygons[2]['poly'],False,(0,255,0))
	drawMap.drawMap('ParkCross',polygons[0]['poly'].intersection(polygons[2]['poly']),False)
	drawMap.drawMap('ParkCross',polygons[1]['poly'].intersection(polygons[2]['poly']),False)

	drawMap.drawMap('OtherCross',polygons[0]['poly'],False,(0,0,0))
	drawMap.drawMap('OtherCross',polygons[1]['poly'],False,(0,0,255))
	drawMap.drawMap('OtherCross',polygons[2]['poly'],False,(0,255,0))
	drawMap.drawMap('OtherCross',polygons[3]['poly'],False,(0,255,255))
	drawMap.drawMap('OtherCross',polygons[0]['poly'].intersection(polygons[3]['poly']),False)
	drawMap.drawMap('OtherCross',polygons[1]['poly'].intersection(polygons[3]['poly']),False)
	drawMap.drawMap('OtherCross',polygons[2]['poly'].intersection(polygons[3]['poly']),False)

	test_poly = polygons[0]['poly']
	drawMap.drawMap('getRegions',polygons[0]['poly'].intersection(polygons[1]['poly']),False)

	print 'Adding Features'
	addFeatures(filename,polygons)
	print 'Done Adding Features'

	drawMap.drawMap('addFeatures',polygons[0]['poly'].intersection(polygons[1]['poly']),False)

	# for poly in polygons:
	# 	#print poly['poly'].bounds
	# 	if poly['poly'].area > 0:
	# 		drawMap.drawMap(filename+'_Allregions',poly['poly'],False)

	# for poly in polygons:
	# 	print poly['poly'].bounds,poly['r_type']
	# 	if poly['poly'].area > 0:
	# 		drawMap.drawMap(filename+'_Allregions2',poly['poly'],False)
	# 	for p2 in polygons:
	# 		if not poly['r_type'] == p2['r_type']:
	# 			x = poly['poly'].intersection(p2['poly'])
	# 			drawMap.drawMap(filename + '_Allregions2',x,False)

	#addBoundaryPoints(polygons)

	#drawMap.drawMap('addBoundary',polygons[0]['poly'].intersection(polygons[1]['poly']),False)

	#This is where markers need to get added
	#This is where features need to get added

	#Create .poly files from the polygons
	poly2polyFile(filename,polygons)

	#for p in polygons:
	#	print p['polyFile']
	#Create the .stl file
	polyFile2stlBin(filename,polygons)

	#Move .stl File?
	#shutil.copy2(directory + 'StlFiles/' + filename + '.stl','./' + params['name'] + '.stl')
	#Create a .poly file for the roadway Polygon
	#poly2polyFile(filename,road_poly,False) #Think this is always false.
	return polygons

def additionalWater(filename):
	import os.path

	osmFile = directory + 'mapFiles/' + filename + '.osm'
	#Make sure the .osm file exists
	if os.path.isfile(osmFile):
		f=open(osmFile);
	else:
		print 'File Read Error',osmFile
		return

	#Use BeautifulSoup to parse the .osm file
	soup = BeautifulSoup(f)

	relations = soup.find_all('relation')
	ways = soup.find_all('way')
	nodes = soup.find_all('node')

	geo_boundary = box(minlon,minlat,maxlon,maxlat)

	geo_polys = []

	for relation in relations:
		tags = relation.find_all('tag')
		for tag in tags:
			if tag.get('k') == 'waterway':
				print 'Waterway Relation'
				members = relation.find_all('member')
				for member in members:
					for way in ways:
						if member['ref'] == way['id']:
							geo_pts = []
							map_pts = []
							print 'Found way listed in relation'
							nds = way.find_all('nd')
							for nd in nds:
								(lon,lat) = getLocfromNode(nodes,int(nd.get('ref')))
								geo_pts.append((lon,lat))

							line = LineString(geo_pts)
							geo_poly = Polygon(line.intersection(geo_boundary))
							#print
							#geo_poly = LineString(geo_pts)

							#geo_poly = geo_poly.intersection(geo_boundary);
							#geo_polys.append(geo_poly)
							drawMap.drawMap('canal',line.intersection(geo_boundary),True)

							#Get the first point in the region
							loc = geo_poly.exterior.coords[0];

							#Collect Map coordinate locations for the point
							x = lon2x(loc[0])
							y = lat2y(loc[1])
							map_pts.append((x,y,base_height + zVal(x,y)))
							#map_pts.append((lon2x(loc[0]),lat2y(loc[1])))

							#iterate through the remaining geographic coordinates of the region
							for coord in geo_poly.exterior.coords[1:]:
								#Collect map coordinate points
								x = lon2x(coord[0])
								y = lat2y(coord[1])
								map_pts.append((x,y,base_height + zVal(x,y)))

								#if MAP_REGIONS_ON:
									#draw.line((lon2img(loc[0]),lat2img(loc[1]),lon2img(coord[0]),lat2img(coord[1])), fill=(0,0,255))						

								#iterate through points for drawing lines	
								loc = coord

							geo_polys.append(Polygon(map_pts))



	return geo_polys

def addBoundaryPoints(regions):
	r_cnt = 0;
	for r in regions:
		print 'adding boundaries to',r['r_type']
		print 'boundary',r['poly'].bounds
		if r['poly'].geom_type == 'MultiPolygon':
			new_polys = []
			change = False
			if r['r_type'] == 1:
				print 'test Road Water Border'
				print r['poly'].intersection(regions[1]['poly'])
			for p in r['poly']:
				#print 'coord count of p',len(p.exterior.coords)
				pflat = Polygon([pt[0:2] for pt in list(p.exterior.coords)])
				#print 'p bounds',p.bounds
				#drawMap.pointCount(p)
				for r2 in regions:
					if not r['r_type'] == r2['r_type']:

						border = pflat.intersection(r2['poly'])
						if r['r_type'] == 1 and r2['r_type'] == 2:
							print 'Road and Water Boundary'
							print border
						#print border
						#print 'adding',r2['r_type']
						#x = pflat.intersection(r2['poly'])
						#pflat = pflat.union(x);
						#drawMap.pointCount(p)
						if 'border_feats' in r2.keys():
							for border_feat in r2['border_feats']:
								L = border_feat['footprint'].intersection(border.buffer(.0000001))
								if not L.is_empty:
									print '***',L.geom_type
									for p in L:
										print list(p.exterior.coords)
								#bf = Polygon([pt[0:2] for pt in list(border_feat.exterior.coords)])
								bf = Polygon([pt[0:2] for pt in list(border_feat['feat'].exterior.coords)])
								#x = p.buffer(.00001).intersection(bf.buffer(.00001))
								#NOTE- Apparently there is some precision loss in either triangle or (more likely)
								#the writing of the .poly file.  Without this buffer setting, points that are distinct
								#to the shapely polygon, appear to be the same when triangulating and stuff breaks
								#x = pflat.intersection(border_feat)
								x = pflat.intersection(bf)
								#x = pflat.intersection(bf.exterior.buffer(0))
								#x = p.buffer(.0000001).intersection(border_feat.buffer(.0000001))
								#Nope. buffer is bad news.
								if not x.is_empty:
									changes = True
									print 'intersection between surface and feature',x
									print 'before union coord cnt:',len(pflat.exterior.coords)
									pflat = pflat.union(x.exterior)
									#print 'after union coord cnt:',len(pflat.exterior.coords)
								#print '**adding',x.geom_type
						#drawMap.pointCount(p)
						new_polys.append(pflat)
			if changes:
				r['poly'] = cascaded_union(new_polys)
			else:
				print "No Features intersected with this Region"
			#print 'after doing stuff'

			if r['poly'].geom_type == 'Polygon':
				r['poly'] = MultiPolygon([r['poly']])
			print 'addBoundary just made a ',r['poly'].geom_type
		else:
			print 'Error','Unexpected region type',r['poly'].geom_type


#creates a json file with relevant parameters
def setParams(name,lon,lat,*locs,**attr):
	import time
	import math

	global minlon, minlat, maxlon, maxlat

	#Add a Timestamp to File Names to avoid overwriting
	filename = name + str(int(math.floor(time.time())))
	print "filename: " + filename

	if len(locs) == 0:
		#if only 1 lat/lon pair is provided, select a default area around the point

		print "Centered at (lon,lat):",lon,lat
		ratio = getLonLatRatio(lat)
		print 'ratio: ' + str(ratio)
		latDelta = dist2LatDegree(.5,lat);
		minlon = lon - latDelta * ratio;
		maxlon = lon + latDelta * ratio;
		minlat = lat - latDelta;
		maxlat = lat + latDelta;
	elif len(locs) == 2:
		#if 2 lat/lon points are provided, set max/min accordingly
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
		#if Elevation attribute is specified
		print "set elevation to: " + str(attr.get('elev'))
		elevation = attr.get('elev')
	else:
		#Adding Elevation Data Defaults Off
		elevation = False
	if 'water' in attr:
		#if Water attribute is specified
		print "set water to: ", attr.get('water')
		water = attr.get('water')
	else:
		#Mapping Water Features Defaults On
		water = True
	if 'park' in attr:
		#if Park attribute is specified
		print "set park to: ", attr.get('park')
		park = attr.get('park')
	else:
		#Mapping Park Features Defaults Off
		park = True
	if 'zmin' in attr:
		#? -- Not sure if this variable is necessary
		zmin = True;
		zminVal = attr.get('zmin')
	else:
		zmin = False;
		zminVal = 0

	print minlon,maxlat,maxlon,minlat
	#road_width is in mm for roadSize = 1
	params = {"name":name, "minlon":minlon, 'minlat':minlat, 'maxlon':maxlon, 'maxlat':maxlat, 
			'elevation':elevation, 'samples':20, 'ratio':1, 'road_width':1, 'road_height':2,
			'imgSize':imgSize, 'base_size':100, 'base_height':1, 'print_size':print_size, 'dir':directory,
			'z_scale':0.2, 'zmin':zmin, 'zminVal':zminVal, 'water':water, 'park':park}
	print params
	json.dump(params, open(directory+ 'ParamFiles/' + filename + '.txt','w'))
	return (filename,params)

#Parse the map into different regions according to specified Features
#Returns a list of Polygons containing the different Feature Regions
#The union of list members should cover the map print region without overlap
def getRegions(filename):
	#Get the parameters for the map
	params = json.load(open(directory + 'ParamFiles/' + filename + '.txt'))

	#Parse the .osm file and return a list of roads
	#roads = getRoads(filename);
	
	#list to hold the regions
	regions = []

	#Get the Road Regions from the Map Data
	roads = getRoads(filename)

	#Turn the road objects into a Polygon
	road_poly = roads2poly(filename,roads)
	#Make Polygons into MultiPolygons for consistancy
	if road_poly.geom_type == "Polygon":
		road_poly = MultiPolygon([road_poly])
	road_obj = {'poly':road_poly,'r_type':1}
	
	regions.append(road_obj)
	#regions.append(road_poly)

	if(params['water']):
		print 'Look for Water Regions'
		#Get the Water Regions from the Map Data
		water_poly = getRegion(filename,waterSelector)
		#Make Polygons into MultiPolygons for consistancy
		add_waters = additionalWater(filename)
		print '*******'
		print add_waters
		drawMap.drawMap('canalTest',add_waters[0],False)
		drawMap.drawMap('canalTest',add_waters[1],False)
		water_poly = water_poly.union(MultiPolygon(add_waters))
		#water_poly = cascaded_union([water_poly,add_waters])
		print add_waters[0].bounds
		if water_poly.geom_type == "Polygon":
			water_poly = MultiPolygon([water_poly])

		water_obj = {'poly':water_poly,'r_type':2}
		#if water_poly.area > 0:
		regions.append(water_obj)
		#regions.append(water_poly)

	if(params['park']):
		print 'Look for Park Regions'
		#Get the Park Regions from the Map Data
		park_poly = getRegion(filename,parkSelector)
		#Make Polygons into MultiPolygons for consistancy
		if park_poly.geom_type == "Polygon":
			park_poly = MultiPolygon([park_poly])
		park_obj = {'poly':park_poly,'r_type':3}
		#if park_poly.area > 0:
		regions.append(park_obj)
		#regions.append(park_poly)

	#Gather all remaining regions
	other_poly = box(0,0,print_size,print_size)
	#other_poly = MultiPolygon([box(0,0,print_size,print_size)])
	other_obj = {'poly':MultiPolygon([other_poly]),'r_type':4}
	regions.append(other_obj)
	#regions.append(box(0,0,print_size,print_size))

	# u_regions = []
	# cnt = 0;
	# for region in regions:
	# 	u_regions.append( Polygon())
	# 	for higher_region in regions[0:cnt]:
	# 		u_regions[cnt] = u_regions[cnt].union(higher_region['poly'])
	# 	cnt += 1
	# 	print 'u_region',cnt,'area',u_regions[cnt-1].area

	# cnt = 0
	# for region in regions:
	# 	regions[cnt]['poly'] = regions[cnt]['poly'].difference(u_regions[cnt])
	# 	cnt += 1

	colors = [(0,0,0),(0,0,255),(0,255,0),(0,255,255)]
	cnt = 0
	for region in regions:
		drawMap.drawMap('rawRegions',region['poly'],False,colors[cnt])
		coord_cnt = []
		for p in region['poly']:
			coord_cnt.append(len(p.exterior.coords))
		cnt += 1
		print 'Region',cnt,'consists of',len(region['poly'].geoms),'polygons'
		print 'Polygon points',list(coord_cnt)

	#Ensure that there are no overlaps in regions
	#The order in which the regions were added defines the heirarchy of region importance
	cnt = 0;
	#iterate through the regions
	for region in regions:
		#iterate through previous regions
		for higher_region in regions[0:cnt]:
			#remove overlap with previous regions
			print 'removing',higher_region['r_type'],'from',region['r_type']
			#print 'before union w/ intersection'
			#for p in higher_region['poly']:
			#	print len(p.exterior.coords)
			higher_region['poly'] = higher_region['poly'].union(higher_region['poly'].intersection(region['poly']))
			if higher_region['poly'].geom_type == 'Polygon':
				higher_region['poly'] = MultiPolygon([higher_region['poly']])
			#region['poly'] = region['poly'].union(region['poly'].intersection(higher_region['poly']))
			#print 'after union w/ intersection'
			#for p in higher_region['poly']:
			#	print len(p.exterior.coords)
			region['poly'] = region['poly'].difference(higher_region['poly']);
			if region['poly'].geom_type == 'Polygon':
				region['poly'] = MultiPolygon([region['poly']])

			pts = []
			for p in higher_region['poly']:
				for pt in p.exterior.coords:
					pts.append((pt[0],pt[1]))
				for hole in p.interiors:
					for pt in hole.coords:
						pts.append((pt[0],pt[1]))
			cross = higher_region['poly'].intersection(region['poly'])
			if not cross.is_empty:
				if cross.geom_type == 'MultiLineString':
					for l in cross:	
						if not (l.coords[0][0],l.coords[0][1]) in pts:
							print 'False 0',l.coords[0]
							#print list(l.coords)
						if not (l.coords[1][0],l.coords[1][1]) in pts:
							print 'False 1',l.coords[1]


		#replace the region in the list with the modified region
		regions[cnt]['poly'] = region['poly'];
		cnt = cnt + 1;


	colors = [(0,0,0),(0,0,255),(0,255,0),(0,255,255)]
	cnt = 0
	for region in regions:
		drawMap.drawMap(filename + '_all',region['poly'],False,colors[cnt])
		coord_cnt = []
		for p in region['poly']:
			coord_cnt.append(len(p.exterior.coords))
		cnt += 1
		print 'Region',cnt,'consists of',len(region['poly'].geoms),'polygons'
		print 'Polygon points',list(coord_cnt)


	return regions

def getRoads(filename):
	import os.path

	#Load parameters for the map
	params = json.load(open(directory + 'ParamFiles/' + filename + '.txt'))

	if MAP_REGIONS_ON:
		#Set up canvas
		map_road = Image.new("RGB",(imgSize,imgSize),"white")
		draw_road = ImageDraw.Draw(map_road)
	
	osmFile = directory + 'mapFiles/' + filename + '.osm'
	#Make sure the .osm file exists
	if os.path.isfile(osmFile):
		f=open(osmFile);
	else:
		print 'File Read Error',osmFile
		return

	#Use BeautifulSoup to parse the .osm file
	soup = BeautifulSoup(f)

	#ways are the features we are interested in
	ways = soup.find_all('way')
	#nodes contain the coordinates of points
	nodes = soup.find_all('node')

	#a list to collect inBound points along a road
	road_pts_g = []
	road_pts_p = []

	#a list to collect complete road objects
	road_objs = []

	#All features (road, water, parks anyway) are specified as ways
	for way in ways:
		#Clear the name value
		name = ''
		#Set the road size to 0
		size = 0
		#Assume the way is not a road
		Road = False

		tags = way.find_all('tag')
		#tags have attributes of the ways
		for tag in tags:

			if str(tag.get('k')) == 'name':
				#store the name of the particular way
				name = tag.get('v')

			if str(tag.get('k')) == 'name:en':
				name = tag.get('v')

			#roads have tags with a 'highway' key
			if str(tag.get('k')) == 'highway':
				#the value of the 'highway' key indicates the road type (e.g. size)
				size = roadSize(str(tag.get('v')))
				if( size > 0):
					#filter out small roads
					Road = True

		#If this way is a valid road, continue processing
		if(Road == True):
			#print name + ' - Road'
			#nd are the nodes (lon,lat) of this particular way
			nds = way.find_all('nd')
			pt1ind = 0
			#we assume the first point (nd) is not located within the bounds of our map
			startInbounds = False
			while startInbounds == False:
				#get the lat/lon of the node listed in the road way
				(lon1,lat1) = getLocfromNode(nodes,int(nds[pt1ind].get('ref')))

				#check if the first point is inBounds
				if(isOutOfBounds(filename,lon1,lat1) == False):
					#first point in bounds, move on
					startInbounds = True
				else:
					#The point is out-of-bounds, get the next point
					(lon2,lat2)= getLocfromNode(nodes,int(nds[pt1ind+1].get('ref')))
					
					#check if the next point is inBounds
					if(isOutOfBounds(filename,lon2,lat2) == False):
						#the 2nd point is inBounds, so move the first point to the boundary border 
						(lon1,lat1) = adjustPoints(filename,lon1,lat1,lon2,lat2)
						#The first 2 points are now inBounds, move on
						startInbounds = True
					else:
						#the first and second point are both out of bounds, ignore point 1 and iterate
						#update the node index to get the next point
						pt1ind = pt1ind + 1
						if(pt1ind == len(nds) - 1):
							print 'Index Error ','No inBound points found'
							return

			#(lon1,lat1) is now inBounds (or on border), pt1ind points to the correspoding node
			#Add the point to the road_pts list
			road_pts_g = [Point(lon1,lat1)]
			road_pts_p = [Point(lon2x(lon1),lat2y(lat1),base_height + zVal(lon2x(lon1),lat2y(lat1)))]

			#The first point is not OOB.  Only continue until the road goes OOB, then stop
			#!!!NOTE - This will mean that a road that goes in and out of bounds will not be fully captured
			OutOfBounds = False

			#iterate through the remaining road nodes (pt1ind + 1 is the next)
			for nd in nds[pt1ind+1:len(nds)]:
				#get the next point
				(lon2,lat2) = getLocfromNode(nodes,int(nd.get('ref')))
				#check if both points are OOB
				if( (isOutOfBounds(filename,lon1,lat1) == True) & (isOutOfBounds(filename,lon2,lat2) == True) ):
					#This should not occurr
					print "ERROR. BOTH OUT OF BOUNDS. Some error?"
					OutOfBounds = True
				#check if the next point is OOB
				elif isOutOfBounds(filename,lon2,lat2) == True:
					#2nd Point is OOB, needs to be plotted to border
					(lon2,lat2) = adjustPoints(filename,lon2,lat2,lon1,lat1)
					#check if the adjusted point is the same as the last point
					if((lon1 == lon2) & (lat1 == lat2)):
						#if the road goes out of bounds, we stop adding points
						#otherwise it will just interpolate it to the edge and we'll have a bunch of identical points
						#might be an issue if the road exits and re-enters the boundaries
						OutOfBounds = True
				#check if only the first point is OOB. Should not happen.
				elif isOutOfBounds(filename,lon1,lat1) == True:
					print "ERROR. Point 1 OOB. error?"
					(lon1,lat1) = adjustPoints(filename,lon1,lat1,lon2,lat2)

				#Check that no OOB condition occurred
				if (OutOfBounds == False):
					#add the next point to the road_pts list
					road_pts_g.append(Point(lon2,lat2))
					x = lon2x(lon2)
					y = lat2y(lat2)
					road_pts_p.append(Point(x,y,base_height + zVal(x,y)))
					#road_pts_p.append(Point(lon2x(lon2),lat2y(lat2),))

					#if Mapping is On, draw the map
					if(MAP_REGIONS_ON):
						draw_road.line((lon2img(lon1),lat2img(lat1),lon2img(lon2),lat2img(lat2)), fill=(0,0,0))
						#draw_road.line((maplon_img(lon1),imgSize - maplat_img(lat1),maplon_img(lon2),imgSize - maplat_img(lat2)), fill=(0,0,0))

				#update the postion and iterate
				lon1 = lon2
				lat1 = lat2

			#create a road object with a name, size and LineString of points
			road_obj = {'name':name,'size':size,'geo_path':LineString(road_pts_g),'print_path':LineString(road_pts_p)}
			#append the road_obj to the road_objs list
			road_objs.append(road_obj)

	#if Mapping is On, save the Map
	if(MAP_REGIONS_ON):
		map_road.save(directory + 'mapImages/' + filename + '_road.png')
	
	print 'found',len(road_objs),'roads'
	#The road_objs contain label information that will probably be needed for labeling purposes.
	return road_objs

def roads2poly(filename,roads=None):
	#load parameters for the map
	params = json.load(open(directory + 'ParamFiles/' + filename + '.txt'))

	#Create Canvas
	map_road = Image.new("RGB",(imgSize,imgSize),"white")
	draw_road = ImageDraw.Draw(map_road)

	grid = box(0,0,print_size,print_size)
	#add a z-Value to grids
	#grid = Polygon([(0,0,base_height),(0,print_size,base_height),(print_size,print_size,base_height),(print_size,0,base_height)])

	#a list for Polygon outlines of roads (there can be disjointed roads)
	road_outlines = []

	#iterate through the road objects passed in
	for road in roads:
		#a road line is widened by (roadSize) * (road_width)
		#road_width is in mm scaled by roadSize (see function for scaling factor)
		#the '0.5' factor is b/c the buffer widens both sides
		fat = road['print_path'].buffer(params['road_width'] * 0.5 * road['size'],16,3,1)
		#buffer() settings - round JointStyle reduced gaps in the road. 
		#CapStyle doesn't seem to match the documentation.  1 is worst. 2,3 seem the same

		#trim the road polygon to fit the print area
		fat = fat.intersection(grid)

		#add the road polygon to the road_outlines list
		road_outlines.append(fat);

	road_poly = cascaded_union(road_outlines)

	if MAP_REGIONS_ON:
		if road_poly.geom_type == 'MultiPolygon':
			color = (255,0,0)
			for single in road_poly:
				bound = single.exterior.coords
				for j in range(1, len(list(bound))):
					draw_road.line((x2img(bound[j-1][0]),y2img(bound[j-1][1]),x2img(bound[j][0]),y2img(bound[j][1])),fill=color,width=2)
				color = (0,0,255)
				drawHoles(filename,single.interiors,draw_road)
		elif road_poly.geom_type == 'Polygon':
			bound = road_poly.exterior.coords
			color = (255,0,0)
			for j in range(1, len(list(bound))):
				draw_road.line((x2img(bound[j-1][0]),y2img(bound[j-1][1]),x2img(bound[j][0]),y2img(bound[j][1])),fill=color,width=2)
			drawHoles(filename,road_poly.interiors,draw_road)

		map_road.save(directory + 'PngFiles/' + filename + '.png')

	return road_poly

def drawHoles(filename,holes,draw):
	#iterate through the holes
	for hole in holes:
		#get a triangle guaranteed to be inside the hole
		p = btpoly_tri.draw_poly(hole.coords[1:]);
		#get an x,y location guaranteed to be inside the hole
		x = p[0][0] + .5*(p[1][0]-p[0][0]) + .25*(p[2][0]-p[0][0])
		y = p[0][1] + .5*(p[1][1]-p[0][1]) + .25*(p[2][1]-p[0][1])
		pp = Point(x,y)
		#draw the interior point
		draw.ellipse((x2img(pp.x)-5,y2img(pp.y)-5,x2img(pp.x)+5,y2img(pp.y)+5), fill=(0,0,255))
		
		#create list of points
		bound = hole.coords
		#iterate starting at 2nd point
		for j in range(1,len(bound)):
			#draw lines
			draw.line((x2img(bound[j-1][0]),y2img(bound[j-1][1]),x2img(bound[j][0]),y2img(bound[j][1])), fill=(0,255,0),width=2)

#getRegion will use the 'featureSelector' function to parse a type of map feature (Water, Park) from .osm data
#'featureSelector' functions must be defined to analyze .osm tags and return 'True/' if the feature matches
#The specifed regions are assumed to be closed loops (i.e. this won't work for roadways or paths)
def getRegion(filename,featureSelector):
	import os.path

	#get map parameters
	params = json.load(open(directory + 'ParamFiles/' + filename + '.txt'))
	
	#set the geographic boundary
	geo_boundary = box(minlon,minlat,maxlon,maxlat)

	if MAP_REGIONS_ON:
		#specify image filepath
		filePath = directory + 'mapImages/' + filename + '_regions.png';
		#check if file exists
		if os.path.isfile(filePath):
			#open if file exists
			mapImg = Image.open(filePath)
		else:
			#create new if file does not exist
			mapImg = Image.new("RGB",(imgSize,imgSize),"white")

		#allow drawing on map
		draw = ImageDraw.Draw(mapImg)

	osmFile = directory + 'mapFiles/' + filename + '.osm'
	#Make sure the .osm file exists
	if os.path.isfile(osmFile):
		f=open(osmFile);
	else:
		print 'File Read Error',osmFile
		return

	#Use BeautifulSoup to parse the .osm file
	soup = BeautifulSoup(f)

	#ways are the features we are interested in
	ways = soup.find_all('way')
	#nodes contain the coordinates of points
	nodes = soup.find_all('node')

	#A list to collect complete region objects
	region_objs = []

	#All features (road, water, parks anyway) are specified as ways
	for way in ways:
		#Clear the name value
		name = ''
		#Set the road size to 0
		size = 0
		#Assume the particular way is not the proper type
		ProperFeature = False

		#Lists to collect region points
		geo_pts = []
		map_pts = []

		tags = way.find_all('tag')
		#tags have attributes of the ways
		for tag in tags:

			if str(tag.get('k')) == 'name':
				#store the name of the particular way
				name = tag.get('v')

			if str(tag.get('k')) == 'name:en':
				print 'name:en Tag Found',tag.get('v')
				name = tag.get('v')

			#Check the tags to see if it matches the specificed Feature Type
			if(featureSelector(tag)):
				#If any tag in the way matches the requirement, the way is the correct type
				ProperFeature = True

		if(ProperFeature == True):
			#nd are the nodes (lon,lat) of this particular way
			nds = way.find_all('nd')
			#iterate through the points of the way
			for nd in nds:
				#get the geographic coordinates
				(lon,lat) = getLocfromNode(nodes,int(nd.get('ref')))
				#store the geographic coordinates
				geo_pts.append((lon,lat))

			#The map feature must be a closed loop (e.g. riverbanks, park areas, ...)
			if geo_pts[0] == geo_pts[-1]:
				
				#Create a Polygon of the geo coordinates.  This may extend beyond the boundary
				geo_poly = Polygon(geo_pts)

				#This may be unnecessary given the closed loop check
				if geo_poly.is_valid:
					print 'Adding map feature:',name

					#Restrict the region to the map boundaries
					geo_poly = geo_poly.intersection(geo_boundary);

					if geo_poly.geom_type == 'MultiPolygon':
						for poly in geo_poly:

							#Get the first point in the region
							loc = poly.exterior.coords[0];

							#Collect Map coordinate locations for the point
							x = lon2x(loc[0])
							y = lat2y(loc[1])
							map_pts.append((x,y,base_height + zVal(x,y)))
							#map_pts.append((lon2x(loc[0]),lat2y(loc[1]),base_height + zVal()))

							#iterate through the remaining geographic coordinates of the region
							for coord in poly.exterior.coords[1:]:
								#Collect map coordinate points
								x = lon2x(coord[0])
								y = lat2y(coord[1])
								map_pts.append((x,y,base_height + zVal(x,y)))
								#map_pts.append((lon2x(coord[0]),lat2y(coord[1])))

								if MAP_REGIONS_ON:
									draw.line((lon2img(loc[0]),lat2img(loc[1]),lon2img(coord[0]),lat2img(coord[1])), fill=(0,0,255))						

								#iterate through points for drawing lines	
								loc = coord

							#Create a Polygon of the region and add it to the Region Object collection
							region_objs.append(Polygon(map_pts))
							map_pts = []

					else:
						#Get the first point in the region
						loc = geo_poly.exterior.coords[0];

						#Collect Map coordinate locations for the point
						x = lon2x(loc[0])
						y = lat2y(loc[1])
						map_pts.append((x,y,base_height + zVal(x,y)))
						#map_pts.append((lon2x(loc[0]),lat2y(loc[1])))

						#iterate through the remaining geographic coordinates of the region
						for coord in geo_poly.exterior.coords[1:]:
							#Collect map coordinate points
							x = lon2x(coord[0])
							y = lat2y(coord[1])
							map_pts.append((x,y,base_height + zVal(x,y)))

							if MAP_REGIONS_ON:
								draw.line((lon2img(loc[0]),lat2img(loc[1]),lon2img(coord[0]),lat2img(coord[1])), fill=(0,0,255))						

							#iterate through points for drawing lines	
							loc = coord

						#Create a Polygon of the region and add it to the Region Object collection
						region_objs.append(Polygon(map_pts))
				else:
					print 'ERROR','invalid polygon from geo coords'
			else:
				print 'Path not polygon'

	if MAP_REGIONS_ON:
		#save the image
		mapImg.save(filePath)

	#Merge the regions into a single Polygon\MultiPolygon object
	region = cascaded_union(region_objs)

	#return the region
	return region

def waterSelector(tag):
	#.osm uses 'waterway'
	if str(tag.get('k')) == 'waterway':
		print '    WaterWay FOUND',tag.get('v')
		#Some 'waterway' tags are paths and will be filtered out by the closed loop requirement in getRegions()
		return True;
		#It seems like the 'riverbanks' tag is unreliable for distinguishing closed loop and path representations
	if str(tag.get('k')) == 'natural':
		#Another marker for water features
		if tag.get('v') == 'water':
			print '     Natural Water',tag.get('v')
			return True;
	return False

def parkSelector(tag):
	#.osm uses 'leisure' key with 'park' value for parks
	if str(tag.get('k')) == 'leisure':
		if tag.get('v') == 'park':
			return True;
	return False

#convert the Feature Polygons to .poly files
def poly2polyFile(filename,regions,inverseRoads=False):
	#set minimum size to prevent a finite math issue
	Hole_Area_Threshold = 10e-20

	#load the map parameters
	params = json.load(open(directory + 'ParamFiles/' + filename + '.txt'))

	#count the current region
	region_cnt = 1;
	for region in regions:
		poly = region['poly']
		if poly.area > 0:
			print 'making poly file number:',region_cnt

			#Count the vertices
			vertices_cnt = 0
			#if poly.geom_type == 'MultiPolygon':
			if poly.geom_type == 'MultiPolygon' or poly.geom_type == 'GeometryCollection':
				for p in poly:
					if p.geom_type == 'Polygon':
						#polygon exteriors have repeated coordinates @ first,last so subtract 1
						vertices_cnt = vertices_cnt + len(p.exterior.coords) - 1
						for hole in p.interiors:
							vertices_cnt = vertices_cnt + len(hole.coords) - 1
			elif poly.geom_type == "Polygon":
				print 'THIS HSOULDL NOT HAPPEN' #delete
			#else:
				vertices_cnt = len(poly.exterior.coords) - 1
				for hole in poly.interiors:
					vertices_cnt = vertices_cnt + len(hole.coords) - 1

			#all the polygons are closed loops so # vertices = # segments
			segment_cnt = vertices_cnt

			#create list for elevation sampled points
			elev_pts = []
			#if elevation is on, additional points should be added
			if params['elevation']:
				#getElevations returns numpy arrays of the sampled lat,lon, elev
				(lats,lons,elev) = getElevations(filename)
				#iterate over the samples 
				for i in range(0,lats.shape[0]):
					for j in range(0,lats.shape[1]):
						p = Point(lon2x(lons[i,j]),lat2y(lats[i,j]))
						#check if a point falls within the region
						if p.intersects(poly):
							elev_pts.append((p.x,p.y))
							vertices_cnt += 1

			#Create .poly file
			polyFile = directory + 'PolyFiles/' + filename + '_' + str(region_cnt)
			region['polyFile'] = polyFile;
			f = open(directory + 'PolyFiles/' + filename + '_' + str(region_cnt) + '.poly','w')

			#Write Header
			f.write("#" + filename + '_' + str(region_cnt) + ".poly\n")
			f.write(str(vertices_cnt) + ' 2 0 1\n')

			v_cnt = 1;
			boundary = 2;

			#Write the vertices
			#if poly.geom_type == 'MultiPolygon':
			if poly.geom_type == 'MultiPolygon' or poly.geom_type == 'GeometryCollection':
				f.write("#the exteriors of the polygon\n")
				ext_cnt = 0;
				#First write the exterior vertices
				for p in poly:

					#?
					if p.geom_type == "Polygon":
						ext_cnt = ext_cnt + 1
						f.write("#exterior " + str(boundary - 1) + "\n")
						(v_cnt,boundary) = writePolygonVertices2Poly(f,p.exterior.coords,v_cnt,boundary)
				f.write("#the end of the exterior polygons\n")

				#Next write the hole vertices
				for p in poly:
					#?
					if p.geom_type == "Polygon":
						for hole in p.interiors:
							f.write("#Hole " + str(boundary - ext_cnt - 1) + "\n")
							(v_cnt,boundary) = writePolygonVertices2Poly(f,hole.coords,v_cnt,boundary)

			else:
				print "I Thought I Killed HTIS"
				print 'type',poly.geom_type
				f.write("#the exterior of the roads polygon\n")
				(v_cnt,boundary) = writePolygonVertices2Poly(f,poly.exterior.coords,v_cnt,boundary)
				f.write("#the first hole\n");
				for hole in poly.interiors:
					f.write("#Hole " + str(boundary - 2) + "\n")
					(v_cnt,boundary) = writePolygonVertices2Poly(f,hole.coords,v_cnt,boundary)

			#if inverseRoads:
			#Add sample points (vs will be empty if elevation is false)
			f.write("#Sample Points\n")
			for i in range(0,len(elev_pts)):
				f.write(str(v_cnt) + ' ' + str(elev_pts[i][0]) + ' ' + str(elev_pts[i][1]) + ' 0\n')
				v_cnt = v_cnt + 1

			#write segment header
			f.write("#the segments\n")
			f.write(str(segment_cnt) + ' 1\n')
			f.write("#exterior segment\n")

			#set values to begin segment portion
			boundary = 2;
			seg_cnt = 1;

			#if poly.geom_type == 'MultiPolygon':
			if poly.geom_type == 'MultiPolygon' or poly.geom_type == 'GeometryCollection':
				first_vert = 1;
				for p in poly:
					if p.geom_type == "Polygon":
						f.write("#exterior seg " + str(boundary - 1) + "\n")
						(seg_cnt,boundary) = writePolygonSegments2Poly(f,p.exterior.coords,seg_cnt,boundary)

				for p in poly:
					if p.geom_type == "Polygon":
						for hole in p.interiors:
							f.write("#interior segment " + str(boundary - ext_cnt - 1) + "\n")
							(seg_cnt,boundary) = writePolygonSegments2Poly(f,hole.coords,seg_cnt,boundary)

			else:
				print 'WTF WHY IS THIS HERE?' #delete
				f.write("#exterior seg 1\n")
				(seg_cnt,boundary) = writePolygonSegments2Poly(f,poly.exterior.coords,seg_cnt,boundary)

				for hole in poly.interiors:
					f.write("#interior segment " + str(boundary - 2) + '\n')
					(seg_cnt,boundary) = writePolygonSegments2Poly(f,hole.coords,seg_cnt,boundary)

			f.write("# X holes\n")

			hole_cnt = 0
			if poly.geom_type == 'MultiPolygon':
				for p in poly:
					hole_cnt = hole_cnt + len(p.interiors)
			elif poly.geom_type == "GeometryCollection":
				#/\?
				hole_cnt = 0
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
								#btpoly_tri triangulates the polygon, returning a single, interior triangle
								pt = btpoly_tri.draw_poly(hole.coords[1:]);
								#find an x,y inside that triangle
								x = pt[0][0] + .5*(pt[1][0]-pt[0][0]) + .25*(pt[2][0]-pt[0][0])
								y = pt[0][1] + .5*(pt[1][1]-pt[0][1]) + .25*(pt[2][1]-pt[0][1])
								#write the interior point
								f.write(str(hole_num) + ' ' + str(x) + ' ' + str(y) +'\n')
								if not (test_poly.contains(Point(x,y))):
									print "ERROR. Point " + str(hole_num) + " is not interior?",pt,x,y
								hole_num = hole_num + 1
							else:
								print "Skipped hole " + str(hole_num) + " is too small: " + str(test_poly.area)

				else:
					for hole in poly.interiors:
						test_poly = Polygon(hole)
						if test_poly.area > Hole_Area_Threshold:
							#btpoly_tri triangulates the polygon, returning a single, interior triangle
							pt = btpoly_tri.draw_poly(hole.coords[1:]);
							#find an x,y inside that triangle
							x = pt[0][0] + .5*(pt[1][0]-pt[0][0]) + .25*(pt[2][0]-pt[0][0])
							y = pt[0][1] + .5*(pt[1][1]-pt[0][1]) + .25*(pt[2][1]-pt[0][1])
							#write the interior point
							f.write(str(hole_num) + ' ' + str(x) + ' ' + str(y) +'\n')
							if not (test_poly.contains(Point(x,y))):
								print "ERROR. Point " + str(hole_num) + " is not interior?",pt,x,y
							hole_num = hole_num + 1
						else:
							print "Skipped hole " + str(hole_num) + " is too small: " + str(test_poly.area)
		region_cnt = region_cnt + 1

def writePolygonVertices2Poly(f,coords,v_cnt,boundary):
	for i in range(1,len(coords)):
		f.write(str(v_cnt) + ' ' + repr(coords[i][0]) + ' ' + repr(coords[i][1]) + ' ' + str(boundary) + '\n')
		v_cnt = v_cnt + 1
	boundary = boundary + 1
	return (v_cnt,boundary)

def writePolygonSegments2Poly(f,coords,seg_cnt,boundary):
	first_vert = seg_cnt
	for i in range(1,len(coords)):
		if(i < len(coords) - 1):
			#segments are incremental 
			f.write(str(seg_cnt) + ' ' + str(seg_cnt) + ' ' + str(seg_cnt + 1) + ' ' + str(boundary) + '\n')
		else:
			#last segment loops back to first
			f.write(str(seg_cnt) + ' ' + str(seg_cnt) + ' ' + str(first_vert) + ' ' + str(boundary) + '\n')
		seg_cnt = seg_cnt + 1
	boundary = boundary + 1
	return (seg_cnt,boundary)

def polyFile2stlBin(filename,regions):
	global height, i2e, samples
	import triangle, triangle.plot
	import elevation2stl
	import numpy as np

	#count the triangles for the .stl file
	facetCnt = 0

	#load map parameters
	params = json.load(open(directory + 'ParamFiles/' + filename + '.txt'))

	samples = params['samples']
	road_offset = params['road_height']
	print_size = params['print_size']

	file_cnt = len(regions)
	
	#the z-value of the bottom of the map
	floor = 0
	#the default z-value of the map surface
	base_height = params['base_height']

	#the edges of the print
	edges = box(0,0,print_size,print_size)

	#points along the base of the triangle that must be included in the triangulation
	edge_vs = [(0,0),(0,print_size),(print_size,0),(print_size,print_size)]

	#write the .stl file header
	elevation2stl.startSTLBin(filename)

	print 'gonna get the polyfiles'
	region_cnt = 0
	for r in regions:
		features = getRegionFeatureParams(r['r_type'])

		print 'region',r['r_type'],r['poly'].geom_type

		r_type = r['r_type']
		if r['poly'].area > 0:

			print 'start stl of',r_type,'facetCnt',facetCnt

			print r['polyFile']
			#load the poly file
			PSLG = triangle.get_data(r['polyFile'])
			#triangulate the surface of the region
			CDT = triangle.triangulate(PSLG,'p')

			#iterate over the region's surface
			for tri in CDT['triangles']:
				#'triangles' refer to the indices of the vertices in the CDT object
				pt1 = Point(CDT['vertices'][tri[0]])
				pt2 = Point(CDT['vertices'][tri[1]])
				pt3 = Point(CDT['vertices'][tri[2]])
				#Poly Files and triangulation only deal w/ the XY plane so Z dimension must be calculated again
				tri_v1 = [pt1.x,pt1.y,zVal(pt1.x,pt1.y) + features['height'] + base_height]
				tri_v2 = [pt2.x,pt2.y,zVal(pt2.x,pt2.y) + features['height'] + base_height]
				tri_v3 = [pt3.x,pt3.y,zVal(pt3.x,pt3.y) + features['height'] + base_height]
				#add the region's surface triangles to the .stl file
				elevation2stl.printTriangleBin(tri_v1,tri_v2,tri_v3,filename)
				facetCnt = facetCnt + 1

			print 'after writing surface of',r_type,'facetCnt',facetCnt
			

			if r['poly'].geom_type == "MultiPolygon":
				poly_cnt = 1
				for poly1 in r['poly']:
					#print 'Working on poly number',poly_cnt,'in region',r_type,poly1.bounds
					poly_cnt += 1
					p = Polygon([pt[0:2] for pt in list(poly1.exterior.coords)])

					#Add edges between adjacent region surfaces
		
					#edges should only be made once per pair
					#check remaining region surfaces for adjacent boundaries
					for r2 in regions[region_cnt + 1:]:
						f2 = getRegionFeatureParams(r2['r_type'])
						if not r['r_type'] == r2['r_type']:
							if r2['poly'].geom_type == "MultiPolygon":
								poly2_cnt = 1
								for p2 in r2['poly']:
									#print 'Finding intersections with poly',poly2_cnt,'in region',r2['r_type']
									poly2_cnt += 1
									p2 = Polygon([pt[0:2] for pt in list(p2.exterior.coords)])
									cross = p.exterior.intersection(p2)
									if cross.geom_type == "LineString":
										cross = MultiLineString([cross])
									if cross.geom_type == 'MultiLineString' or cross.geom_type == 'GeometryCollection':
										for l in cross:
											if l.geom_type == 'LineString':
												
												loc = l.coords[0]

												for coord in l.coords[1:]:
													tri_v1 = [loc[0],loc[1],zVal(loc[0],loc[1]) + base_height + features['height']]
													tri_v2 = [loc[0],loc[1],zVal(loc[0],loc[1]) + base_height + f2['height']]
													tri_v3 = [coord[0],coord[1],zVal(coord[0],coord[1]) + base_height + features['height']]
													tri_v4 = [coord[0],coord[1],zVal(coord[0],coord[1]) + base_height + f2['height']]
													elevation2stl.printTriangleBin(tri_v1,tri_v2,tri_v3,filename)
													facetCnt = facetCnt + 1
													elevation2stl.printTriangleBin(tri_v3,tri_v4,tri_v2,filename)
													facetCnt = facetCnt + 1
													loc = coord
											else:
												print 'ignoring an intersection forming a',l.geom_type
									else:
										print 'region',r_type,'intersects',r2['r_type'],'forming:',cross.geom_type
									#print '*|**',cross
						#print 'after writing edge of',r_type,'and',r2['r_type'],'facetCnt',facetCnt

					#Add edges between surface and adjacent region features

					#check all regions for bordering features
					for r2 in regions:
						if 'border_feats' in r2.keys() and not r2['r_type'] == r_type:
							for bf in r2['border_feats']:
								#flatten the polygon to avoid weird Z-axis issues
								#bf2d = bf
								#bf2d = Polygon([pt[0:2] for pt in list(bf.exterior.coords)])
								cross = p.intersection(bf)
								#cross = p.exterior.intersection(bf.exterior)
								if cross.geom_type == "LineString":
									cross = MultiLineString([cross])
								if cross.geom_type == 'MultiLineString' or cross.geom_type == "GeometryCollection":
									#print 'feature coords',list(bf.exterior.coords)
									#print 'polygon is 3d?',p.has_z
									for l in cross:
										if l.geom_type == 'LineString':
											#print 'Feature intersection line:',list(l.coords)
											loc = l.coords[0]
											#There seems to be an odd error with shapely intersection that can give the wrong z value
											#it seems pretty rare.  Don't really get why.  This work around seems to fix it.
											trueFpoint = bf.intersection(Point(loc[0],loc[1]))
											#print 'correct z:',trueFpoint.z,'returned z:',loc[2]
											for coord in l.coords[1:]:
												tri_v1 = [loc[0],loc[1],zVal(loc[0],loc[1]) + base_height + features['height']]
												#tri_v2 = [loc[0],loc[1],loc[2]]
												tri_v2 = [loc[0],loc[1],trueFpoint.z]

												#There seems to be an odd error with shapely intersection that gives the wrong z values for the line
												trueFpoint = bf.intersection(Point(coord[0],coord[1]))
												#print 'correct z:',trueFpoint.z,'returned z:',coord[2]

												tri_v3 = [coord[0],coord[1],zVal(coord[0],coord[1]) + base_height + features['height']]
												#tri_v4 = [coord[0],coord[1],coord[2]]
												tri_v4 = [coord[0],coord[1],trueFpoint.z]
												elevation2stl.printTriangleBin(tri_v1,tri_v2,tri_v3,filename)
												facetCnt = facetCnt + 1
												elevation2stl.printTriangleBin(tri_v3,tri_v4,tri_v2,filename)
												facetCnt = facetCnt + 1
												loc = coord
										else:
											print 'ignoring an intersection forming a',l.geom_type
								else:
									print 'region',r_type,'intersects features in',r2['r_type'],'forming:',cross.geom_type

					for hole in poly1.interiors:
						h = Polygon([pt[0:2] for pt in list(hole.coords)])
						#edges should only be made once per pair
						#check remaining region surfaces for adjacent boundaries
						for r2 in regions[region_cnt + 1:]:
							f2 = getRegionFeatureParams(r2['r_type'])
							if not r['r_type'] == r2['r_type']:
								if r2['poly'].geom_type == "MultiPolygon":
									poly2_cnt = 1
									for p2 in r2['poly']:
										#print 'Finding intersections with poly',poly2_cnt,'in region',r2['r_type']
										poly2_cnt += 1
										p2 = Polygon([pt[0:2] for pt in list(p2.exterior.coords)])
										cross = h.exterior.intersection(p2)
										if cross.geom_type == "LineString":
											cross = MultiLineString([cross])
										if cross.geom_type == 'MultiLineString' or cross.geom_type == 'GeometryCollection':
											for l in cross:
												if l.geom_type == 'LineString':
													#print 'intersection line:',list(l.coords)
													loc = l.coords[0]
													for coord in l.coords[1:]:
														tri_v1 = [loc[0],loc[1],zVal(loc[0],loc[1]) + base_height + features['height']]
														tri_v2 = [loc[0],loc[1],zVal(loc[0],loc[1]) + base_height + f2['height']]
														tri_v3 = [coord[0],coord[1],zVal(coord[0],coord[1]) + base_height + features['height']]
														tri_v4 = [coord[0],coord[1],zVal(coord[0],coord[1]) + base_height + f2['height']]
														elevation2stl.printTriangleBin(tri_v1,tri_v2,tri_v3,filename)
														facetCnt = facetCnt + 1
														elevation2stl.printTriangleBin(tri_v3,tri_v4,tri_v2,filename)
														facetCnt = facetCnt + 1
														loc = coord
												else:
													print 'ignoring an intersection forming a',l.geom_type
										else:
											print 'region',r_type,'intersects',r2['r_type'],'forming:',cross.geom_type
										#print '*|**',cross
							#print 'after writing edge of',r_type,'and',r2['r_type'],'facetCnt',facetCnt

						#Add edges between sruface and adjacent region features

						#check all regions for bordering features
						for r2 in regions:
							if 'border_feats' in r2.keys() and not r2['r_type'] == r_type:
								for bf in r2['border_feats']:
									#flattn the polygon to avoid weird Z-axis issues
									#bf2d = bf
									#bf2d = Polygon([pt[0:2] for pt in list(bf.exterior.coords)])
									cross = h.intersection(bf)
									#cross = h.exterior.intersection(bf.exterior)
									if cross.geom_type == "LineString":
										cross = MultiLineString([cross])
									if cross.geom_type == 'MultiLineString' or cross.geom_type == "GeometryCollection":
										#print 'feature coords',list(bf.exterior.coords)
										#print 'polygon is 3d?',p.has_z
										for l in cross:
											if l.geom_type == 'LineString':
												#print 'Feature intersection line:',list(l.coords)
												loc = l.coords[0]
												#There seems to be an odd error with shapely intersection that can give the wrong z value
												#it seems pretty rare.  Don't really get why.  This work around seems to fix it.
												trueFpoint = bf.intersection(Point(loc[0],loc[1]))
												#print 'correct z:',trueFpoint.z,'returned z:',loc[2]
												for coord in l.coords[1:]:
													tri_v1 = [loc[0],loc[1],zVal(loc[0],loc[1]) + base_height + features['height']]
													#tri_v2 = [loc[0],loc[1],loc[2]]
													tri_v2 = [loc[0],loc[1],trueFpoint.z]

													#There seems to be an odd error with shapely intersection that gives the wrong z values for the line
													trueFpoint = bf.intersection(Point(coord[0],coord[1]))
													#print 'correct z:',trueFpoint.z,'returned z:',coord[2]

													tri_v3 = [coord[0],coord[1],zVal(coord[0],coord[1]) + base_height + features['height']]
													#tri_v4 = [coord[0],coord[1],coord[2]]
													tri_v4 = [coord[0],coord[1],trueFpoint.z]
													elevation2stl.printTriangleBin(tri_v1,tri_v2,tri_v3,filename)
													facetCnt = facetCnt + 1
													elevation2stl.printTriangleBin(tri_v3,tri_v4,tri_v2,filename)
													facetCnt = facetCnt + 1
													loc = coord
											else:
												print 'ignoring an intersection forming a',l.geom_type
									else:
										print 'region',r_type,'intersects features in',r2['r_type'],'forming:',cross.geom_type



					#check for intersection with the print edges
					cross = p.exterior.intersection(edges.exterior)
					if cross.geom_type == "LineString":
						cross = MultiLineString([cross])
					if cross.geom_type == 'MultiLineString':
						for l in cross:
							#print 'intersection line:',list(l.coords)
							loc = l.coords[0]
							edge_vs.append((loc[0],loc[1]))
							for coord in l.coords[1:]:
								tri_v1 = [loc[0],loc[1],zVal(loc[0],loc[1]) + base_height + features['height']]
								tri_v2 = [loc[0],loc[1],floor]
								tri_v3 = [coord[0],coord[1],zVal(coord[0],coord[1]) + base_height + features['height']]
								tri_v4 = [coord[0],coord[1],floor]
								elevation2stl.printTriangleBin(tri_v1,tri_v2,tri_v3,filename)
								facetCnt = facetCnt + 1
								elevation2stl.printTriangleBin(tri_v3,tri_v4,tri_v2,filename)
								facetCnt = facetCnt + 1
								loc = coord
								edge_vs.append((coord[0],coord[1]))
					else:
						'Irregular intersection between',r_type,'and the edge',cross.geom_type

			else:
				print 'Non-MultiPolygon Region?',r['r_type']
		

			#check for intersection of features with print edges
			if 'border_feats' in r.keys():
				for bf in r['border_feats']:
					cross = bf.intersection(edges)
					#cross = bf.exterior.intersection(edges.exterior)
					if cross.geom_type == "LineString":
						cross = MultiLineString([cross])
					if cross.geom_type == 'MultiLineString':
						for l in cross:
							#print 'Edge/Feature intersection line:',list(l.coords)
							loc = l.coords[0]
							trueFpoint = bf.intersection(Point(loc[0],loc[1]))
							edge_vs.append((loc[0],loc[1]))
							for coord in l.coords[1:]:
								tri_v1 = [loc[0],loc[1],trueFpoint.z]
								tri_v2 = [loc[0],loc[1],floor]
								trueFpoint = bf.intersection(Point(coord[0],coord[1]))
								tri_v3 = [coord[0],coord[1],trueFpoint.z]
								tri_v4 = [coord[0],coord[1],floor]
								elevation2stl.printTriangleBin(tri_v1,tri_v2,tri_v3,filename)
								facetCnt = facetCnt + 1
								elevation2stl.printTriangleBin(tri_v3,tri_v4,tri_v2,filename)
								facetCnt = facetCnt + 1
								loc = coord
								edge_vs.append((coord[0],coord[1]))
					else:
						'Irregular intersection between',r_type,'and the edge',cross.geom_type
					#print '*|**',cross

					#check if the border feature intersects any other border features
					#NOTE- THIS IS NOT QUITE CORRECT
					#the feature triangles have already been calculated, so if there is a misalignment between points
					#a manifold or gap issue will arise.  This code fills the gap nicely, but does not match points properly
					for r2 in regions[region_cnt + 1:]:
						if 'border_feats' in r2.keys():
							for bf2 in r2['border_feats']:
								cross = bf.exterior.intersection(bf2.exterior)
								if cross.geom_type == "LineString":
									cross = MultiLineString([cross])
								if cross.geom_type == 'MultiLineString':
									print 'First Feature',list(bf.exterior.coords)
									print 'Second Feature',list(bf2.exterior.coords)
									for l in cross:
										print 'Feature/Feature intersection line:',list(l.coords)
										loc = l.coords[0]
										trueFpoint1 = bf.intersection(Point(loc[0],loc[1]))
										trueFpoint2 = bf2.intersection(Point(loc[0],loc[1]))
										edge_vs.append((loc[0],loc[1]))
										for coord in l.coords[1:]:
											tri_v1 = [loc[0],loc[1],trueFpoint1.z]
											tri_v2 = [loc[0],loc[1],trueFpoint2.z]
											trueFpoint1 = bf.intersection(Point(coord[0],coord[1]))
											trueFpoint2 = bf2.intersection(Point(coord[0],coord[1]))
											tri_v3 = [coord[0],coord[1],trueFpoint1.z]
											tri_v4 = [coord[0],coord[1],trueFpoint2.z]
											elevation2stl.printTriangleBin(tri_v1,tri_v2,tri_v3,filename)
											facetCnt = facetCnt + 1
											elevation2stl.printTriangleBin(tri_v3,tri_v4,tri_v2,filename)
											facetCnt = facetCnt + 1
											loc = coord
											edge_vs.append((coord[0],coord[1]))
								else:
									'Irregular intersection between',r_type,'and the edge',cross.geom_type
								#print '*|**',cross
									print 'after writing feature edges of',r_type,'and',r2['r_type'],'facetCnt',facetCnt


		if 'add_tris' in r.keys():
			for tri in r['add_tris']:
				#print tri
				pt1 = [tri[0].x, tri[0].y, tri[0].z]
				pt2 = [tri[1].x, tri[1].y, tri[1].z]
				pt3 = [tri[2].x, tri[2].y, tri[2].z]
					#pt1 = [tri[0].x, tri[0].y, tri[0].z + base_height]
					#pt2 = [tri[1].x, tri[1].y, tri[1].z + base_height]
					#pt3 = [tri[2].x, tri[2].y, tri[2].z + base_height]
				elevation2stl.printTriangleBin(pt1,pt2,pt3,filename)
				facetCnt += 1

		print 'after writing features of',r_type,'facetCnt',facetCnt

		region_cnt += 1

	#build the base
	#turn the edge_vs into an array for converting to triangle object for triangulation
	edge_array = np.array(edge_vs)
	base = dict(vertices=edge_array)
	#triangulate the base
	baseT = triangle.triangulate(base)

	#iterate over the base triangles
	for tri in baseT['triangles']:
		pt1 = Point(baseT['vertices'][tri[0]][0],baseT['vertices'][tri[0]][1])
		pt2 = Point(baseT['vertices'][tri[1]][0],baseT['vertices'][tri[1]][1])
		pt3 = Point(baseT['vertices'][tri[2]][0],baseT['vertices'][tri[2]][1])
		tri_v1 = [pt1.x,pt1.y,floor]
		tri_v2 = [pt2.x,pt2.y,floor]
		tri_v3 = [pt3.x,pt3.y,floor]

		elevation2stl.printTriangleBin(tri_v1,tri_v2,tri_v3,filename)
		facetCnt = facetCnt + 1

	#insert the correct number of triangles and close the .stl file
	elevation2stl.endSTLBin(filename,facetCnt)
	#?
	return edge_vs

#def get_height(pt):
#	global height, i2e, samples
#	return height((samples-1) - i2e(pt.x),i2e(pt.y))[0]

# def getRegionOffset(r_type):
# 	if r_type == 1:
# 		#Road
# 		return road_height;
# 	elif r_type == 2:
# 		#Water
# 		return 0;
# 	elif r_type == 3:
# 		#Park
# 		return 0;
# 	elif r_type ==4:
# 		#filler
# 		return 0;
# 	return 0

def getRegionFeatureParams(r_type):
	if r_type == 1:
		#Roads
		features = {'height':road_height,'type':'Flat'}
		return features;
	elif r_type == 2:
		#Water
		features = {'height':0, 'type':'Flat'}
		#features = {'height':.5,'type':'Dome','x_space':27,'y_space':25,'size':2}
		return features;
	elif r_type == 3:
		#Park
		#features = {'height':0,'type':'Flat'}
		features = {'height':0,'type':'Diamond','x_space':25,'y_space':25,'size':2}
		return features;
	else :
		#Other
		features = {'height':0.5,'type':'Flat'}
		return features;


#return the geographic coordinate of the node specified by ref
def getLocfromNode(nodes,ref):
	for node in nodes:
		if ref == int(node.get('id')):
			lon = float(node.get('lon'))
			lat = float(node.get('lat'))
			return (lon,lat)

#will return the intersection point between a line from (lon1,lat1) to (lon_old,lat_old) and the boundary
def adjustPoints(filename,lon_old,lat_old,lon1,lat1):
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
	return (lon_old,lat_old)

#check if lon,lat are within the boundary
def isOutOfBounds(filename,lon,lat):	
	if lon > maxlon:
		return True
	if lon < minlon:
		return True
	if lat > maxlat:
		return True
	if lat < minlat:
		return True
	return False

#return scale factor for different road types
def roadSize(x):
	return{
		'motorway' : 3,
		'trunk' : 3,
		'primary' : 3,
		'secondary' : 2.5,
		'tertiary' : 2.5,
		'unclassified' : 0,
		'residential' : 2, #0, #2
 		'service' : 1 #0  #1
		}.get(x,-1)

#convert geographic coordinate to print coordinate
def lon2x(lon):
	scale = interp1d([minlon,maxlon],[0,print_size])
	return scale(lon)

#convert geographic coordinate to print coordinate
def lat2y(lat):
	scale = interp1d([minlat,maxlat],[0,print_size])
	return scale(lat)

#convert geographic coordinate to image coordinate
def lon2img(lon):
	return x2img(lon2x(lon))

#convert geographic coordinate to image coordinate
def lat2img(lat):
	return -1*y2img(lat2y(lat))

#convert print coordinate to image coordinate
def x2img(x):
	return x*imgSize/print_size

#convert print coordinate to image coordinate
def y2img(y):
	return -1*y*imgSize/print_size

#convert print coordinate to image coordinate
def x2lon(x):
	scale = interp1d([0,print_size],[minlon,maxlon])
	return scale(x)

#convert print coordinate to image coordinate
def y2lat(y):
	scale = interp1d([0,print_size],[minlat,maxlat])
	return scale(y)

#This should return a default latitude to map a consistent area size when only given a point location
#NOTE - this works roughly for Pittsburgh's Latitude
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

def getElevations(filename):
	import csv
	import numpy as np

	#create lists to hold the elevation data from the .csv file
	dataList = []
	latList = []
	lonList = []

	#read the .csv file
	with open(directory + 'elevationFiles/' + filename + '.csv','rb') as csvfile:
		csvreader = csv.reader(csvfile,delimiter=',')
		for row in csvreader:
			#ignore the line number (entry 1)
			#store lat of start and end of each row (entry 2-5)
			latList.append([float(row[1]),float(row[3])])
			lonList.append([float(row[2]),float(row[4])])
			#the remaining points are sampled elevations (entry 6+)
			dataList.append(row[5:]);

	X = len(dataList)		#number of lines in .csv
	Y = len(dataList[1])	#number of entries per line
  
	#create arrays for lat, lon, height data
	lats = np.zeros((X,Y))
	lons = np.zeros((X,Y))
	height = np.zeros((X,Y))

	for i in range(0,X):
		for j in range(0,Y):
			#elevation csv only lists (lat,lon) at the endpoints of a row
			#(lat,lon) must be interpolated
			lats[i,j] = latList[i][0] + (float(i)*(latList[i][0] - latList[i][1]))/(Y - 1)
			lons[i,j] = lonList[j][0] - (float(j)*(lonList[j][0] - lonList[j][1]))/(X - 1)
			#height values are per row
			height[i,j] = float(dataList[i][j])
			#NOTE - it's possible that the mapping here will not be minlat/lon -> maxlat/lon

	return lats,lons,height


def setzVal(filename):
	import numpy as np
	global zVal
	#load map parameters
	params = json.load(open(directory + 'ParamFiles/' + filename + '.txt'))

	#get elevation points for interpolating the height
	if params['elevation']:
		print 'Adding elevation data'
		(lats,lons,elevs) = getElevations(filename)

		#try print coords instead of geo coords
		lats = print_size * (lats - lats.min()) / (lats.max() - lats.min())
		lons = print_size * (lons - lons.min()) / (lons.max() - lons.min())
	else:
		print 'No elevation data'
		#elevation is set to 0
		elevs = np.zeros((2,2))
		#lat and lon arrays are created to match interpolation
		lons = np.array([[minlon,maxlon]])
		lats = np.array([[minlat],[maxlat]])

		print 'print_size in setzVal',print_size
		#try print coords instead of geo coords
		lons = np.array([[0,print_size]])
		lats = np.array([[0],[print_size]])

	z_scale = float(params['z_scale'])

	if params['zmin']:
		#zmin will set a minimum elevation to treat as the low point rather than using the lowest elevation in the data
		#This will ensure that tiled maps use a consistent z-axis reference
		print 'Using common Z reference',params['zminVal']
		zVal = interp2d(lons[0,:],lats[:,0],z_scale*(elevs - float(params['zminVal'])))
	else:
		print 'Scaling height relative to local minimum:',elevs.min(),'max:',elevs.max()
		v = z_scale*(elevs - elevs.min());
		#interpolate the elevation over the latitude/longitude range
		zVal = interp2d(lons[0,:],lats[:,0],v)
		#height = interp2d(range(0,len(elevs)),range(0,len(elevs[0])),v);
	

def addFeatures(filename,regions):
	from shapely.affinity import translate
	#print 'getting zVal',zVal(-80.025,40.425)

	for r in regions:
		if not r['poly']:
			print 'no polygon',r['r_type']
		else:
			#get the feature specifications for the region
			print 'r_type',r['r_type']
			features = getRegionFeatureParams(r['r_type'])
			print 'features',features
			#adjust the elevation by the feature height
			r['poly'] = translate(r['poly'],0,0,features['height'])

			if not features['type'] == 'Flat':
				r['border_feats'] = []
				r['add_tris'] = []
				#Iterate over the separate polygons in the region
				poly_cnt = 0
				if r['poly'].geom_type == 'MultiPolygon':
					new_polys = []
					for p in r['poly']:
						#run the tesselation
						(neg, tris, border_features) = tesselateArea(p,features,regions)
						#(neg, tris, border_features) = tesselateArea(p,features)
						r['border_feats'] += border_features;
						r['add_tris'] += tris
						#Can't access the Polygons within the MultiPolygon via indexing, so store new ones
						new_polys.append(neg)
					#Set the region polygon to the union of the new polygons
					r['poly'] = cascaded_union(new_polys)
				elif r['poly'].geom_type == 'Polygon':
					print 'THIS SHOULD NOT HAPPEN'  #delete
					#run the tesselation
					(neg, tris, border_features) = tesselateArea(r['poly'],features,regions)
					r['border_feats'] += border_features;
					r['add_tris'] += tris
				
					#Set the region polygon to the union of the new polygons
					r['poly'] = neg

def tesselateArea(poly,features,regions):
#def tesselateArea(poly,features):
	#poly = r['poly']
	global test_poly
	#r_type = r['r_type']

	import math
	import numpy as np

	#the specified feature type
	f_type = features['type']
	#the specified features size
	f_size = features['size']
	#the region height adjustment
	f_height = features['height']

	#feature spacing in mm
	tess_width = features['x_space'];
	tess_height = features['y_space'];

	#region boundary
	bounds = poly.bounds;
	#print 'tesselate bounds',bounds

	#Create the polygon to triangulate (the negative space around extruded features)
	neg_space = poly;

	#Create a list to hold all the triangles that compose the features
	feat_tris = []

	#calculate the number of tesselation in x,y
	tess_x_cnt = int(math.ceil(abs(bounds[2] - bounds[0] ) / tess_width));
	tess_y_cnt = int(math.ceil(abs(bounds[3] - bounds[1] ) / tess_height));

	#create an array to track whether the tesselation is wholely, partially or not inside the region
	tess_status = np.zeros((tess_x_cnt,tess_y_cnt))

	#create a list to contain polygons of the features on the border
	border_feats = []

	x0 = bounds[0];
	y0 = bounds[1];
	for i in range(0,tess_x_cnt):
		for j in range(0,tess_y_cnt):
			tess = box(x0 + tess_width*i,y0 + tess_height*j,x0 + tess_width*(i+1), y0 + tess_height*(j+1))
			center = Point((tess.bounds[0] + tess.bounds[2])/2,(tess.bounds[1] + tess.bounds[3])/2)
			border_feat = Polygon() #empty
			if poly.contains(tess):
				#print 'i,j',i,j,'Tess bounds',tess.bounds
				#print 'feature center',center
				tess_status[i,j] = 1;
				#get the feature for the tesselation
				#(tris,outer) = getFeature(center,f_type);
				z_adj = base_height + f_height + zVal(center.x,center.y)
				(tris,outer) = getFeature(center,f_type,f_size,z_adj);
				feat_tris += tris
				#remove the feature region from the negative space
				pts = []		#Maybe should not use Point() in makeDome? just a list of coords?
				for o in outer:
					pts.append((o.x,o.y))
				neg_space = neg_space.difference(Polygon(pts))
			elif poly.intersects(tess):
				tess_status[i,j] = 2;
				#get the feature area for the tesselation
				z_adj = base_height + f_height + zVal(center.x,center.y)
				(tris,outer) = getFeature(center,f_type,f_size,z_adj);
				#feat_tris += tris
				#clip the feature area to the portion within the polygon region
				pts = []		#Maybe should not use Point() in makeDome? just a list of coords?
				for o in outer:
					pts.append((o.x,o.y))

				feat_area = Polygon(pts)
				#the feature area might be smaller than the tesselation and wholely contained
				if poly.contains(feat_area):
					neg_space = neg_space.difference(feat_area)
					feat_tris += tris
				#the feature area may not overlap the polygon (even though its tesselation grid does)
				elif feat_area.intersects(poly):

					for r in regions:
						if feat_area.intersects(r['poly']):
							if not poly.within(r['poly']):
								print '***Feature intersects', r['r_type']
								temp = r['poly'].union(r['poly'].intersection(feat_area))
								#higher_region['poly'] = higher_region['poly'].union(higher_region['poly'].intersection(region['poly']))
								coord_cnt = []
								for p in temp:
									#print p
									coord_cnt.append(len(p.exterior.coords))
								print 'Polygon points',list(coord_cnt)
								r['poly'] = temp


					print 'border Feature',feat_area.bounds
					print feat_area
					#print 'BEF',neg_space.intersection(test_poly)
					neg_space = neg_space.difference(feat_area)
					test_poly = test_poly.difference(feat_area)
					#print 'AFT',neg_space.intersection(test_poly)
					drawMap.drawMap('testPoly',neg_space.intersection(test_poly),False)
					for tri in tris:
						newTris = trimTri(tri,poly)
						feat_tris += newTris

						for nT in newTris:
							#triPoly = Polygon([(nT[0][0],nT[0][1],nT[0][2]),(nT[1][0],nT[1][1],nT[1][2]),(nT[2][0],nT[2][1],nT[2][2])])
							triPoly = Polygon([(nT[0].x,nT[0].y,nT[0].z),(nT[1].x,nT[1].y,nT[1].z),(nT[2].x,nT[2].y,nT[2].z)])
							#print triPoly
							border_feat = border_feat.union(triPoly)
						#print 'trimmed tri:',t
						#for newTri in newTris:
						#	feat_tris += t
			if not border_feat.is_empty:
				#border_feats.append({'feat':border_feat,'footprint':feat_area})
				border_feats.append(border_feat)
			
	return neg_space, feat_tris, border_feats


def getFeature(location,type,size,z_adj):
	import makeDome
	if(type == 'Dome'):
		(tri,outer) = makeDome.makeCircle(size,location,2,z_adj);
		return (tri,outer);
	elif(type == 'Diamond'):
		(tri,outer) = makeDome.makeCircle(size,location,1,z_adj);
		return (tri,outer);
	else:
		print 'feature type error'
		return [];

def trimTri(tri,poly):
	poly = Polygon([pt[0:2] for pt in list(poly.exterior.coords)])
	nTri = Polygon([(tri[0].x,tri[0].y,tri[0].z),(tri[1].x,tri[1].y,tri[1].z),(tri[2].x,tri[2].y,tri[2].z)])
	nTri = nTri.intersection(poly)

	if not nTri.is_empty:
		if len(nTri.exterior.coords) == 3:
			print 'first:',nTri.exterior.coords[0]
			print 'Trimming tRiangle has 3 points!'
			print list(nTri.exterior.coords)
			return [[Point(nTri.exterior.coords[0]),Point(nTri.exterior.coords[1]),Point(nTri.exterior.coords[2])]]
			#return [list(nTri.exterior.coords)]
		elif len(nTri.exterior.coords) == 4:
			#print 'Trimming tRiangle has 4 points!'
			return [[Point(nTri.exterior.coords[0]),Point(nTri.exterior.coords[1]),Point(nTri.exterior.coords[2])]]
			#return [list(nTri.exterior.coords[0:3])]
		elif len(nTri.exterior.coords) == 5:
			#print 'Trimming tRiangle has 5 points!'
			return [[Point(nTri.exterior.coords[0]),Point(nTri.exterior.coords[1]),Point(nTri.exterior.coords[2])],
					[Point(nTri.exterior.coords[2]),Point(nTri.exterior.coords[3]),Point(nTri.exterior.coords[4])]]
			#return [nTri.exterior.coords[0:3]] + [nTri.exterior.coords[2:5]]
		else:
			print 'Trimming tirangle has',len(nTri.exterior.coords),'points!!'
			return []
	else:
		return []


#This may be a problem.  It seems like manually trimming the individual triangles (rather than use some intersection formula)
#may lead to slight discrepancies such that the trimmed triangles to not 'intersect' with the lines forming the polygon when
#tested later.
def trimTri2(tri,poly):
	#Having a z-dimension does weird things w/ intersection in shapely
	#The operations are done in the XY plane, but the Z value is some averaging(?) of the present values
	#to preserve the Z of the triangle being trimmed, remove the z dimension of the polygon
	poly = Polygon([pt[0:2] for pt in list(poly.exterior.coords)])
	#print 'tri in trimTri',tri
	#tris is a set of triagnles in 3d space
	good_pts = [];
	bad_pts = [];
	if tri[0].intersects(poly):
		good_pts.append(tri[0])
	else:
		bad_pts.append(tri[0])
	if tri[1].intersects(poly):
		good_pts.append(tri[1])
	else:
		bad_pts.append(tri[1])
	if tri[2].intersects(poly):
		good_pts.append(tri[2]);
	else:
		bad_pts.append(tri[2]);

	if(len(good_pts) == 3):
		#The triangle is wholely inbounds
		return [tri];
	if(len(good_pts) == 0):
		#The triangle is wholely outOfBounds
		return []
	if(len(good_pts) == 1):
		#the triangle has one good point
		pt1 = good_pts[0]
		#find the point along the first leg that intersects with the Polygon
		l1 = LineString([good_pts[0],bad_pts[0]])
		#print 'in triTrim line:',list(l1.coords)
		#print 'poly z',
		l = l1.intersection(poly)
		#print 'l1.intersection(poly):',list(l.coords)
		if l.geom_type == 'MultiLineString':
			print 'Irregular polygon path splitting feature triangle'
			#the line will start at good_pts[0] (==coord[0]), coord[1] will be next intersection
			pt2 = Point(l[0].coords[1])
		else:
			#the line will start at good_pts[0] (==coord[0]), coord[1] will be next intersection
			pt2 = Point(l.coords[1])
		#find the point along the second leg that intersects with the Polygon
		l2 = LineString([good_pts[0],bad_pts[1]])
		l = l2.intersection(poly)
		if l.geom_type == 'MultiLineString':
			print 'Irregular polygon path splitting feature triangle'
			#the line will start at good_pts[0] (==coord[0]), coord[1] will be next intersection
			pt3 = Point(l[0].coords[1])
		else:
			#the line will start at good_pts[0] (==coord[0]), coord[1] will be next intersection
			pt3 = Point(l.coords[1])
		return [[pt1,pt2,pt3]]
	if(len(good_pts) == 2):
		#print 'good,bad', good_pts[0],bad_pts[0]
		l1 = LineString([good_pts[0],bad_pts[0]])
		l1 = l1.intersection(poly)
		if l1.geom_type == 'MultiLineString':
			print 'Irregular polygon path splitting feature triangle'
			#the line will start at good_pts[0] (==coord[0]), coord[1] will be next intersection
			new_pt1 = Point(l1[0].coords[1])
		else:
			#the line will start at good_pts[0] (==coord[0]), coord[1] will be next intersection
			new_pt1 = Point(l1.coords[1])
		l2 = LineString([good_pts[1],bad_pts[0]])
		l2 = l2.intersection(poly)
		if l2.geom_type == 'MultiLineString':
			print 'Irregular polygon path splitting feature triangle'
			#the line will start at good_pts[0] (==coord[0]), coord[1] will be next intersection
			new_pt2 = Point(l2[0].coords[1])
		else:
			#the line will start at good_pts[0] (==coord[0]), coord[1] will be next intersection
			new_pt2 = Point(l2.coords[1])
		return [[good_pts[0],new_pt1,good_pts[1]],[good_pts[1],new_pt1,new_pt2]]


