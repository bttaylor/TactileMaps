import Image, ImageDraw
import os.path
import MapMaker3d
import json
from scipy.interpolate import interp1d

directory = MapMaker3d.directory;
imgSize = 800;

def drawMap(filename,poly,geo=True,color=(0,0,0)):

	#params = json.load(open(directory + 'ParamFiles/' + filename + '.txt'))
	params = json.load(open(directory + 'ParamFiles/' + 'v2test1431710278' + '.txt'))
	#global minlat, minlon, maxlat, maxlon, 
	minlat = params['minlat']
	minlon = params['minlon']
	maxlat = params['maxlat']
	maxlon = params['maxlon']
	#print minlon, maxlat, maxlon, minlat
	#imgSize = params['imgSize']
	maplat_img = interp1d([minlat,maxlat],[0,imgSize])
	maplon_img = interp1d([minlon,maxlon],[0,imgSize])

	filePath = directory + 'mapImages/' + filename + '.png';
	if os.path.isfile(filePath):
		mapImg = Image.open(filePath)
	else:
		mapImg = Image.new("RGB",(imgSize,imgSize),"white")

	draw_map = ImageDraw.Draw(mapImg)

	if poly.geom_type == 'MultiPolygon':
		for p in poly:
			drawPoly(draw_map,p,geo,color)

	elif poly.geom_type == "GeometryCollection":
		#print 'GeometryCollection'
		for p in poly:
			if p.geom_type == "Polygon":
				drawPoly(draw_map,p,geo,color)
			elif p.geom_type == "LineString":
				#print 'LineString'
				drawLine(draw_map,p,geo)

	elif poly.geom_type == "Polygon":
		drawPoly(draw_map,poly,geo,color)

	elif poly.geom_type == "LineString":
		drawLine(draw_map,poly,geo)
	elif poly.geom_type == "MultiLineString":
		for l in poly:
			drawLine(draw_map,l,geo)
	else:
		print 'missed',poly.geom_type

	mapImg.save(filePath)

def drawPoly(draw_map,poly,geo,color):
	minlat = 35.654
	minlon = 139.786
	maxlat = 35.668
	maxlon = 139.804
	imgSize = 800
	maplat_img = interp1d([minlat,maxlat],[0,imgSize])
	maplon_img = interp1d([minlon,maxlon],[0,imgSize])

	loc = poly.exterior.coords[0]
	for coord in poly.exterior.coords[1:]:
		if geo:
			draw_map.line((maplon_img(loc[0]),imgSize - maplat_img(loc[1]),maplon_img(coord[0]),imgSize - maplat_img(coord[1])), fill=(0,0,0))
		else:
			draw_map.line((10*loc[0],800-10*loc[1],10*coord[0],800-10*coord[1]), fill=color)
		loc = coord
	for hole in poly.interiors:
		loc = hole.coords[0]
		for coord in hole.coords[1:]:
			if geo:
				draw_map.line((maplon_img(loc[0]),imgSize - maplat_img(loc[1]),maplon_img(coord[0]),imgSize - maplat_img(coord[1])), fill=(0,0,0))
			else:
				draw_map.line((10*loc[0],800-10*loc[1],10*coord[0],800-10*coord[1]), fill=color)
			loc = coord

def drawLine(draw_map,line,geo):
	minlat = 35.64 #35.6604499
	minlon = 139.77 #139.7921080792733 
	maxlat = 35.68 #35.6676499
	maxlon = 139.83 #139.80093092072667
	imgSize = 800
	maplat_img = interp1d([minlat,maxlat],[0,imgSize])
	maplon_img = interp1d([minlon,maxlon],[0,imgSize])
	loc = line.coords[0]
	for coord in line.coords[1:]:
		if geo:
			draw_map.line((maplon_img(loc[0]),imgSize - maplat_img(loc[1]),maplon_img(coord[0]),imgSize - maplat_img(coord[1])), fill=(255,0,0))
		else:
			draw_map.line((10*loc[0],800-10*loc[1],10*coord[0],800-10*coord[1]), fill=(255,0,0),width=2)
		loc = coord

def drawMapGrid(filename,poly,X,Y,geo=True):

	#params = json.load(open(directory + 'ParamFiles/' + filename + '.txt'))
	params = json.load(open(directory + 'ParamFiles/' + 'v2test1431710278' + '.txt'))
	#global minlat, minlon, maxlat, maxlon, 
	minlat = params['minlat']
	minlon = params['minlon']
	maxlat = params['maxlat']
	maxlon = params['maxlon']
	#print minlon, maxlat, maxlon, minlat
	#imgSize = params['imgSize']
	maplat_img = interp1d([minlat,maxlat],[0,imgSize])
	maplon_img = interp1d([minlon,maxlon],[0,imgSize])

	filePath = directory + 'mapImages/' + filename + '.png';
	if os.path.isfile(filePath):
		mapImg = Image.open(filePath)
	else:
		mapImg = Image.new("RGB",(imgSize*12,imgSize*7),"white")

	draw_map = ImageDraw.Draw(mapImg)

	if poly.geom_type == 'MultiPolygon':
		for p in poly:
			loc = p.exterior.coords[0]
			for coord in p.exterior.coords[1:]:
				if geo:
					draw_map.line((maplon_img(loc[0]),imgSize - maplat_img(loc[1]),maplon_img(coord[0]),imgSize - maplat_img(coord[1])), fill=(0,0,0))
				else:
					x1 = (X-1)*imgSize + 10*loc[0];
					y1 = imgSize - (Y-1)*imgSize + 10*loc[1];
					x2 = (X-1)*imgSize + 10*coord[0];
					y2 = imgSize - (Y-1)*imgSize + 10*coord[1];
					draw_map.line((x1,y1,x2,y2), fill=(0,0,0))
				loc = coord
			for hole in p.interiors:
				loc = hole.coords[0]
				for coord in hole.coords[1:]:
					if geo:
						draw_map.line((maplon_img(loc[0]),imgSize - maplat_img(loc[1]),maplon_img(coord[0]),imgSize - maplat_img(coord[1])), fill=(0,0,0))
					else:
						x1 = (X-1)*imgSize + 10*loc[0];
						y1 = imgSize - (Y-1)*imgSize + 10*loc[1];
						x2 = (X-1)*imgSize + 10*coord[0];
						y2 = imgSize - (Y-1)*imgSize + 10*coord[1];
						draw_map.line((x1,y1,x2,y2), fill=(0,0,0))
						#draw_map.line(((X-1)*imgSize + 10*loc[0],(Y-1)*imgSize + 10*loc[1],(X-1)*imgSize + 10*coord[0],(Y-1)*imgSize + 10*coord[1]), fill=(0,0,0))
					loc = coord
	else:
		loc = poly.exterior.coords[0]
		for coord in poly.exterior.coords[1:]:
			if geo:
				draw_map.line((maplon_img(loc[0]),imgSize - maplat_img(loc[1]),maplon_img(coord[0]),imgSize - maplat_img(coord[1])), fill=(0,0,0))
			else:
				x1 = (X-1)*imgSize + 10*loc[0];
				y1 = imgSize - (Y-1)*imgSize + 10*loc[1];
				x2 = (X-1)*imgSize + 10*coord[0];
				y2 = imgSize - (Y-1)*imgSize + 10*coord[1];
				draw_map.line((x1,y1,x2,y2), fill=(0,0,0))
				#draw_map.line(((X-1)*imgSize + 10*loc[0],(Y-1)*imgSize + 10*loc[1],(X-1)*imgSize + 10*coord[0],(Y-1)*imgSize + 10*coord[1]), fill=(0,0,0))
			loc = coord
		for hole in poly.interiors:
			loc = hole.coords[0]
			for coord in hole.coords[1:]:
				if geo:
					draw_map.line((maplon_img(loc[0]),imgSize - maplat_img(loc[1]),maplon_img(coord[0]),imgSize - maplat_img(coord[1])), fill=(0,0,0))
				else:
					x1 = (X-1)*imgSize + 10*loc[0];
					y1 = imgSize - (Y-1)*imgSize + 10*loc[1];
					x2 = (X-1)*imgSize + 10*coord[0];
					y2 = imgSize - (Y-1)*imgSize + 10*coord[1];
					draw_map.line((x1,y1,x2,y2), fill=(0,0,0))
					#draw_map.line(((X-1)*imgSize + 10*loc[0],(Y-1)*imgSize + 10*loc[1],(X-1)*imgSize + 10*coord[0],(Y-1)*imgSize + 10*coord[1]), fill=(0,0,0))
				loc = coord

	mapImg.save(filePath)

def pointCount(poly):
	count = 0
	if poly.geom_type == "MultiPolygon":
		for p in poly:
			if p.geom_type == "Polygon":
				count += len(p.exterior.coords)
				for hole in p.interiors:
					count += len(hole.coords)
			elif p.geom_type == "LineString":
				print 'Line'
				count += len(p.coords)
	elif poly.geom_type == "Polygon":
		count += len(poly.exterior.coords)
		for hole in poly.interiors:
			count += len(hole.coords)
	else:
		print 'other type',poly.geom_type
	print 'Polygon Count:',count