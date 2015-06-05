import urllib
import MapMaker3d

def downloadMap(lon1,lat1,lon2,lat2,filename):
	f = MapMaker3d.directory + 'mapFiles/' + filename
	api_string = 'http://api.openstreetmap.org/api/0.6/map?bbox='
	d=urllib.urlretrieve(api_string + str(lon1) + "," + str(lat2) + 
			"," + str(lon2) + "," + str(lat1), f + '.osm')
	print d[0]

