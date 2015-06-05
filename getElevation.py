
import simplejson
import urllib
import MapMaker3d
import json

#samples = MapMaker3d.samples
ratio = 1  #force square
ELEVATION_BASE_URL = 'https://maps.googleapis.com/maps/api/elevation/json'
CHART_BASE_URL = 'http://chart.apis.google.com/chart'

def getElevation(lon1,lat1,lon2,lat2,filename):
	print MapMaker3d.directory
	params = json.load(open(MapMaker3d.directory + 'ParamFiles/' + filename + '.txt'))
	samples = params['samples']
	print '***' + str(samples)
	#gets Elevation data and saves as a .csv file
	#Samples evenly along straight lines from lon1-lon2 and lat1-lat2
	#each row of the csv file is along a single latitude
	#there are N rows (where N = samples)
	#
	#filename.csv format:
	#index, latitude1, longitude1, latitude2, longitude2, elevation1,...,elevationN

	#store the files in the proper folder
	filename = MapMaker3d.directory + "elevationFiles/" + filename;

  	for i in range(0,int(samples*ratio)):
	    print "i: " +str(i)
	    # interpolate latitude for the i'th row
	    lat = lat1 - (float(i)*(lat1-lat2))/(samples - 1)
	    # format string for the path along the row 
	    row = str(lat) + "," + str(lon1) + "|" + str(lat) + "," + str(lon2);
	    # save index, lat, lon info
	    f = open(filename + '.csv','a')
	    f.write(str(i) + "," + str(lat) + "," + str(lon1) + "," + str(lat) + "," + str(lon2))
	    f.close()
	    # get the elvation along the row
	    error = True
	    while(error):
	    	error = getElevationPath(filename,row,str(samples))


def getElevationPath(filename,path="40.414116,-79.988460|40.702337,-80.384687",samples="100", **elvtn_args):
	elvtn_args.update({
	  'path': path,
	  'samples': samples
	})

	url = ELEVATION_BASE_URL + '?' + urllib.urlencode(elvtn_args)
	response = simplejson.load(urllib.urlopen(url))

	# Create a dictionary for each results[] object
	elevationArray = []

	for resultset in response['results']:
		elevationArray.append(resultset['elevation'])

	# Create the chart passing the array of elevation data
	return getChart(filename,chartData=elevationArray,samples=float(samples))


def getChart(filename,chartData,samples=20, chartDataScaling="-500,5000", chartType="lc",chartLabel="Elevation in Meters",chartSize="500x160",chartColor="orange", **chart_args):
    chart_args.update({
      'cht': chartType,
      'chs': chartSize,
      'chl': chartLabel,
      'chco': chartColor,
      'chds': chartDataScaling,
      'chxt': 'x,y',
      'chxr': '1,-500,5000'
    })

    dataString = 't:' + ','.join(str(x) for x in chartData)\

    f = open(filename + '.csv','a')
    print len(chartData)
    if len(chartData) == samples:
		for x in chartData:
			f.write(',' + str(x))
		f.write('\n')
		return False

    return True

    #chart_args['chd'] = dataString.strip(',')

    #chartUrl = CHART_BASE_URL + '?' + urllib.urlencode(chart_args)
