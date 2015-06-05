
from shapely.geometry import box, Polygon, LineString, Point, MultiLineString, MultiPolygon
import math
import elevation2stl

# def makeDome(r,center,coarseness):

# 	(tri,outer) = makeCircle(r,center,coarseness);
# 	for i in range(0,len(tri)):
# 		    p1 = [tri[i][0].x,tri[i][0].y,tri[i][0].z]
# 		    p2 = [tri[i][1].x,tri[i][1].y,tri[i][1].z]
# 		    p3 = [tri[i][2].x,tri[i][2].y,tri[i][2].z]
# 		    elevation2stl.printTriangle(p1,p2,p3,fname)

def makeCircle(r,center,coarseness,z_adj=0):
	inner_pts = [Point(center.x, center.y, r + z_adj)];
	#print 'Row 0: ' 
	#print inner

	triangles = [];

	row_cnt = coarseness + 1;
	phi = (math.pi/2)/coarseness;

	for i in range(0,coarseness):
		#inner_r = r * math.sin(i * phi);
		#inner_z = r * math.cos(i * phi);
		outer_r = r * math.sin((i+1) * phi);
		outer_z = r * math.cos((i+1) * phi);

		tri_row = i + 1;
		#print 'Row ' + str(tri_row) + ': '

		#calc points for inner row
		theta = (math.pi/2)/tri_row;

		#first point (along the axis)
		x1 = center.x + outer_r * math.sin(0 * theta)
		y1 = center.y + outer_r * math.cos(0 * theta)
		outer_pts = [Point(x1,y1,outer_z + z_adj)];
		#print [x1, y1, outer_z]
		for j in range(0,tri_row):
			#next j points
			x2 = center.x + outer_r * math.sin((j+1) * theta);
			y2 = center.y + outer_r * math.cos((j+1) * theta);
			outer_pts.append(Point(x2,y2,outer_z + z_adj))
			#print [x2, y2, outer_z]

		triangles.append([inner_pts[0],outer_pts[0],outer_pts[1]])

		tri_count = 1 + (len(inner_pts)-1)*2;
		for j in range(1,len(inner_pts)):
			#one quandrant
			triangles.append([inner_pts[j-1],	outer_pts[j],	inner_pts[j]]);
			triangles.append([inner_pts[j],		outer_pts[j],	outer_pts[j+1]])

		inner_pts = outer_pts;
	all_outer = outer_pts[0:len(outer_pts)-1];

	rotate = math.pi/2;
	inner_pts = [Point(center.x, center.y, r + z_adj)];
	for i in range(0,coarseness):
		#inner_r = r * math.sin(i * phi);
		#inner_z = r * math.cos(i * phi);
		outer_r = r * math.sin((i+1) * phi);
		outer_z = r * math.cos((i+1) * phi);

		tri_row = i + 1;
		#print 'Row ' + str(tri_row) + ': '

		#calc points for inner row
		theta = (math.pi/2)/tri_row;

		#first point (along the axis)
		x1 = center.x + outer_r * math.sin(0 * theta + rotate)
		y1 = center.y + outer_r * math.cos(0 * theta + rotate)
		outer_pts = [Point(x1,y1,outer_z + z_adj)];
		#print [x1, y1, outer_z]
		for j in range(0,tri_row):
			#next j points
			x2 = center.x + outer_r * math.sin((j+1) * theta + rotate);
			y2 = center.y + outer_r * math.cos((j+1) * theta + rotate);
			outer_pts.append(Point(x2,y2,outer_z + z_adj))
			#print [x2, y2, outer_z]

		triangles.append([inner_pts[0],outer_pts[0],outer_pts[1]])

		tri_count = 1 + (len(inner_pts)-1)*2;
		for j in range(1,len(inner_pts)):
			#one quandrant
			triangles.append([inner_pts[j-1],	outer_pts[j],	inner_pts[j]]);
			triangles.append([inner_pts[j],		outer_pts[j],	outer_pts[j+1]])

		inner_pts = outer_pts;
	all_outer = all_outer + outer_pts[0:len(outer_pts)-1];
	#there is overlap w/ the points and the mirror technique places the duplicate in
	#either first or last index depending on how it was mirrored

	rotate = math.pi
	inner_pts = [Point(center.x, center.y, r + z_adj)];
	for i in range(0,coarseness):
		#inner_r = r * math.sin(i * phi);
		#inner_z = r * math.cos(i * phi);
		outer_r = r * math.sin((i+1) * phi);
		outer_z = r * math.cos((i+1) * phi);

		tri_row = i + 1;
		#print 'Row ' + str(tri_row) + ': '

		#calc points for inner row
		theta = (math.pi/2)/tri_row;

		#first point (along the axis)
		x1 = center.x + outer_r * math.sin(0 * theta + rotate)
		y1 = center.y + outer_r * math.cos(0 * theta + rotate)
		outer_pts = [Point(x1,y1,outer_z + z_adj)];
		#print [x1, y1, outer_z]
		for j in range(0,tri_row):
			#next j points
			x2 = center.x + outer_r * math.sin((j+1) * theta + rotate);
			y2 = center.y + outer_r * math.cos((j+1) * theta + rotate);
			outer_pts.append(Point(x2,y2,outer_z + z_adj))
			#print [x2, y2, outer_z]

		triangles.append([inner_pts[0],outer_pts[0],outer_pts[1]])

		tri_count = 1 + (len(inner_pts)-1)*2;
		for j in range(1,len(inner_pts)):
			#one quandrant
			triangles.append([inner_pts[j-1],	outer_pts[j],	inner_pts[j]]);
			triangles.append([inner_pts[j],		outer_pts[j],	outer_pts[j+1]])

		inner_pts = outer_pts;
	all_outer = all_outer + outer_pts[0:len(outer_pts)-1];

	rotate = 3*math.pi/2;
	inner_pts = [Point(center.x, center.y, r + z_adj)];
	for i in range(0,coarseness):
		#inner_r = r * math.sin(i * phi);
		#inner_z = r * math.cos(i * phi);
		outer_r = r * math.sin((i+1) * phi);
		outer_z = r * math.cos((i+1) * phi);

		tri_row = i + 1;
		#print 'Row ' + str(tri_row) + ': '

		#calc points for inner row
		theta = (math.pi/2)/tri_row;

		#first point (along the axis)
		x1 = center.x + outer_r * math.sin(0 * theta + rotate)
		y1 = center.y + outer_r * math.cos(0 * theta + rotate)
		outer_pts = [Point(x1,y1,outer_z + z_adj)];
		#print [x1, y1, outer_z]
		for j in range(0,tri_row):
			#next j points
			x2 = center.x + outer_r * math.sin((j+1) * theta + rotate);
			y2 = center.y + outer_r * math.cos((j+1) * theta + rotate);
			outer_pts.append(Point(x2,y2,outer_z + z_adj))
			#print [x2, y2, outer_z]

		triangles.append([inner_pts[0],outer_pts[0],outer_pts[1]])

		tri_count = 1 + (len(inner_pts)-1)*2;
		for j in range(1,len(inner_pts)):
			#one quandrant
			triangles.append([inner_pts[j-1],	outer_pts[j],	inner_pts[j]]);
			triangles.append([inner_pts[j],		outer_pts[j],	outer_pts[j+1]])

		inner_pts = outer_pts;
	all_outer = all_outer + outer_pts[0:len(outer_pts)-1];

	return (triangles, all_outer)


