'''
Created by Michael Jendryke 2016
0 is Inputfile
1 is Height field
2 is Weight field
3 is Multiplier Min
4 is Multiplier Max
5 is Minimum Total Points Min
6 is Minimum Total Points Max
7 is Output folder
'''



import csv
import numpy as np
import networkx as nx
import arcpy as ap
import arcpy.mapping as am
import os


def main():
    inFC = ap.GetParameterAsText(0)
    heigthField = ap.GetParameterAsText(1)
    weightField = ap.GetParameterAsText(2)
    mtpMin = int(ap.GetParameterAsText(3))
    mtpMax = int(ap.GetParameterAsText(4))
    multMin = float(ap.GetParameterAsText(5))
    multMax = float(ap.GetParameterAsText(6))
    outdir = ap.GetParameterAsText(7)


    # Is data projected?
    desc = ap.Describe(inFC)
    ap.AddMessage(desc.SpatialReference.type)
    if desc.SpatialReference.type != 'Projected':
        ap.AddError('Your data does not seem to be projected')
        quit()


    # Get data to work with
    data = []
    for row in ap.da.SearchCursor(inFC, ["OID@", "SHAPE@XY", heigthField,weightField]):
        id = int(row[0])
        x, y = row[1]
        z = row[2]
        w = row[3]
        data.append([id, x, y, z, w]);


    # Loop through mtp and weight, then export Shape and add it to
    for m in range(mtpMin, mtpMax+1, 1):
        for w in range(int(multMin*10), int(multMax*10)+1, 1):
            r = archsphere(data, m, w/10)
            outFC = createshape(inFC, outdir, data, m, w)
    ap.AddMessage("The End.")


def archsphere(Data,mtp,multiplier):
    combinations = soapbubbles(Data, multiplier) #second is the multiplier
    cluster = magic2(combinations)
    result = addclustertodata(Data, cluster, mtp)
    return result


def createshape(infile, outdir, d, m, w):
    outfile = os.path.join(outdir, "Result_m%s_w%s.shp" % (str(m), str(w)))
    ap.AddMessage("Creating file: %s" % outfile)
    ap.Copy_management(infile, outfile, "Shapefile")
    ap.AddField_management(outfile, "cluster", "SHORT","8")
    cursor = ap.da.UpdateCursor(outfile, ['OID@', 'cluster'])
    # Update the road buffer distance field based on road type.
    # Road type is either 1,2,3,4  Distance is in meters.
    for row in cursor:
        # ap.AddMessage("{0}, {1}".format(row[0],d[row[0]][5])) #column 5 holds the cluster
        row[1] = d[row[0]][5]
        cursor.updateRow(row)

    # Delete cursor and row objects
    addlayertotoc(outfile)
    del cursor, row


def addlayertotoc(o):
    # Set up the dataframe to display the shapes
    mxd = am.MapDocument("CURRENT")
    # get the data frame
    df = am.ListDataFrames(mxd, "*")[0]
    # create a new layer
    newlayer = am.Layer(o)
    # add the layer to the map at the bottom of the TOC in data frame 0
    ap.mapping.AddLayer(df, newlayer, "BOTTOM")
    ap.RefreshActiveView()
    ap.RefreshTOC()
    del mxd, df, newlayer


def magic2(mat):
    G = nx.from_numpy_matrix(np.array(mat))
    G = nx.connected_components(G)
    return G


def soapbubbles(data, multiplier):
    #print('rows in data: ' + str(np.size(data,0)))
    nbrr = np.size(data,0)
    #nbrr = 10
    combinations = np.zeros([nbrr,nbrr])
    for i in range(nbrr):
        for j in range(i,nbrr,1):
            combinations[i][j] = 1
            if i !=j :
                #print('combining: i=' + str(i) + ' with: ' + str(j))
                distance = distance3D(data[i][1],data[i][2],data[i][3],data[j][1],data[j][2],data[j][3])
                #print('distance: ' + str(distance))
                twobubbles = bubbles(data, multiplier, i, j) #distance of two bubble together
                #print('twobubbles: ' + str(twobubbles))
                if distance < twobubbles:
                    combinations[i][j] = 1
                else:
                    combinations[i][j] = 0
    return np.array(combinations)


def bubbles(pts, ratio, ids, idd):  # bubble diameter for ids(ource) and idd(estination)
    # get max diameter of ids
    if float(pts[ids][3]) > float(pts[ids][4]):
        mds = float(pts[ids][3])
    else:
        mds = float(pts[ids][4])
    # get max diameter of idd
    if float(pts[idd][3]) > float(pts[idd][4]):
        mdd = float(pts[idd][3])
    else:
        mdd = float(pts[idd][4])
    b = (pts[ids][4] * ratio) + (pts[idd][4] * ratio)
    return b


def distance3D(xs,ys,zs,xd,yd,zd): #coordinates from s(ource) to d(estination)
    d = np.sqrt((xs - xd) ** 2 + (ys - yd) ** 2 + (zs - zd) ** 2)
    return d

	
def addclustertodata(d,c,m):
    clusterID = 0
    totalids = 0
    totalnoise = 0
    for cc in c:
        #print(len(cc))
        totalids = totalids + len(cc)
        if len(cc) <= m:
            marker = 0
            totalnoise = totalnoise + len(cc)
        else:
            clusterID += 1
            marker = clusterID
        #clusterID += 1
        #marker = clusterID
        for r in cc:
            d[r] = d[r] + [marker]
            #print(d[r])
    ap.AddMessage(' number of clusters: ' + str(clusterID) + ' noise: ' + str(totalnoise))
    ap.AddMessage(' and % 6.2f percent are clustered' % ((1-(totalnoise/totalids))*100))
    return d


main()
