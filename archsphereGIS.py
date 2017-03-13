'''
Created by Michael Jendryke 2016
http://stackoverflow.com/questions/41112073/point-cloud-cluster-analysis-in-python-identifying-clusters-from-binary-matrix
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
#import shapefile
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

def magic(mat):
    # Make the undirected version of the graph (no self loops)
    A = (mat + mat.T) * (1 - np.eye(mat.shape[0]))
    # Make the degree matrix
    D = np.diag(A.sum(axis=1) + A.sum(axis=0)) / 2
    # thats all we need to define the laplacian
    L = D - A

    # The number of zeros eigenvalues of the Laplacian is exactly the number of CCs
    np.isclose(np.linalg.eigvals(L), 0).sum()

    # The connected compoments themselves are identified by rows that have the same nullspace vector
    u, s, vh = np.linalg.svd(L)
    ns = vh[(s >= 1e-13).sum():].conj().T

    # the following is a little numpy trick to find unique rows
    # chopping off the last few decimal places to account for numerical errors
    ns_ = np.ascontiguousarray(np.round(ns, 8)).view(np.dtype((np.void, ns.dtype.itemsize * ns.shape[1])))
    ns_basis, row_to_cc_id = np.unique(ns_, return_inverse=True)
    # Finally we can just use this to convert to the desired output format
    groups = [[] for _ in range(len(ns_basis))]
    for row, id in enumerate(row_to_cc_id):
        groups[id].append(row)
    return groups


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


def soapfoam(s):
    r = 0
    idx = []
    while r < np.size(s,0):
        print('row: ' + str(r) + '\t\t' + str(s[r]))
        row = []
        for i in range(np.size(s,1)):
            if s[r][i] == 1:
                row = row + [i]
        r = r + 1
        #print(row)
        idx = idx + [row]
    #print(idx)
    return idx


def distance3D(xs,ys,zs,xd,yd,zd): #coordinates from s(ource) to d(estination)
    d = np.sqrt((xs - xd) ** 2 + (ys - yd) ** 2 + (zs - zd) ** 2)
    return d


def writeData(data, clusterlists, ratio, minpts):
    file = 'E:\\Dropbox\\DATA\\research\\Papers(int)\\Archaeological_Science\\DBSCAN\\DBSCAN_clone\\results\\SBC_ratio' + str(
        ratio) + '_minpts' + str(minpts) + '.csv'

    with open(file, 'a') as fileout:
        cluster = 0
        for row in clusterlists:
            print(row)
            for ID in row:
                for item in data[ID]:
                    fileout.write("{},".format(item))
                if np.size(row, 0) < minpts:
                    fileout.write(str(9999))
                    cluster = cluster - 1
                else:
                    if cluster < 0:
                        cluster = 1
                    fileout.write(str(cluster))
                fileout.write("\n")

            cluster = cluster + 1
        fileout.close()
    return


def getData(dataPath):
    Data = []

    with open(dataPath, 'r') as filein:
        reader = csv.reader(filein)
        c = 0
        for row in reader:
            if c == 0:
                c = c + 1
                continue
            else:
                # row = re.split(r'\t+',row[0])
                Data.append([float(row[0]), float(row[1]), float(row[2]), float(row[3]), float(row[4]), float(row[5]),
                             str(row[6])])
        filein.close()
    return [Data]


def parse(line):
    data = line.split(" ")
    return [int(data[0]), int(data[1])]


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


# def exportshape(points, mtp, weight):
#     w = shapefile.Writer()
#     w.shapeType = 1 #see 'shape type' at  http://www.esri.com/library/whitepapers/pdfs/shapefile.pdf
#     #w.autoBalance = 1
#     w.field('Object')
#     w.field('Z')
#     w.field('NS')
#     w.field('WE')
#     w.field('Height')
#     w.field('cluster')
#
#
#     for p in points:
#         #print(p)
#         w.point(p[0], p[1],p[2],0)
#         w.record(p[6], p[2],p[3],p[4],p[5],p[7])
#
#
#     w.save('result_mtp_' + str(mtp) + '_weight_' + str(weight*10))
#
#     print(' done.')

main()
