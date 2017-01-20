'''
Created by Michael Jendryke 2016
http://stackoverflow.com/questions/41112073/point-cloud-cluster-analysis-in-python-identifying-clusters-from-binary-matrix
'''

import csv
import numpy as np
import networkx as nx
import shapefile


dataPath   = '/Users/mac/Documents/Michael/DBSCAN/DBSCAN_clone/input/data.csv'
#dataPath   = 'input/data.csv'

def main():
    print("Soap Bubble Clustering")
    #print(configPath)
    Data = getData()[0]
    #Data = [1,1,1,1,1],[1,1,2,1,1],[2,2,2,1,1],[3,3,3,1,1],[4,4,4,1,1],[5,5,5,1,1],[50,50,50,1,1],[95,95,95,1,1],[96,96,96,1,1],[97,97,97,1,1],[98,98,98,1,1],[99,99,99,1,1],[2,2,3,1,1],[2,2,1,1,1],[2,2,4,1,1]
    print("These are the input points")
    print(Data)
    multiplier = 3.0
    combinations = soapbubbles(Data, multiplier) #second is the ratio
    print('These are all the combinatioins: ')
    print(combinations)

    print("Let the magic happen")
    #cluster = magic(combinations)
    #print(cluster)
    cluster = magic2(combinations)
    print('good')
    data = addclustertodata(Data,cluster)
    exportshape(data, multiplier)


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


def soapbubbles(data, ratio):
    print('rows in data: ' + str(np.size(data,0)))
    nbrr = np.size(data,0)
    #nbrr = 10
    combinations = np.zeros([nbrr,nbrr])
    for i in range(nbrr):
        for j in range(i,nbrr,1):
            combinations[i][j] = 1
            if i !=j :
                #print('combining: i=' + str(i) + ' with: ' + str(j))
                distance = distance3D(data[i][0],data[i][1],data[i][2],data[j][0],data[j][1],data[j][2])
                #print('distance: ' + str(distance))
                twobubbles = bubbles(data, ratio, i, j) #distance of two bubble together
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
    b = (mds * ratio) + (mdd * ratio)
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


def getData():
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


def addclustertodata(d,c):
    clusterID = 0
    totalids = 0
    for cc in c:
        #print(len(cc))
        totalids = totalids + len(cc)
        if len(cc) == 1:
            marker = 0
        else:
            clusterID += 1
            marker = clusterID
        #clusterID += 1
        #marker = clusterID
        for r in cc:
            d[r] = d[r] + [marker]
            #print(d[r])
    print(totalids)
    return d


def exportshape(points,m):
    w = shapefile.Writer()
    w.shapeType = 1 #see 'shape type' at  http://www.esri.com/library/whitepapers/pdfs/shapefile.pdf
    #w.autoBalance = 1
    w.field('Object')
    w.field('Z')
    w.field('NS')
    w.field('WE')
    w.field('Height')
    w.field('cluster')


    for p in points:
        print(p)
        w.point(p[0], p[1],p[2],0)
        w.record(p[6], p[2],p[3],p[4],p[5],p[7])


    w.save('result_mult' + str(m*10))

    print('done.')

main()
