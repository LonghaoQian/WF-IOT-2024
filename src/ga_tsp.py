import numpy as np
import matplotlib.pyplot as plt

def CalcDistance(p1, p2):
    return np.linalg.norm(p1 - p2)


if __name__ == '__main__':

    popSize = 2500 # population size
    genNum = 1000 # number of generation size
    mutSiteP = 0.01 # probability of exchange 2 random sites in the path (per gene, per genration)
    mutpathP = 0.02 # probability of exchange 2 peices of path (per gene, per genration)
    flipPathP = 0.08 # probability of flip random pieces of path
    
    startLocation = np.array([0, 0])
    location = np.array([[12, 4],
                         [33, 0],
                         [5, 4],
                         [9, 16],
                         [8, 0],
                         [10, 15],
                         [7, 9],
                         [4, 21],
                         [14, 10],
                         [20, 18]])/10
    
    numSites = location.shape[0]
    # distance between each site
    distM = np.zeros((numSites, numSites))
    # distance from the starting point to each site
    distS = np.zeros((numSites, 1))
    for i in range(numSites):
        distS[i] = np.linalg.norm(startLocation - location[i])
        for j in range(numSites):
            distM[i][j] = np.linalg.norm(location[i] - location[j])

    print(distM)


    # plt.figure(figsize=(5,5), dpi= 100)
    