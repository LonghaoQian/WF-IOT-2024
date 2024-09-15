import numpy as np
import ga_lib as gl
import matplotlib.pyplot as plt

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
    # children
    child_genes = np.zeros((popSize, numSites))
    for i in range(numSites):
        distS[i] = np.linalg.norm(startLocation - location[i])
        for j in range(numSites):
            distM[i][j] = np.linalg.norm(location[i] - location[j])

    # print(distM)
    # generate genes
    genes = np.zeros((popSize, numSites))
    for i in range(popSize):
        genes[i] = np.random.permutation(numSites)

    # total path length
    pathLength = np.zeros((popSize, 1))
    # gene probability to be selected
    probP = np.zeros((popSize, 1))
    for ig in range(genNum):
        # step 1, find total path distance including the start location
        for ip in range(popSize):
            gc = genes[ip]
            pt = 0 # total distance
            for j in range(numSites - 1):
                pt = pt + distM[gc[j]][gc[j + 1]]
            pt = pt + distS[gc[0]] + distS[gc[-1]]
            probP[ip] = 1 / pt
        pp = np.sum(probP)
        probP = probP / pp
        # find the best plan
        best_gene = genes[genes.index(max(probP))]
        # perform crossover
        cross_over_indices = gl.get_roulette_wheel_indexes(popSize, probP)
        # reset children
        child_genes.fill(0)



    # print(genes)
    # plt.figure(figsize=(5,5), dpi= 100)
    