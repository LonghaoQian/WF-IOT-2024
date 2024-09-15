import numpy as np
import ga_lib as gl
import matplotlib.pyplot as plt
import random

def perform_cross_over(popSize, popProb, num_of_sites, genes):
    cross_over_indices = gl.get_roulette_wheel_indexes(popSize, popProb)
    child_genes = np.zeros((popSize, numSites))
    for i in range(popSize / 2):
        # perform cross_over for each pair of parent genes
        i1 = 2 * i
        i2 = 2 * i + 1
        # get the index of genes
        i1_g = cross_over_indices[i1]
        i2_g = cross_over_indices[i2]
        # generate a crossover point
        cp = random.randint(0,  num_of_sites - 2)
        # generate 2 children
        child_genes[i1] = gl.insert_at_beginning(genes[i1_g], genes[i2_g], cp)
        child_genes[i2] = gl.insert_at_beginning(genes[i2_g], genes[i1_g], cp)
    return child_genes


def perform_mutation(popSize, num_of_sites, genes, mute_sites, mute_path, mute_flip):
    # 1. randomly exchange 2 sites
    for i in range(popSize):
        if random.uniform(0.0, 1.0) < mute_sites:
            # random number of sites
            rns = random.randint(0, num_of_sites - 1)
            # use the first chunk of the random sequence
            rnss = np.random.permutation(num_of_sites)
            ctp = rnss[0 : rns]
            g_list = genes[i].tolist()
            gt = [g_list[j] for j in ctp]
            # hash the sites
            sp = np.random.permutation(rns)
            gt2 = [gt[j] for j in sp]
            for j in range(ctp):
                genes[i][j] = gt2[j]

    # 2. randomly exchange 2 segments of path
    for i in range(popSize):
        if random.uniform(0.0, 1.0) < mute_path:
            cp = random.randint(1, num_of_sites - 2)
            g1 = genes[i][0: cp]
            g2 = genes[i][cp + 1:]
            genes[i] = np.concatenate(g2, g1)

    # 3. randomly flip a piece of path
    for i in range(popSize):
        if random.uniform(0.0, 1.0) < mute_flip:
            n1 = random.randint(0, num_of_sites - 1)
            n2 = random.randint(n1, num_of_sites - 1)
            genes[i][n1 : n2 + 1] = np.flip(genes[i][n1 : n2 + 1])


if __name__ == '__main__':

    popSize = 2500 # population size (should be an even number)
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

    # initiation plot
    plt.ion()  # turning interactive mode on
    fig, ax = plt.subplots()
    ax.axis('equal')
    for i in range(numSites):
        ax.scatter(location[i][0], location[i][1], c='g')

    fig.canvas.draw()
    fig.canvas.flush_events()

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
        # choose parents to crossover
        # set children to population
        genes = perform_cross_over(popSize, probP, numSites, genes)
        # perform mutation
        perform_mutation(popSize, numSites, genes, mutSiteP, mutpathP, flipPathP )
        # elitism
        genes[0] = best_gene


        # print(genes)
    # plt.figure(figsize=(5,5), dpi= 100)
    