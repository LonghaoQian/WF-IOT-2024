import numpy as np
import matplotlib.pyplot as plt

def CalcDistance(p1, p2):
    return np.linalg.norm(p1 - p2)


if __name__ == '__main__':
    # define the starting location
    start = np.array([0, 0])
    num = 10 # number of sites
    dist = np.zeros((num + 1, num + 1))
    location = np.zeros((num, 2))
    # load the locations
    location[0, :] = np.array([10, 12])
    

    print(location)

    # plt.figure(figsize=(5,5), dpi= 100)
    