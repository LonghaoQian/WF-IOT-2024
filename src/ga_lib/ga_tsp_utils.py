import numpy as np
import random


def get_site_dist(p1, p2):
    return np.linalg.norm(p1 - p2)


def get_roulette_wheel_indexes(m, prob_list):
    # get m samples of indices form probList
    cpr = np.cumsum(prob_list)
    n = prob_list.size
    return random.choices(range(n), prob_list, m)


def insert_at_beginning(g1, g2, cp):
    # insert gene from g2 (0 : cp) to the beginning of g1(0: g1)
    res = g1
    for i in range(cp):
        if g2[i] != res[i]:
            idx = np.where(res == g2[i])
            # swap the elements
            res[idx] = res[i]
            res[i] = g2[i]

    return res