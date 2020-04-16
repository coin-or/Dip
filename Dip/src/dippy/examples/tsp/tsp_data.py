from math import sqrt

# 2d Euclidean TSP

# x,y coords of cities
CITY_LOCS = [(0, 2), (0, 4), (1, 2), (1, 4), \
             (4, 1), (4, 4), (4, 5),  (5, 0), \
             (5, 2), (5, 5)]
CITIES = list(range(len(CITY_LOCS)))

ARCS = [] # list of arcs (no duplicates)
ARC_COSTS = {} # distance

# for each city, list of arcs into/out of
CITY_ARCS = [[] for i in CITIES]

# use 2d euclidean distance
def dist(x1, y1, x2, y2):
    return sqrt((x1-x2)**2 + (y1-y2)**2)

# construct list of arcs
for i in CITIES:
    i_x, i_y = CITY_LOCS[i]
    for j in CITIES[i+1:]:
        j_x, j_y = CITY_LOCS[j]
        ARC_COSTS[(i,j)] = dist(i_x, i_y, j_x, j_y)
        ARCS.append((i, j))
        CITY_ARCS[i].append((i, j))
        CITY_ARCS[j].append((i, j))

