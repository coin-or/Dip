from random import randint, seed

# The requirements for the products
DEMAND = {
    1 : 7,
    2 : 5,
    3 : 3,
    4 : 2,
    5 : 2
}

# Set of all products
PRODUCTS = list(DEMAND.keys())
PRODUCTS.sort()

# Costs of the facilities
FIXED_COST = {
    1 : 1,
    2 : 1,
    3 : 1, 
    4 : 1, 
    5 : 1
}

# Set of facilities
LOCATIONS = list(FIXED_COST.keys())
LOCATIONS.sort()

ASSIGNMENTS = [(i, j) for i in LOCATIONS for j in PRODUCTS]

seed(3)
ASSIGNMENT_COSTS = {}
for a in ASSIGNMENTS:
    ASSIGNMENT_COSTS[a] = randint(1, 10)

# The capacity of the facilities

CAPACITY = 8
