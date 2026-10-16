from random import randint

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
    1 : 10,
    2 : 20,
    3 : 16, 
    4 : 1, 
    5 : 2
}

# Set of facilities
LOCATIONS = list(FIXED_COST.keys())
LOCATIONS.sort()

# The capacity of the facilities

CAPACITY = 8
