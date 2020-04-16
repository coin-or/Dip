from .cutting_stock_func import cross

length = {
"9cm": 9,
"7cm": 7,
"5cm": 5
}

ITEMS = list(length.keys())

demand = {
"9cm": 3,
"7cm": 2,
"5cm": 2
}

total_patterns = sum(demand[i] for i in ITEMS)

total_length = {}
for p in range(total_patterns):
    total_length[p] = 20
PATTERNS = list(total_length.keys())

CUTS = cross(PATTERNS, ITEMS)

