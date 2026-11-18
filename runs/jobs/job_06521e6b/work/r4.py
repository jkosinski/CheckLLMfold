from lib import *
import json, sys
tag = sys.argv[1] if len(sys.argv) > 1 else 'r4'
c = json.load(open('r3.json'))
c = [tuple(x) for x in c]
add = [(7, 28), (11, 34), (66, 102), (65, 107), (50, 73), (71, 77), (72, 88), (12, 15), (12, 16), (23, 26), (24, 27),
       (24, 28), (83, 86), (83, 87), (104, 108), (105, 108), (43, 50), (44, 50), (11, 51), (11, 52), (50, 74),
       (2, 5), (3, 6), (1, 4)]
rem = [(68, 81), (58, 107), (4, 107), (4, 104), (6, 103)]
c = [x for x in c if x not in rem] + add
X, r, c = run(c, tag)
json.dump(c, open(tag + '.json', 'w'))
