from lib import *
import pickle
cs = []
# hairpin: strand 1-7 / strand 11-17, turn 8-10
cs += anti(18, 1, 7, (-2, -1, 0, 1, 2))
cs += [(7, 10), (8, 11), (6, 9)]
# helix 23-34
cs += helix(22, 35)
# strand 40-46 with 66-72 antiparallel
cs += anti(112, 40, 46, (-1, 0, 1))
# hairpin 43-45 / 48-50
cs += anti(93, 43, 45, (-1, 0, 1))
# parallel 1-7 with 65-71
cs += par(64, 1, 7, (-1, 0, 1))
# turn / short helix 56-60, local turns
cs += helix(55, 61) + [(51, 54), (52, 55), (62, 65), (18, 21), (19, 22), (36, 39), (37, 40), (37, 41), (72, 75)]
# helix hydrophobic face onto strands and other hydrophobics
cs += [(3, 26), (5, 30), (15, 26), (17, 23), (15, 30), (26, 43), (30, 43), (30, 41), (27, 43), (23, 56), (26, 56),
       (26, 67), (30, 69), (56, 61), (61, 67), (50, 59), (43, 50), (3, 61), (17, 56), (23, 61), (5, 67), (3, 23)]
cs, X, r = run(cs, 'round1'); pickle.dump(cs, open('c1.pkl', 'wb'))
