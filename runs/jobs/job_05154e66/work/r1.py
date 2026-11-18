from common import *

c = []
c += helix(32, 52)
c += anti(37, 6, 17)      # strand 6-16 with strand 21-31, hairpin
c += anti(73, 7, 16)      # strand 7-16 with strand 57-66, antiparallel
c += anti(142, 59, 67)    # strand 59-66 with strand 76-83
c += anti(175, 76, 86)    # strand 76-86 with strand 89-99
# helix laid on one face of the strands
c += [(35, 7), (35, 30), (35, 66), (39, 9), (39, 64), (39, 28), (42, 11), (42, 62), (43, 26), (46, 13), (46, 60),
      (47, 24), (50, 15), (50, 58), (49, 56), (31, 34), (30, 33), (52, 55), (53, 56), (51, 54)]
# loops
c += [(67, 70), (68, 71), (70, 73), (72, 75), (98, 101), (99, 102), (101, 104)]
X, r, c = run(c, 'r1')
