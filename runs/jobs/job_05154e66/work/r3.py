from common import *

c = []
c += helix(32, 52)
c += anti(37, 6, 17)      # strand 6-16 with strand 21-31, hairpin
c += anti(73, 7, 16)      # strand 7-16 with strand 57-66, antiparallel
c += anti(142, 59, 67)    # strand 59-66 with strand 76-83
c += anti(175, 76, 86)    # strand 76-86 with strand 89-99
# helix laid diagonally on the face carrying V9 C11 F13 W15 / L58 F60 Y64 I66 / V78 M82
c += [(35, 7), (35, 30), (34, 30), (38, 28), (38, 9), (39, 9), (39, 64), (42, 11), (42, 62), (42, 64),
      (43, 62), (46, 13), (46, 60), (46, 62), (47, 60), (47, 80), (50, 15), (50, 58), (49, 56),
      (50, 60), (47, 82), (44, 78), (43, 78), (40, 78), (36, 66), (39, 66), (35, 66), (39, 65), (40, 65),
      (31, 34), (52, 55), (51, 54), (50, 56), (50, 57), (46, 61), (4, 32), (5, 32), (5, 31)]
# loop 67-75
c += [(67, 70), (68, 71), (70, 73), (72, 75), (68, 74), (67, 73), (72, 35), (72, 36), (72, 66), (72, 77)]
# TEST: opposite-face residues of the two outer strands (Y23 V27 V29 / Y92 Y96) brought together, parallel register
c += par(69, 23, 29)
# cross contacts over that face between non-adjacent strands (V12 V16 A61 A63 F65 F79 L81)
c += [(12, 81), (16, 79), (27, 63), (23, 61), (29, 65), (12, 94), (63, 96), (61, 92), (27, 81), (14, 92), (25, 79)]
# tail
c += [(100, 96), (99, 102), (101, 104), (103, 31), (105, 76), (107, 97), (107, 78), (104, 29)]
X, r, c = run(c, 'r3')
