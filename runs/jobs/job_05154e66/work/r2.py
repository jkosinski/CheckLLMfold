from common import *

c = []
c += helix(32, 52)
c += anti(37, 6, 17)      # strand 6-16 with strand 21-31, hairpin
c += anti(73, 7, 16)      # strand 7-16 with strand 57-66, antiparallel
c += anti(142, 59, 67)    # strand 59-66 with strand 76-83
c += anti(175, 76, 86)    # strand 76-86 with strand 89-99
# helix laid diagonally on the face carrying V9 C11 F13 W15 / L58 F60 Y64 I66 / V78 M82
c += [(35, 7), (35, 30), (34, 30), (38, 28), (38, 9), (39, 9), (39, 64), (39, 28), (42, 11), (42, 62), (42, 64),
      (43, 62), (46, 13), (46, 60), (46, 62), (47, 60), (47, 80), (50, 15), (50, 58), (49, 58), (49, 56),
      (50, 60), (51, 82), (47, 82), (44, 78), (43, 78), (40, 76), (36, 66), (39, 66), (35, 9),
      (31, 34), (30, 33), (52, 55), (53, 56), (51, 54), (50, 56), (50, 57)]
# loop 67-75: V71 / I72 tucked against the start of the helix and F65 / F79
c += [(67, 70), (68, 71), (70, 73), (72, 75), (68, 74), (72, 35), (71, 32), (72, 65), (71, 36), (72, 77)]
# C-terminal tail folded back over the opposite face (V12 V16 A61 A63 F65 F79 L81 Y92 Y96)
c += [(98, 101), (99, 102), (101, 104), (100, 96), (102, 94), (103, 79), (104, 92), (104, 81), (106, 81),
      (107, 63), (107, 81), (107, 12), (106, 61), (105, 79), (103, 96)]
# opposite-face residues of the outer strands brought towards each other (Y23 V27 V29 with tail)
c += [(107, 14), (105, 63)]
c += [(2, 5), (3, 6), (4, 31), (4, 33)]
X, r, c = run(c, 'r2')
