from lib import *
L = 89
c = []
c += helix(4, 12)
c += helix(14, 27)
c += helix(73, 88)
# hairpin 31-38 / 45-52
c += ladder_anti(32, 51, 8)   # 32-51 ... 39-44
# hairpin 45-51 / 54-60
c += ladder_anti(45, 60, 7)   # 45-60 ... 51-54
# hairpin 60-66 / 69-72
c += ladder_anti(62, 73, 4)   # 62-73.. 65-70
# helix 14-26 packs on sheet
c += [(26,31),(23,33),(23,35),(19,35),(22,34),(19,48),(23,48),(18,36),(15,37)]
# C-term helix packs onto sheet/helix
c += [(76,61),(79,61),(83,64),(76,58),(80,47),(83,45),(86,18),(83,22),(79,24),(86,14)]
# N-term helix
c += [(6,18),(7,22),(6,14),(12,60),(7,58)]
c = clean(c, L)
print(len(c))
X = fold(L, c); r = report(X, c)
print(near(r, c)[:60])
