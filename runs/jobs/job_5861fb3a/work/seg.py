import numpy as np, sys
r = np.load(sys.argv[1])
segs = {'N': (1, 4), 'hN': (5, 12), 'h1a': (14, 20), 'h1b': (21, 27), 'b1': (31, 38), 'b2': (45, 51), 'b3': (54, 59),
        'b4': (60, 66), 'b5': (69, 73), 'hCa': (75, 82), 'hCb': (83, 89)}
names = list(segs)
print('min distance between segments')
print('%5s' % '' + ''.join('%6s' % n for n in names))
for a in names:
    row = ''
    for b in names:
        i0, i1 = segs[a]; j0, j1 = segs[b]
        row += '%6.1f' % r[i0 - 1:i1, j0 - 1:j1].min()
    print('%5s' % a + row)
