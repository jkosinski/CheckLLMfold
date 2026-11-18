from lib import *
import r3

REMOVE = {(88,99),(90,98),(86,103)}           # helix 95-107 register against strand 85-91 was too steep
ADD = {(86,107),(87,106),(81,107),            # 103->87, 106/107->86/85 instead
       (3,49),(3,46),(5,42),                  # helix 34-49 C-terminal half along the start of strand 3-8
       (41,56),(42,56),(25,45)}               # close in the round-3 model, middle of helix over middle of strands

def contacts():
    s = set(norm(r3.contacts())) - set(norm(REMOVE))
    return s | set(norm(ADD))

if __name__ == '__main__':
    c, r = run(contacts(), 'r4')
    lim = [(i, j) for i, j in c if j - i >= 6 and r[i-1, j-1] >= 7.8]
    print('long-range contacts held only at the limit (>=7.8 A): %d' % len(lim), lim)
    # per-residue count of predicted contacts
    cnt = np.zeros(L + 1, int)
    for i, j in c: cnt[i] += 1; cnt[j] += 1
    print('residues with <=3 contacts:', [(k, SEQ[k-1], int(cnt[k])) for k in range(1, L + 1) if cnt[k] <= 3])
    print('residues with >=11 contacts:', [(k, SEQ[k-1], int(cnt[k])) for k in range(1, L + 1) if cnt[k] >= 11])
