"""Idealised CB layout from my secondary-structure assignment, used only to write a geometrically
self-consistent contact list (strands 4.8 A apart, 3.3 A rise; helices 10 A over the sheet, 1.5 A rise)."""
import numpy as np

def layout(S12=84.0, S23=105.0, S34=119.5, h1_phase=120.0, hC_phase=300.0, hC_dir=(-0.45, 0.89), hC_start=(15.0, -2.0),
           h1_x=2.4, h1_y26=0.0, hN=True, hN_x=11.0, hN_side=-1, h1_side=-1, hC_side=1):
    P = {}
    def side(up):
        return 1.5 if up else -1.5
    # b2 44..52 ; face X (z>0) = residues with parity such that 51-54 direct pair shares side with 54 (even in b3 up)
    b3_up = lambda i: i % 2 == 0
    b2_up = lambda k: b3_up(int(round(S23)) - k)
    b1_up = lambda i: b2_up(int(round(S12)) - i)
    b4_up = lambda i: i % 2 == 0
    for k in range(44, 53): P[k] = (4.8, 3.3 * (52 - k), side(b2_up(k)))
    for i in range(31, 41): P[i] = (0.0, 3.3 * (52 - S12 + i), side(b1_up(i)))
    for i in range(54, 59): P[i] = (9.6, 3.3 * (52 - S23 + i), side(b3_up(i)))
    for i in range(60, 67): P[i] = (14.4, 3.3 * (52 - S23 + S34 - i), side(b4_up(i)))
    def hel(a, b, p0, u, z, phase0, ref, toward):
        u = np.array([u[0], u[1], 0.0]); u /= np.linalg.norm(u)
        v = np.cross(u, [0, 0, 1.0])            # in-plane perpendicular
        w = np.array([0, 0, 1.0])
        for k in range(a, b + 1):
            ax = np.array([p0[0], p0[1], z]) + 1.5 * (k - a) * u
            ang = np.radians(100.0 * (k - ref) - phase0)     # ang = 0 -> pointing toward the sheet
            P[k] = tuple(ax + 3.3 * (np.cos(ang) * toward * w + np.sin(ang) * v))
    # helix C 70..89 on the X side
    hel(70, 89, hC_start, hC_dir, 10.0 * hC_side, hC_phase, 76, -hC_side)
    # helix 1: residue 26 at the b1 start, running back up to 14
    hel(14, 26, (h1_x, h1_y26 + 1.5 * 12), (0.0, -1.0), 10.0 * h1_side, h1_phase, 14, -h1_side)
    if hN:
        hel(5, 12, (hN_x, 9.0), (0.0, 1.0), 10.0 * hN_side, 20.0, 6, -hN_side)
    return P

def contacts_from(P, cut=8.0, sep=3):
    ks = sorted(P); out = {}
    for a in ks:
        for b in ks:
            if b - a >= sep:
                d = float(np.linalg.norm(np.array(P[a]) - np.array(P[b])))
                if d < cut: out[(a, b)] = d
    return out
