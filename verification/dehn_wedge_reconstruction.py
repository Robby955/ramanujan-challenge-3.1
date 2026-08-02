"""
Independent reconstruction of the Dehn-invariant "sign check" in solution.tex sec 5:

  "the certificate forms z_i wedge (1-z_i) from the exponent vectors of all 24
   signed-coordinate monomials.  Every wedge involving b,c,e cancels, and the sole
   surviving coefficient is +2 on M wedge L."   ==>  nu = 2 M wedge L    (eq:dehn)

WHAT IS AND IS NOT NEEDED FROM THE (unprinted) TRIANGULATION
------------------------------------------------------------
For an N=2 (SL2) Ptolemy variety, each tetrahedron carries 6 edge coordinates
c_1100,c_0011,c_1010,c_0101,c_1001,c_0110 and ONE Ptolemy relation, the 3-term
Pluecker identity
        P1 - P2 + P3 = 0,   P1=c_1100 c_0011,  P2=c_1010 c_0101,  P3=c_1001 c_0110,
where P1,P2,P3 are the three products of OPPOSITE edge pairs.  The four printed
equations (solution.tex eq:ptolemy-system) each have exactly three monomial terms;
those three terms ARE the three opposite-edge products P1,P2,P3 of that tetrahedron
(after face-identifications express the 24 coords as signed monomials in M,L,b,c,e).

The Garoufalidis-Thurston-Zickert cross-ratio (public formula):
        z   = - P1/P2 ,       1 - z = - P1/P3 .
So the wedge contribution of one tetrahedron, mod 2-torsion (d log(-1)=0), is
        {z} wedge {1-z} = (u1-u2) wedge (u1-u3),   u_k = log(monomial P_k).
Algebra:  (u1-u2)wedge(u1-u3) = u1^u2 + u2^u3 + u3^u1  (cyclic, fully symmetric under
3-cycles, sign flip under a transposition = choice of orientation z<->1-z).

Two robustness facts (proved inline below) make the reconstruction independent of the
NOT-printed vertex ordering / edge labeling:
  (A) common monomial rescaling of a whole equation (clearing c_1100,0=1 etc.) leaves
      S_i = u1^u2+u2^u3+u3^u1 unchanged;
  (B) which of the three terms is P1/P2/P3 only changes S_i by an overall sign
      (tetrahedron orientation eps_i in {+1,-1}).
Hence nu = sum_i eps_i S_i, and the ONLY datum not fixed by the printed equations is
the 4 orientation signs eps_i.  We solve for them by demanding the b,c,e parts cancel
and read off the M^L coefficient.  A clean +/-2 M^L with a unique sign pattern (up to
global orientation) is a genuine reproduction of the claim.
"""
import itertools

GENS = ['M', 'L', 'b', 'c', 'e']
IDX = {g: i for i, g in enumerate(GENS)}

def vec(**kw):
    v = [0]*5
    for g, p in kw.items():
        v[IDX[g]] += p
    return tuple(v)

# The three monomial terms (as log-exponent vectors) of each printed Ptolemy relation.
# m = 1/M, so log m = -M.
# tet0:  -M b - m c - b e = 0     -> terms  M b ,  M^-1 c ,  b e
# tet1:  -L m b^2 - m - e = 0     -> terms  L M^-1 b^2 ,  M^-1 ,  e
# tet2:   L m^4 b + L m b e - c=0 -> terms  L M^-4 b ,  L M^-1 b e ,  c
# tet3:   M c^2 - M e^2 + c = 0   -> terms  M c^2 ,  M e^2 ,  c
TETS = {
    0: [vec(M=1, b=1),          vec(M=-1, c=1),           vec(b=1, e=1)],
    1: [vec(L=1, M=-1, b=2),    vec(M=-1),                vec(e=1)],
    2: [vec(L=1, M=-4, b=1),    vec(L=1, M=-1, b=1, e=1), vec(c=1)],
    3: [vec(M=1, c=2),          vec(M=1, e=2),            vec(c=1)],
}

# ---------- exterior algebra on 5 generators (coeffs over Z), keyed by (i<j) ----------
def wedge(u, w):
    """u wedge w in ^2(Z^5); returns dict {(i,j): coeff} with i<j."""
    out = {}
    for i in range(5):
        if u[i] == 0:
            continue
        for j in range(5):
            if w[j] == 0 or i == j:
                continue
            a, b = (i, j) if i < j else (j, i)
            s = 1 if i < j else -1
            out[(a, b)] = out.get((a, b), 0) + s*u[i]*w[j]
    return {k: v for k, v in out.items() if v}

def add(d1, d2, scale=1):
    out = dict(d1)
    for k, v in d2.items():
        out[k] = out.get(k, 0) + scale*v
    return {k: v for k, v in out.items() if v}

def S_of_tet(terms):
    """S = u1^u2 + u2^u3 + u3^u1 (cyclic symmetric)."""
    u1, u2, u3 = terms
    return add(add(wedge(u1, u2), wedge(u2, u3)), wedge(u3, u1))

def show(d):
    if not d:
        return "0"
    parts = []
    for (i, j), v in sorted(d.items()):
        sign = '+' if v > 0 else '-'
        mag = abs(v)
        parts.append(f"{sign}{mag}*{GENS[i]}^{GENS[j]}")
    return " ".join(parts)

# ---- sanity: (A) common-shift invariance and (B) 3-cycle symmetry, on tet0 ----
def check_invariances():
    t = TETS[0]
    base = S_of_tet(t)
    shift = vec(M=2, L=1, c=3)                     # arbitrary common monomial factor
    tshift = [tuple(a+b for a, b in zip(term, shift)) for term in t]
    same_shift = (S_of_tet(tshift) == base)
    # 3-cycle: rotate the three terms
    rot = [t[1], t[2], t[0]]
    same_rot = (S_of_tet(rot) == base)
    # transposition flips sign
    swp = [t[1], t[0], t[2]]
    flip = (S_of_tet(swp) == {k: -v for k, v in base.items()})
    return same_shift, same_rot, flip

sh, rot, flip = check_invariances()
print("Invariance sanity (make reconstruction independent of unprinted labeling):")
print(f"  (A) common monomial rescaling leaves S_i unchanged : {sh}")
print(f"  (B) 3-cycle of terms leaves S_i unchanged          : {rot}")
print(f"  (B) transposition of terms flips sign of S_i       : {flip}")
print()

S = {i: S_of_tet(TETS[i]) for i in range(4)}
print("Per-tetrahedron S_i = z_i wedge (1-z_i)  (from opposite-edge products):")
for i in range(4):
    print(f"  S_{i} = {show(S[i])}")
print()

ML = (IDX['M'], IDX['L'])
bce_keys = [(i, j) for i in range(5) for j in range(i+1, 5)
            if GENS[i] in 'bce' or GENS[j] in 'bce']

print("Searching orientation signs eps in {+1,-1}^4 that cancel every b,c,e wedge:")
solutions = []
for eps in itertools.product([1, -1], repeat=4):
    tot = {}
    for i in range(4):
        tot = add(tot, S[i], scale=eps[i])
    if all(tot.get(k, 0) == 0 for k in bce_keys):    # no b,c,e survivors
        solutions.append((eps, tot))

for eps, tot in solutions:
    ml = tot.get(ML, 0)
    others = {k: v for k, v in tot.items() if k != ML}
    print(f"  eps={eps}:  total = {show(tot)}   |  M^L coeff = {ml:+d}"
          f"   (other terms: {show(others) if others else 'none'})")
print()

# verdict
clean = [(eps, tot) for eps, tot in solutions
         if all(k == ML for k in tot) and abs(tot.get(ML, 0)) == 2]
if clean:
    eps, tot = clean[0]
    coeff = tot[ML]
    print(f"RESULT: with orientation signs eps={eps} (and its global-sign twin),")
    print(f"        nu = sum eps_i (z_i wedge (1-z_i)) = {coeff:+d} * M^L, ")
    print(f"        every b,c,e wedge cancels.  Matches solution.tex eq:dehn nu = 2 M^L.")
    print("VERDICT: wedge sign-check REPRODUCED (sign fixed by orientation convention).")
else:
    print("VERDICT: no sign pattern yields a clean +/-2 M^L; claim NOT reproduced.")
