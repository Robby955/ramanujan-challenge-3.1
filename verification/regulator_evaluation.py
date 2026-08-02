"""
Independent evaluation attempt of the extended Rogers / Neumann regulator S at the
two endpoints of solution.tex (Problem 3.1), targets

    S_alpha = -26 pi^2 / 15      (= -1.73333... pi^2)   mod pi^2/2
    S_beta  = -86 pi^2 / 51      (= -1.68627... pi^2)   mod pi^2/2
    S_beta - S_alpha = 4 pi^2/85                        mod pi^2/2.

WHAT THIS SCRIPT CAN AND CANNOT DO -- read before trusting any "PASS".
------------------------------------------------------------------------
Neumann's extended Rogers sum on the extended Bloch group is
        S = sum_j  R( z_j ; p_j , q_j )   ( + a constant 6-torsion correction ),
        R(z;p,q) = Li_2(z) + (1/2) log z log(1-z)
                   + (pi i /2)( p log(1-z) + q log z ) - pi^2/6 ,
where (z_0..z_3) are the four tetrahedron cross-ratios and (p_j,q_j) are INTEGER
log-branch / flattening data.  To evaluate S one needs BOTH:
   (i)  the four explicit signed cross-ratios z_j at each endpoint, and
   (ii) the integer flattenings (p_j,q_j), i.e. the strong flattening that solves the
        edge- and Dehn-filling relations of the 5x12 PGL gluing matrix.

From the printed paper (solution.tex) + the Dehn-wedge reconstruction we have ONLY:
   - M,L,b,c,e at each endpoint (exact, high precision)  -> yes
   - the three opposite-edge PRODUCTS per tetrahedron as SIGN-FREE monomials  -> yes
   - the explicit 24 SIGNED coordinates                                  -> NO (unprinted)
   - the ordered (z_0,z_1,z_2,z_3) primary-shape assignment              -> NO (unprinted)
   - the 5x12 PGL gluing matrix that fixes (p_j,q_j)                     -> NO (unprinted)

So this script does the honest maximum: it reconstructs the three equation-terms per
tetrahedron (which DO sum to zero, the Plucker relation), forms every cross-ratio the
printed data leaves open, evaluates the PRINCIPAL-branch (p=q=0) Neumann sum for each,
and reports whether ANY branch-0 assignment reproduces the target -- WITHOUT claiming
that a match constitutes an independent certification (it cannot, because the branch
integers are exactly the unprinted data).
"""
import itertools
import mpmath as mp

mp.mp.dps = 80
PI = mp.pi
PI2 = PI**2

# ---------------------------------------------------------------- endpoints
def fa(t):
    return (t**12-3*t**11+4*t**10-5*t**9+6*t**8-7*t**7+7*t**6-7*t**5
            +6*t**4-5*t**3+4*t**2-3*t+1)
def fb(r):
    return (r**16-7*r**15+22*r**14-48*r**13+87*r**12-133*r**11+178*r**10-211*r**9
            +223*r**8-211*r**7+178*r**6-133*r**5+87*r**4-48*r**3+22*r**2-7*r+1)

t = mp.findroot(fa, mp.mpf('0.59098942867025644049'))
r = mp.findroot(fb, mp.mpf('0.40681308133678976238'))

def bce(Mv, Lv):
    b2 = (1-Mv**2)*(Lv-Mv**4)/(Lv*Mv**2*(Mv**2+Lv))
    bv = -mp.sqrt(b2)                     # negative continuous branch (paper)
    ev = -(Lv*b2+1)/Mv
    cv = bv*(Lv*b2+1-Mv**2)
    return bv, cv, ev

# per-tetrahedron three signed terms of the printed Ptolemy equations (eq:ptolemy-system)
# each equation's three terms sum to 0 : this is the Plucker relation of that tetrahedron.
def tet_terms(Mv, Lv):
    mv = 1/Mv
    bv, cv, ev = bce(Mv, Lv)
    b2 = bv*bv
    T = [
        [-Mv*bv,            -mv*cv,               -bv*ev            ],   # tet0
        [-Lv*mv*b2,         -mv,                  -ev               ],   # tet1
        [ Lv*mv**4*bv,       Lv*mv*bv*ev,         -cv               ],   # tet2
        [ Mv*cv**2,         -Mv*ev**2,             cv               ],   # tet3
    ]
    return T, (bv, cv, ev)

# ---------------------------------------------------------------- Neumann principal-branch R
def R_principal(z):
    """Neumann extended Rogers with p=q=0 (principal logs): Li2(z)+1/2 log z log(1-z) - pi^2/6."""
    return mp.polylog(2, z) + mp.mpf('0.5')*mp.log(z)*mp.log(1-z) - PI2/6

def cross_ratios_from_terms(terms):
    """
    terms sum to 0: T0+T1+T2+T3? no -- three terms, T0+T1+T2 = 0.
    A cross ratio uses one term as denominator: z = -T_a / T_c , 1-z = -T_b / T_c,
    with {a,b,c} a permutation of {0,1,2}. Because T_a+T_b = -T_c, z+(1-z)=1 holds.
    Return the set of distinct z-values reachable (the 6 anharmonic images).
    """
    zs = []
    for c in range(3):
        rest = [i for i in range(3) if i != c]
        for a, b in (rest, rest[::-1]):
            z = -terms[a]/terms[c]
            zs.append(z)
    return zs

def evaluate(name, Mv, Lv):
    print(f"\n================  {name}   (M,L)=({mp.nstr(Mv,12)},{mp.nstr(Lv,12)})")
    T, (bv, cv, ev) = tet_terms(Mv, Lv)
    print(f"   b={mp.nstr(bv,10)}  c={mp.nstr(cv,10)}  e={mp.nstr(ev,10)}   (all <0 as required)")
    # verify each equation's three terms sum to zero (Plucker relation holds numerically)
    for j in range(4):
        s = sum(T[j])
        print(f"   tet{j}: terms sum = {mp.nstr(s,4)}   {'OK (=0)' if abs(s)<mp.mpf(10)**-40 else 'NONZERO!'}")
    # for each tetrahedron collect the candidate cross ratios (6 anharmonic images each)
    cand = [cross_ratios_from_terms(T[j]) for j in range(4)]
    print("   candidate primary cross-ratios per tet (6 anharmonic images):")
    for j in range(4):
        uniq = []
        for z in cand[j]:
            if not any(abs(z-u) < mp.mpf(10)**-30 for u in uniq):
                uniq.append(z)
        print(f"     tet{j}: " + ", ".join(mp.nstr(z,8) for z in uniq))
    return T, cand

TA, candA = evaluate("ALPHA  (-1,2) filling", t**2, t)
TB, candB = evaluate("BETA   (-1,1) filling", r, r)

# ---------------------------------------------------------------- constrained scan, branch p=q=0
def scan(name, cand, target_frac):
    """
    The endpoint representations are REAL (totally real fields), so the geometric
    volume is 0.  In the raw Rogers sum Ssum = Vol + i CS this means Re(Ssum) = 0, and
    the paper's i(Vol+iCS) convention gives S = i*Ssum, so Re(S)/pi^2 = -Im(Ssum)/pi^2
    is the number to compare with the target.  We keep only principal-branch (p=q=0)
    assignments whose reconstructed volume Re(Ssum) vanishes, and report the smallest
    volume residual actually achievable on the principal branch.
    """
    print(f"\n----  {name}: principal-branch (p=q=0), assignments with reconstructed Vol=0")
    tgt = mp.mpf(target_frac[0])/mp.mpf(target_frac[1])
    voltol = mp.mpf(10)**-30
    hits = {}          # rounded CS/pi^2 mod 1/2  ->  count
    target_reached = False
    min_imres = None   # smallest |Im(Ssum)| = smallest reconstructed "volume residual"
    for combo in itertools.product(range(6), repeat=4):
        Ssum = mp.mpc(0)
        ok = True
        for j in range(4):
            z = cand[j][combo[j]]
            if abs(z) < mp.mpf(10)**-20 or abs(1-z) < mp.mpf(10)**-20:
                ok = False; break
            Ssum += R_principal(z)
        if not ok:
            continue
        volres = abs(mp.re(Ssum))/PI2         # reconstructed volume residual (must be 0)
        if min_imres is None or volres < min_imres:
            min_imres = volres
        if abs(mp.re(Ssum)) > voltol:         # volume must vanish for a real rep
            continue
        cs = mp.im(Ssum)/PI2                   # CS/pi^2
        s_over_pi2 = -cs                       # i(Vol+iCS): Re(S)/pi^2 = -Im(Ssum)/pi^2
        # reduce mod 1/2 into [0,1/2)
        red = s_over_pi2 - mp.floor(s_over_pi2/mp.mpf('0.5'))*mp.mpf('0.5')
        key = mp.nstr(red, 10)
        hits[key] = hits.get(key, 0) + 1
        d = s_over_pi2 - tgt
        d = d - mp.nint(d/mp.mpf('0.5'))*mp.mpf('0.5')
        if abs(d) < mp.mpf(10)**-10:
            target_reached = True
    tgt_red = tgt - mp.floor(tgt/mp.mpf('0.5'))*mp.mpf('0.5')
    print(f"   target Re(S)/pi^2 = {target_frac[0]}/{target_frac[1]} = {mp.nstr(tgt,10)}   (mod 1/2 -> {mp.nstr(tgt_red,10)})")
    print(f"   smallest volume residual |Re(Ssum)|/pi^2 over ALL principal-branch assignments: {mp.nstr(min_imres,6)}")
    print(f"   # principal-branch assignments with Vol=0: {sum(hits.values())}")
    print(f"   distinct S/pi^2 values reached (mod 1/2), with multiplicity:")
    for k in sorted(hits, key=lambda x: float(x)):
        mark = "  <-- TARGET" if abs(mp.mpf(k)-tgt_red) < mp.mpf(10)**-8 else ""
        print(f"       {k}   (x{hits[k]}){mark}")
    print(f"   -> target {'IS' if target_reached else 'is NOT'} among the reachable principal-branch values")
    return target_reached

bA = scan("ALPHA", candA, (-26, 15))
bB = scan("BETA",  candB, (-86, 51))

print("\n================  DIFFERENCE  S_beta - S_alpha vs 4/85")
print(f"   4/85 = {mp.nstr(mp.mpf(4)/85,12)}")
print("   (a principal-branch difference is only meaningful if BOTH endpoints matched;")
print("    otherwise the branch integers are exactly the missing data.)")
