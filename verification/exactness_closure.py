"""
exactness_closure.py -- the FINAL exactness step for Ramanujan Challenge 2026,
Problem 3.1 (7_2 knot, regulator integral = 4 pi^2 / 85).

The by-hand regulator (byhand_flattening.py) reproduces (S_beta - S_alpha)/pi^2 = 4/85
to ~160 digits from our own shapes/flattenings.  160 digits of agreement is NOT a
proof of exact equality.  This script upgrades the numeric agreement to a rigorous
exact statement, resting only on published theorems plus our own rigorous
computation -- no dependence on the private SnapPy Chern-Simons certificate.

The upgrade has two independent halves.

(A) FINITE-LATTICE BOUND  [rigorous field theory + CITED theorems]
    The two endpoint characters have totally-real trace fields: the reciprocal
    fixed fields F_alpha (deg 6) and F_beta (deg 8) of the printed A-polynomial
    factors.  For totally real F, Borel's rank theorem [Borel 1977] makes
    K_3^ind(F) FINITE (rank r_2 = 0), so the Cheeger-Chern-Simons / extended-Rogers
    regulator of each endpoint class is a rational multiple of pi^2.  The
    denominator is controlled by the torsion order w_2(F) = #(roots of unity of
    the relevant twist), bounded elementarily by the degree criterion
        w_2(F) = 2 * prod_p p^{nu_p},   nu_p = max{ nu : phi(p^nu) | 2*deg F }.
    With the 2-torsion change-of-lift between the SL/PSL flattening conventions
    [Zickert 2009; Neumann 2004], the difference (S_beta - S_alpha)/pi^2 (mod 1/2)
    lies on the lattice (1/(2N))Z, N = lcm(2 w_2(F_alpha), 2 w_2(F_beta)).
    We COMPUTE this lattice here (spacing 1/17821440) and CITE Borel/Zickert/Neumann
    for the finiteness + normalization; we do not reprove them.

(B) CERTIFIED INTERVAL  [fully rigorous, our own computation]
    Re-evaluate S_alpha, S_beta on OUR OWN by-hand shapes (the exact real-positive
    algebraic solutions of the by-hand gluing system, reused from byhand_flattening)
    in RIGOROUS INTERVAL / BALL ARITHMETIC (Sage RealBallField / ComplexBallField =
    arb, which carries certified error radii).  This yields a certified enclosure of
    (S_beta - S_alpha)/pi^2 mod 1/2 whose diameter is ~5e-300, i.e. ~1e291 times
    finer than the lattice spacing.  The enclosure is shown to sit strictly inside
    the open lattice cell around 4/85 (both neighbouring lattice points rigorously
    excluded).  Since the true value is a lattice point (A) and the only lattice
    point in the enclosure is 4/85, the difference equals EXACTLY 4/85 mod 1/2.

The regulator's real part (the Chern-Simons value we reduce mod pi^2/2) is, for a
REAL shape z>0, independent of the flattening integers: the p,q terms are purely
imaginary (they carry the volume, which is 0 here).  So Re R(z) is unambiguous and
we compute it with certified bounds via the dilogarithm inversion formula, keeping
every polylog argument strictly OFF the branch cut [1, infty).

Run:  sage -python exactness_closure.py     (HOME=/tmp is set automatically)
"""
import os
import shutil
import sys
from pathlib import Path


def enter_sage_if_needed():
    try:
        import sage.all  # noqa: F401
    except ModuleNotFoundError:
        if os.environ.get("RAMANUJAN31_SAGE_RUNTIME") == "1":
            raise SystemExit("SageMath could not be imported from the selected Sage runtime")
        candidates = [
            os.environ.get("SAGE"),
            shutil.which("sage"),
            "/Applications/SageMath-10-9.app/Contents/Frameworks/"
            "Sage.framework/Versions/Current/local/bin/sage",
            "/usr/local/bin/sage",
        ]
        sage = next((Path(p) for p in candidates if p and Path(p).is_file()), None)
        if sage is None:
            raise SystemExit("SageMath 10.9 is required; set SAGE=/path/to/sage")
        env = os.environ.copy()
        env["RAMANUJAN31_SAGE_RUNTIME"] = "1"
        env.setdefault("SAGE_TMPDIR", "/tmp")
        os.execve(str(sage), [str(sage), "-python", str(Path(__file__).resolve())], env)

    try:
        import snappy  # noqa: F401
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "SnapPy 3.3.2 is missing from Sage's Python. Install it with "
            "`$SAGE -python -m pip install snappy==3.3.2`."
        ) from exc


enter_sage_if_needed()

sys.path.insert(0, str(Path(__file__).resolve().parent))
import byhand_flattening as bh  # noqa: E402  (reuse the by-hand shape machinery)
from sage.all import (PolynomialRing, QQ, ZZ, NumberField, AA,  # noqa: E402
                      RealBallField, ComplexBallField,
                      euler_phi, next_prime, lcm, factor)

PREC = 1000  # working bits for the ball arithmetic (~300 decimal digits)


# ============================================================================
# (A)  FINITE-LATTICE BOUND
# ============================================================================
def reciprocal_transform(P):
    """Palindromic P of degree 2n -> monic g of degree n with g(x + 1/x) = P(x)/x^n.
    g is the minimal polynomial of the totally-real reciprocal fixed field."""
    n = P.degree() // 2
    Ry = PolynomialRing(QQ, "y")
    y = Ry.gen()
    s = [Ry(2), y]                       # s_k = x^k + x^{-k} in terms of y = x+1/x
    for k in range(2, n + 1):
        s.append(y * s[k - 1] - s[k - 2])
    c = P.list()
    g = Ry(c[n])                         # middle coefficient (k = 0 term)
    for k in range(1, n + 1):
        g += Ry(c[n + k]) * s[k]
    return g


def w2_degree_bound(d):
    """Elementary upper bound on the torsion order w_2(F) for [F:Q] = d, via
        nu_p = max{ nu >= 0 : phi(p^nu) | 2 d },   w_2(F) | 2 * prod_p p^{nu_p}.
    Rationale: mu_{p^nu} subset F(mu_{p^nu}) and [F(mu):F] | 2 force
    phi(p^nu) = [Q(mu):Q] | [F(mu):Q] | 2 [F:Q] = 2 d.  Returns (bound, {p:nu}).
    Reproduces the classical w_2(Q) = 24 at d = 1."""
    limit = 2 * d
    fac = {}
    p = 2
    while p - 1 <= limit:                # p can contribute only if phi(p)=p-1 <= 2d
        nu = 0
        while euler_phi(p ** (nu + 1)) <= limit and limit % euler_phi(p ** (nu + 1)) == 0:
            nu += 1
        if nu > 0:
            fac[p] = nu
        p = next_prime(p)
    prod = 1
    for pp, nu in fac.items():
        prod *= pp ** nu
    return 2 * prod, fac


def part_A():
    print("=" * 78)
    print("(A) FINITE-LATTICE BOUND  -- torsion arithmetic of the totally-real")
    print("    reciprocal trace fields   [rigorous field theory + cited theorems]")
    print("=" * 78)

    Rx = PolynomialRing(QQ, "x")
    palindromic = {
        "alpha": Rx(bh.ENDPOINTS["alpha"]["coeffs"]),
        "beta":  Rx(bh.ENDPOINTS["beta"]["coeffs"]),
    }
    orders = {}
    for lbl in ("alpha", "beta"):
        P = palindromic[lbl]
        assert P.list() == P.list()[::-1], f"{lbl} A-poly factor not palindromic"
        g = reciprocal_transform(P)
        F = NumberField(g, "a")
        sig = F.signature()
        disc = F.discriminant()
        d = F.degree()
        assert sig[1] == 0, f"F_{lbl} is NOT totally real: signature {sig}"
        w, fac = w2_degree_bound(d)
        order = 2 * w                     # 2-torsion change-of-lift (Zickert/Neumann)
        orders[lbl] = order
        print(f"\n  F_{lbl}:  g = {g}")
        print(f"    reciprocal (x+1/x) transform of the deg-{P.degree()} A-poly factor")
        print(f"    degree            = {d}")
        print(f"    signature         = {sig}   => TOTALLY REAL (r_2 = 0): "
              f"{sig[1] == 0}")
        print(f"    discriminant      = {disc} = {factor(disc)}")
        print(f"    Borel [1977]: K_3^ind(F_{lbl}) is FINITE (rank r_2 = 0)  [CITED]")
        print("    degree criterion  nu_p = max{{nu : phi(p^nu) | {}}}:".format(2 * d))
        for pp in sorted(fac):
            nu = fac[pp]
            phis = ", ".join(f"phi({pp}^{k})={euler_phi(pp**k)}" for k in range(1, nu + 1))
            print(f"        p={pp:>2}: nu={nu}  ({phis})  -> {pp}^{nu} = {pp**nu}")
        pieces = " * ".join(f"{pp}^{fac[pp]}" for pp in sorted(fac))
        print(f"    w_2(F_{lbl}) | 2 * ({pieces}) = {w}")
        print(f"    order (with 2-torsion change-of-lift) = 2 * {w} = {order}  [Zickert/Neumann]")

    N = lcm(orders["alpha"], orders["beta"])
    spacing = QQ(1) / (2 * N)
    print("\n  " + "-" * 74)
    print(f"  N = lcm(order_alpha, order_beta) = lcm({orders['alpha']}, "
          f"{orders['beta']}) = {N}")
    print(f"    = {factor(N)}")
    print(f"  LATTICE SPACING of (S_beta - S_alpha)/pi^2 mod 1/2  =  1/(2N) = 1/{2*N}")
    assert 2 * N == 17821440, 2 * N
    # 4/85 must itself be a lattice point:  4/85 = k/(2N)  <=>  8N/85 in Z
    k = QQ(8 * N) / 85
    assert k.denominator() == 1, "4/85 is NOT a lattice point!"
    print(f"  4/85 is a lattice point:  4/85 = {ZZ(k)}/(2N)   (8N/85 = {ZZ(k)} in Z): OK")
    return N, spacing


# ============================================================================
# (B)  CERTIFIED INTERVAL
# ============================================================================
def solve_shapes_exact(coeffs, Mfun, Lfun, hint, matrix):
    """Exact real-positive algebraic shapes (elements of AA) of the by-hand gluing
    system -- same system as byhand_flattening.solve_shapes, but returned as exact
    algebraic reals so their ball enclosures are rigorous."""
    Rx = PolynomialRing(QQ, "x")
    poly = Rx(coeffs)
    rts = poly.roots(AA, multiplicities=False)
    hv = float(hint)
    root = min(rts, key=lambda r: abs(float(r) - hv))
    assert abs(float(root) - hv) < 1e-12, "endpoint field root far from printed hint"
    F = NumberField(poly, "a", embedding=root)     # embed into AA at the real root
    g = F.gen()
    M, L = Mfun(g), Lfun(g)
    R = PolynomialRing(F, ["z0", "z1", "z2", "z3"])
    zs = R.gens()
    edge = matrix[:4]
    mer, lon, _ = bh.cusp_rows_public()
    eqs = [bh._laurent_poly(R, F, zs, edge[i], F(1)) for i in range(4)]
    eqs.append(bh._laurent_poly(R, F, zs, mer, M ** -2))
    eqs.append(bh._laurent_poly(R, F, zs, lon, L ** -2))
    I = R.ideal(eqs)
    assert I.dimension() == 0, I.dimension()
    real_pos = []
    for sol in I.variety():
        coords = [AA(sol[zs[j]]) for j in range(4)]
        if all(c > 0 for c in coords):
            real_pos.append(coords)
    assert len(real_pos) == 1, f"expected a unique real-positive shape, got {len(real_pos)}"
    return real_pos[0]


def re_R_ball(z_AA, RBF, CBF, pi2):
    """Certified enclosure of Re R_neumann(z; p, q) for a REAL shape z>0.

    For real z>0 the flattening (p,q) contributes only to the imaginary part
    (the volume, = 0 here), so
        Re R(z) = Re Li_2(z) + (1/2) log(z) log|1-z| - pi^2/6.
    Every dilog argument is kept strictly OFF the cut [1, infty):
      * 0<z<1 : Li_2(z) is real, evaluate directly (arb, certified).
      * z>1   : inversion  Re Li_2(z) = pi^2/3 - (1/2)(log z)^2 - Li_2(1/z),
                and Li_2(1/z) has 1/z in (0,1), off the cut (arb, certified)."""
    z = RBF(z_AA)                                   # rigorous ball enclosing the algebraic z
    log_z = z.log()
    if z_AA > 1:
        re_li2 = pi2 / 3 - (log_z ** 2) / 2 - CBF(1 / z).polylog(2).real()
        log_abs_1mz = (z - 1).log()                 # log|1 - z| = log(z - 1)
    else:
        re_li2 = CBF(z).polylog(2).real()
        log_abs_1mz = (1 - z).log()                 # 1 - z > 0
    return re_li2 + log_z * log_abs_1mz / 2 - pi2 / 6


def reduce_mod_half_ball(ball):
    """Reduce a real ball into the (-1/4, 1/4] representative of (value mod 1/2),
    by subtracting the exact rational (nearest integer)*(1/2).  Rigorous: the
    subtracted constant is exact, so the returned ball still encloses the true value."""
    k = ZZ(int(round(float(ball.mid()) / 0.5)))
    return ball - QQ(k) / 2


def part_B(N, spacing):
    print("\n" + "=" * 78)
    print("(B) CERTIFIED INTERVAL  -- regulator in rigorous ball arithmetic")
    print(f"    on our own by-hand exact shapes   [arb, {PREC} bits ~ "
          f"{int(PREC*0.301)} digits]")
    print("=" * 78)

    RBF = RealBallField(PREC)
    CBF = ComplexBallField(PREC)
    pi2 = RBF.pi() ** 2

    reduced = {}
    for lbl in ("alpha", "beta"):
        spec = bh.ENDPOINTS[lbl]
        matrix, _, _ = bh.build_matrix(spec["slope"])
        shapes = solve_shapes_exact(spec["coeffs"], spec["Mfun"], spec["Lfun"],
                                    spec["hint"], matrix)
        S = RBF(0)
        for z in shapes:
            S = S + re_R_ball(z, RBF, CBF, pi2)
        over = S / pi2
        red = reduce_mod_half_ball(over)
        reduced[lbl] = (over, red)
        num, den = spec["target"]
        tgt = reduce_mod_half_ball(RBF(QQ(num) / den))
        sits = (red - tgt).contains_zero()
        rad = red.rad()
        print(f"\n  endpoint {lbl}  (target S/pi^2 = {num}/{den}):")
        print(f"    Re(S)/pi^2 mod 1/2 (certified) = {red}")
        print(f"    printed rep {num}/{den} mod 1/2      = {tgt}   "
              f"[enclosure brackets it: {bool(sits)}]")
        print(f"    enclosure radius               ~ {RBF(rad)}")

    # ---- the difference, and the exact selection --------------------------------
    print("\n  " + "-" * 74)
    print("  DIFFERENCE  (S_beta - S_alpha)/pi^2  mod 1/2   vs   4/85")
    raw_diff = reduced["beta"][0] - reduced["alpha"][0]
    D = reduce_mod_half_ball(raw_diff)
    q = QQ(4) / 85
    left = q - spacing
    right = q + spacing

    diam = RBF(D.rad()) * 2
    print(f"    certified enclosure D          = {D}")
    print(f"    4/85                           = {RBF(q)}")
    print(f"    enclosure diameter             ~ {diam}")
    print(f"    lattice spacing 1/(2N)         = 1/{2*N} ~ {RBF(spacing)}")
    ratio = RBF(spacing) / diam
    print(f"    spacing / diameter             ~ {ratio}   (>> 1 required)")

    # Rigorous selection certificate:
    #   D strictly inside the OPEN lattice cell (4/85 - s, 4/85 + s):
    #     every point of D exceeds the left neighbour  AND  is below the right one.
    #   That open cell contains exactly one lattice point, 4/85.  The true value is
    #   a lattice point (part A) lying in D, hence equals 4/85.
    above_left = (D - RBF(left)).lower() > 0      # D > 4/85 - spacing  (rigorous)
    below_right = (RBF(right) - D).lower() > 0     # D < 4/85 + spacing  (rigorous)
    contains_q = (D - RBF(q)).contains_zero()      # 4/85 lies inside D
    print(f"    left  neighbour 4/85 - 1/(2N) rigorously EXCLUDED (D above it): "
          f"{bool(above_left)}")
    print(f"    right neighbour 4/85 + 1/(2N) rigorously EXCLUDED (D below it): "
          f"{bool(below_right)}")
    print(f"    4/85 itself lies inside the enclosure D:                        "
          f"{bool(contains_q)}")

    selected = bool(above_left and below_right and contains_q)
    assert selected, "certified interval did NOT isolate 4/85 -- do not claim closure"
    print("\n    => D sits strictly inside the open lattice cell around 4/85, whose")
    print("       only lattice point is 4/85.  The true value is a lattice point (A)")
    print("       inside D, so  (S_beta - S_alpha)/pi^2 = 4/85  EXACTLY (mod 1/2).")
    return selected


def main():
    import sage.all as sg
    print(f"sage={sg.version()}  arb ball precision={PREC} bits")
    N, spacing = part_A()
    ok = part_B(N, spacing)
    print("\n" + "=" * 78)
    print("VERDICT")
    print("=" * 78)
    if ok:
        print("CLOSED-FULLY (modulo cited published theorems).")
        print("  (A) The totally-real reciprocal trace fields F_alpha (deg 6, disc 5^3*3881)")
        print("      and F_beta (deg 8, disc 17^7) give, via the elementary degree criterion")
        print("      + Borel finiteness + the Zickert/Neumann 2-torsion normalization, a")
        print("      lattice of spacing 1/17821440 for (S_beta - S_alpha)/pi^2 mod 1/2.")
        print("  (B) Our own by-hand exact shapes, evaluated in certified arb ball arithmetic")
        print("      with every dilog off the branch cut, enclose the difference to diameter")
        print("      ~5e-300 and select 4/85 as the unique lattice point (both neighbours")
        print("      rigorously excluded).")
        print("  Rests only on published theorems (Borel 1977 finiteness; Zickert 2009 /")
        print("  Neumann 2004 extended-Bloch regulator normalization) plus our own exact")
        print("  field theory and interval computation.  NO private SnapPy CS certificate.")
        print("  Honest caveats: the degree criterion is an UPPER bound on w_2 (=> a")
        print("  possibly-finer, still valid lattice); the exact power of 2 in the spacing")
        print("  follows the cited SL/PSL normalization and is not re-derived here.  Neither")
        print("  matters for the conclusion: the enclosure is ~1e291x finer than the spacing,")
        print("  so 4/85 is selected under any spacing >= 1/(2N).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
