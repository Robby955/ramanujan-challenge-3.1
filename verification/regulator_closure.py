"""
regulator_closure.py -- independent reproduction of the Section-4 Chern-Simons
regulator values S_alpha, S_beta for Ramanujan Challenge 2026, Problem 3.1
(the 7_2 regulator integral = 4 pi^2 / 85).

TARGETS (mod pi^2/2):
    S_alpha         = -26 pi^2 / 15      (= -1.7333333... pi^2)
    S_beta          = -86 pi^2 / 51      (= -1.6862745... pi^2)
    S_beta - S_alpha=   4 pi^2 / 85

WHAT MAKES THIS INDEPENDENT OF THE CERTIFICATE'S OWN CHERN-SIMONS COMPUTATION
----------------------------------------------------------------------------
The blocker recorded in independent_verification/LEDGER.md (Sec 4, CERTIFICATE-
BOUND) was that evaluating S needs two ingredients the paper delegates to SnapPy
and never prints: the integer flattenings (p_j, q_j) and the 5x12 PGL gluing
matrix of the four-tetrahedron 7_2 triangulation.

This script:
  1. Solves the enhanced-Ptolemy system (b,c,e) exactly over each endpoint field
     (our own Groebner solve; it equals the printed eq:generic-ptolemy formulas).
  2. Reads the two delegated ingredients from SnapPy's shipped extended-Ptolemy
     interface: the flattening triples (z_j, p_j, q_j) and the 5x12 PGL gluing
     matrix. This is SnapPy's public triangulation data, re-runnable by anyone
     with SnapPy 3.3.2 -- NOT a private oracle for the answer.
  3. INDEPENDENT BRANCH CHECK: verifies the lifted flattening logs satisfy every
     edge row and the filling row of the 5x12 matrix (residual ~ 0). This proves
     (p_j, q_j) are the correct branch integers WITHOUT calling SnapPy's
     complex_volume / verified_complex_volume_from_lifted_ptolemys.
  4. OUR EVALUATOR: evaluates the extended-Rogers sum S = sum_j R(z_j;p_j,q_j) in
     mpmath (Neumann's function), a transcendental-dilog implementation entirely
     separate from pari's dilog and SnapPy's _L_function. SnapPy's own
     complex_volume is computed only as a labelled cross-check, never as the
     result.

So the topological data (which shape is z, the branch integers, the gluing
matrix) comes from SnapPy; the regulator arithmetic is ours. Verdict target:
CLOSED-VIA-SNAPPY.
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

import mpmath as mp
import snappy
from sage.all import PolynomialRing, QQ, NumberField, RealField
from sage.env import SAGE_VERSION
from snappy.dev.extended_ptolemy import extended
from snappy.ptolemy.coordinates import PtolemyCoordinates
from snappy import pari

assert snappy.__version__ == "3.3.2", snappy.__version__
assert SAGE_VERSION == "10.9", SAGE_VERSION

PREC_DECIMALS = 180          # pari / RealField working precision
pari.set_real_precision(PREC_DECIMALS)
RF = RealField(600)
mp.mp.dps = 150              # our mpmath evaluator precision

PI = mp.pi
PI2 = PI ** 2
PIi = mp.mpc(0, 1) * PI

ENDPOINTS = {
    "alpha": dict(
        slope=(-1, 2),
        coeffs=(1, -3, 4, -5, 6, -7, 7, -7, 6, -5, 4, -3, 1),
        hint="0.59098942867025644049",
        M=lambda t: t ** 2, L=lambda t: t,
        target=(-26, 15),
    ),
    "beta": dict(
        slope=(-1, 1),
        coeffs=(1, -7, 22, -48, 87, -133, 178, -211, 223, -211, 178,
                -133, 87, -48, 22, -7, 1),
        hint="0.40681308133678976238",
        M=lambda r: r, L=lambda r: r,
        target=(-86, 51),
    ),
}


def solve_bce(field, M, L):
    """Our own exact Groebner solve of the enhanced-Ptolemy system in b,c,e
    (eq:ptolemy-system).  Independent of the certificate."""
    m = 1 / M
    R = PolynomialRing(field, names=("b", "c", "e"))
    b, c, e = R.gens()
    eqs = [-M * b - m * c - b * e,
           -L * m * b ** 2 - m - e,
           L * m ** 4 * b + L * m * b * e - c,
           M * c ** 2 - M * e ** 2 + c]
    G = R.ideal(eqs).groebner_basis()
    vals = {}
    for var, name in [(b, "c_1010_0"), (c, "c_1001_0"), (e, "c_0101_0")]:
        matches = [g for g in G
                   if g.monomial_coefficient(var) == 1 and g - var in field]
        assert len(matches) == 1
        vals[name] = field(-(matches[0] - var))
    sub = {b: vals["c_1010_0"], c: vals["c_1001_0"], e: vals["c_0101_0"]}
    assert [q.subs(sub) for q in eqs] == [0, 0, 0, 0]
    # cross-check against the printed closed form eq:generic-ptolemy
    b2 = (1 - M ** 2) * (L - M ** 4) / (L * M ** 2 * (M ** 2 + L))
    assert vals["c_1010_0"] ** 2 == b2
    assert vals["c_0101_0"] == -(L * b2 + 1) / M
    assert vals["c_1001_0"] == vals["c_1010_0"] * (L * b2 + 1 - M ** 2)
    vals.update({"c_1100_0": field(1), "M": M, "L": L, "m": m, "l": 1 / L})
    return vals


# ---------- OUR extended-Rogers evaluator (Neumann), pure mpmath ----------
def R_neumann(z, p, q):
    """Neumann's extended Rogers function
        R(z;p,q) = Li_2(z) + (1/2) log z log(1-z)
                   + (pi i/2)(p log(1-z) + q log z) - pi^2/6.
    Principal branches (mpmath), independent of pari/SnapPy."""
    z = mp.mpc(z)
    lz = mp.log(z)
    l1z = mp.log(1 - z)
    return (mp.polylog(2, z) + mp.mpf('0.5') * lz * l1z
            + (PIi / 2) * (p * l1z + q * lz) - PI2 / 6)


def L_snappy_convention(z, p, q):
    """SnapPy's _L_function convention (for reconciliation only):
        L(z;p,q) = Li_2(z) + (log z + p pi i)(log(1-z) + q pi i)/2 - pi^2/6.
    Differs from R_neumann by -p q pi^2/2 (integer multiple of pi^2/2)."""
    z = mp.mpc(z)
    lz = mp.log(z)
    l1z = mp.log(1 - z)
    return (mp.polylog(2, z) + (lz + p * PIi) * (l1z + q * PIi) / 2 - PI2 / 6)


def reduce_mod_half(x):
    """Reduce a real number into (-1/4, 1/4] representative of x mod 1/2."""
    x = mp.mpf(x)
    return x - mp.nint(x / mp.mpf('0.5')) * mp.mpf('0.5')


def run_endpoint(label):
    spec = ENDPOINTS[label]
    print("=" * 78)
    print(f"ENDPOINT {label}   filling {spec['slope']}   target "
          f"S/pi^2 = {spec['target'][0]}/{spec['target'][1]}")

    # 1. exact (b,c,e)
    Rx = PolynomialRing(QQ, "x")
    poly = Rx(spec["coeffs"])
    F = NumberField(poly, label[0])
    g = F.gen()
    M = spec["M"](g)
    L = spec["L"](g)
    vals = solve_bce(F, M, L)
    print(f"  enhanced-Ptolemy (b,c,e) solved exactly over deg-{F.degree()} field "
          f"and matched eq:generic-ptolemy")

    # 2. signed coords + flattenings + gluing matrix from SnapPy
    man = snappy.ManifoldHP("7_2")
    man.dehn_fill(spec["slope"])
    ideal, full_signed = extended.ptolemy_ideal_for_filled(
        man, nonzero_cond=False, return_full_var_dict=True, notation="full")
    ring_sub = {ideal.ring()(k): v for k, v in vals.items()}
    assert all(gen.subs(ring_sub) == 0 for gen in ideal.gens())
    signed = {k: F(v.subs(ring_sub))
              for k, v in full_signed.items() if k.startswith("c_")}
    assert len(signed) == 24

    root = [rt for rt in poly.roots(RF, multiplicities=False)
            if abs(rt - RF(spec["hint"])) < RF("1e-20")][0]
    emb = F.hom([root], RF, check=False)
    num_signed = {k: emb(v) for k, v in signed.items()}
    P = PtolemyCoordinates(num_signed, is_numerical=True,
                           manifold_thunk=lambda: man)
    cr = P.cross_ratios()
    cr.check_against_manifold(man, epsilon=1e-100)
    flat = P.flattenings_numerical()

    ntet = man.num_tetrahedra()
    triples = []
    print("  extracted flattening triples (z_j, p_j, q_j)  [SnapPy topological data]:")
    for j in range(ntet):
        z, p, q = flat.get_zpq_triple(f"z_0000_{j}")
        triples.append((mp.mpf(str(pari(z))), int(p), int(q)))
        print(f"    tet{j}:  p={p}  q={q}   z={mp.nstr(mp.mpf(str(pari(z))), 22)}")

    eqns = man.gluing_equations_pgl(2, equation_type="all")
    cols = eqns.explain_columns
    rows = eqns.explain_rows
    print(f"  5x{len(cols)} PGL gluing matrix rows={rows}:")
    for r in range(len(rows)):
        print(f"    {rows[r]:>13}: "
              f"{[int(eqns.matrix[r, c]) for c in range(len(cols))]}")

    # 3. INDEPENDENT branch check: flattening logs satisfy every gluing row
    #    (uses the stored lifted logs w0,w1,w2 and the public gluing matrix;
    #     no call to complex_volume).
    wlog = {}
    for j in range(ntet):
        wlog[f"z_0000_{j}"] = mp.mpc(str(pari(flat[f"z_0000_{j}"][0]).real()),
                                     str(pari(flat[f"z_0000_{j}"][0]).imag()))
        wlog[f"zp_0000_{j}"] = mp.mpc(str(pari(flat[f"zp_0000_{j}"][0]).real()),
                                      str(pari(flat[f"zp_0000_{j}"][0]).imag()))
        wlog[f"zpp_0000_{j}"] = mp.mpc(str(pari(flat[f"zpp_0000_{j}"][0]).real()),
                                       str(pari(flat[f"zpp_0000_{j}"][0]).imag()))
    max_resid = mp.mpf(0)
    for r in range(len(rows)):
        resid = sum(int(eqns.matrix[r, c]) * wlog[cols[c]]
                    for c in range(len(cols)))
        max_resid = max(max_resid, abs(resid))
    print(f"  INDEPENDENT gluing-equation residual on flattening logs "
          f"(all 4 edges + filling): max |row| = {mp.nstr(max_resid, 4)}")
    assert max_resid < mp.mpf(10) ** -60, "flattening does not satisfy gluing eqns!"

    # 4. OUR evaluator: extended-Rogers sum
    S_ours = sum(R_neumann(z, p, q) for (z, p, q) in triples)
    S_snappy_conv = sum(L_snappy_convention(z, p, q) for (z, p, q) in triples)

    return dict(label=label, spec=spec, triples=triples,
                S_ours=S_ours, S_snappy_conv=S_snappy_conv, flat=flat, man=man)


def main():
    print(f"snappy={snappy.__version__}  sage={SAGE_VERSION}  "
          f"mpmath_dps={mp.mp.dps}")
    results = {lbl: run_endpoint(lbl) for lbl in ("alpha", "beta")}

    print("\n" + "=" * 78)
    print("EXTENDED-ROGERS SUM S  (OUR mpmath evaluator), reduced mod pi^2/2")
    print("=" * 78)
    S = {}
    for lbl in ("alpha", "beta"):
        res = results[lbl]
        num, den = res["spec"]["target"]
        s_over_pi2 = mp.re(res["S_ours"]) / PI2      # real reps: Vol=0
        s_full = res["S_ours"] / PI2
        S[lbl] = s_over_pi2
        tgt = mp.mpf(num) / den
        red_ours = reduce_mod_half(s_over_pi2)
        red_tgt = reduce_mod_half(tgt)
        diff = reduce_mod_half(s_over_pi2 - tgt)
        agree = int(-mp.log10(abs(diff))) if diff != 0 else mp.mp.dps
        # cross-check: snappy-convention value reduces to the same thing mod 1/2
        red_snappy = reduce_mod_half(mp.re(res["S_snappy_conv"]) / PI2)
        print(f"\n  {lbl}:")
        print(f"    S/pi^2 (raw, ours)          = {mp.nstr(s_full, 30)}")
        print(f"    S/pi^2 mod 1/2 (ours)       = {mp.nstr(red_ours, 30)}")
        print(f"    target {num}/{den} mod 1/2      = {mp.nstr(red_tgt, 30)}")
        print(f"    |difference| mod 1/2        = {mp.nstr(abs(diff), 6)}  "
              f"-> ~{agree} digits agreement")
        print(f"    (snappy-conv mod 1/2, xchk) = {mp.nstr(red_snappy, 30)}  "
              f"[matches ours: {abs(red_ours-red_snappy) < mp.mpf(10)**-60}]")

    print("\n" + "=" * 78)
    print("DIFFERENCE  S_beta - S_alpha  vs  4/85   (mod 1/2)")
    print("=" * 78)
    d = S["beta"] - S["alpha"]
    tgt = mp.mpf(4) / 85
    dred = reduce_mod_half(d - tgt)
    agree = int(-mp.log10(abs(dred))) if dred != 0 else mp.mp.dps
    print(f"  (S_beta - S_alpha)/pi^2 mod 1/2 = {mp.nstr(reduce_mod_half(d), 30)}")
    print(f"  4/85 mod 1/2                    = {mp.nstr(reduce_mod_half(tgt), 30)}")
    print(f"  |difference| mod 1/2            = {mp.nstr(abs(dred), 6)}  "
          f"-> ~{agree} digits agreement")

    # labelled cross-check: OUR mpmath sum vs SnapPy's own pari _L_function on the
    # SAME triples (internal consistency of the two dilog implementations).
    from snappy.ptolemy.coordinates import _L_function
    print("\n" + "=" * 78)
    print("CROSS-CHECK ONLY: our mpmath R-sum vs SnapPy pari _L_function (same triples)")
    print("=" * 78)
    pari_pi2 = pari("Pi")**2
    for lbl in ("alpha", "beta"):
        res = results[lbl]
        S_snappy_pari = sum(_L_function(t) for t in
                            [res["flat"].get_zpq_triple(f"z_0000_{j}")
                             for j in range(res["man"].num_tetrahedra())])
        snappy_over_pi2 = S_snappy_pari / pari_pi2
        ours_over_pi2 = res["S_snappy_conv"] / PI2   # same convention as _L_function
        re_gap = abs(mp.mpf(str(snappy_over_pi2.real())) - mp.re(ours_over_pi2))
        im_gap = abs(mp.mpf(str(snappy_over_pi2.imag())) - mp.im(ours_over_pi2))
        print(f"  {lbl}: SnapPy pari sum_L/pi^2  = {snappy_over_pi2}")
        print(f"        our mpmath sum_L/pi^2   = {mp.nstr(ours_over_pi2, 30)}")
        print(f"        |re gap|={mp.nstr(re_gap,4)}  |im gap|={mp.nstr(im_gap,4)}"
              f"  -> implementations agree: {re_gap < mp.mpf(10)**-60 and im_gap < mp.mpf(10)**-60}")

    print("\nDONE.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
