"""
ANCHOR CHECK for Ramanujan Challenge 2026, Problem 3.1.

Direct numerical verification of the challenge identity, independent of the
Ptolemy / K-theory / Dehn machinery in solution.tex.  We integrate

    eta = log x * d(log y) - log y * d(log x)

along the real branch y = y(x) of A_{7_2}(x,y) = 0 from x = alpha to x = beta,
and confirm it equals 4*pi^2/85.

Coordinate dictionary (matches solution.tex Sec. "The curve and its two
representations"):
    x = M  (meridian eigenvalue),   y = L  (longitude eigenvalue).
The challenge branch runs from (M,L) = (alpha, sqrt(alpha)) to (beta, beta),
with  alpha = t^2 , y(alpha) = t ,  and  beta = r , y(beta) = r ,
where t is the real root of f_alpha near 0.59098942867 and r the real root of
f_beta near 0.40681308133.

y(x) is obtained by solving A_{7_2}(x,y)=0 for the correct real root at each x
(continuity-tracked); dy/dx by implicit differentiation dy/dx = -A_x / A_y.

The A-polynomial is transcribed exactly from
ptolemy_independent_check.py / solution.tex.
"""

import sympy as sp
import mpmath as mp

mp.mp.dps = 60  # >= 40 digits requested; use 60 for margin

# --------------------------------------------------------------------------
# A_{7_2}(M, L) transcribed exactly (M = x, L = y).
# --------------------------------------------------------------------------
M, L = sp.symbols('M L')
A72 = (L**5
 + L**4*(M**14 - M**12 + 3*M**4 + 4*M**2 - 2)
 + L**3*(-2*M**18 + 5*M**16 + M**14 - 4*M**12 + 6*M**8 + 5*M**6 + 2*M**4 - 4*M**2 + 1)
 + L**2*(M**22 - 4*M**20 + 2*M**18 + 5*M**16 + 6*M**14 - 4*M**10 + M**8 + 5*M**6 - 2*M**4)
 + L*(-2*M**22 + 4*M**20 + 3*M**18 - M**10 + M**8) + M**22)

Ax = sp.diff(A72, M)   # partial wrt x
Ay = sp.diff(A72, L)   # partial wrt y

A_f  = sp.lambdify((M, L), A72, 'mpmath')
Ax_f = sp.lambdify((M, L), Ax,  'mpmath')
Ay_f = sp.lambdify((M, L), Ay,  'mpmath')

# Coefficients of A_{7_2} as a polynomial in L (highest degree first),
# each a function of M, for continuity-independent root enumeration.
polyL = sp.Poly(A72, L)
coeff_funcs = [sp.lambdify(M, c, 'mpmath') for c in polyL.all_coeffs()]

def y_roots(x):
    return mp.polyroots([f(x) for f in coeff_funcs], maxsteps=200, extraprec=200)

# --------------------------------------------------------------------------
# Endpoints.
# --------------------------------------------------------------------------
fa = lambda z: (z**12 - 3*z**11 + 4*z**10 - 5*z**9 + 6*z**8 - 7*z**7 + 7*z**6
                - 7*z**5 + 6*z**4 - 5*z**3 + 4*z**2 - 3*z + 1)
fb = lambda z: (z**16 - 7*z**15 + 22*z**14 - 48*z**13 + 87*z**12 - 133*z**11
                + 178*z**10 - 211*z**9 + 223*z**8 - 211*z**7 + 178*z**6
                - 133*z**5 + 87*z**4 - 48*z**3 + 22*z**2 - 7*z + 1)

t = mp.findroot(fa, mp.mpf('0.59098942867025644049'))
r = mp.findroot(fb, mp.mpf('0.40681308133678976238'))
alpha = t**2
beta  = r

print("=" * 74)
print("ANCHOR: direct integral of eta along the A_{7_2} real branch")
print("=" * 74)
print(f"t (root f_alpha) = {mp.nstr(t, 40)}   f_alpha(t) = {mp.nstr(fa(t),3)}")
print(f"r (root f_beta ) = {mp.nstr(r, 40)}   f_beta(r)  = {mp.nstr(fb(r),3)}")
print(f"alpha = t^2      = {mp.nstr(alpha, 40)}")
print(f"beta  = r        = {mp.nstr(beta, 40)}")
print(f"A_72(alpha, t)   = {mp.nstr(A_f(alpha, t), 3)}")
print(f"A_72(beta,  r)   = {mp.nstr(A_f(beta,  r), 3)}")

# --------------------------------------------------------------------------
# Branch selection at the endpoints: list all roots, confirm which is chosen.
# --------------------------------------------------------------------------
def real_positive_roots(x):
    out = []
    for z in y_roots(x):
        if abs(mp.im(z)) < mp.mpf(10)**(-30):
            zr = mp.re(z)
            if zr > 0:
                out.append(zr)
    return sorted(out)

print("\nBranch selection subtlety (all real positive roots of A_72(x,.)=0):")
print(f"  at x=alpha={mp.nstr(alpha,12)}: {[mp.nstr(z,14) for z in real_positive_roots(alpha)]}")
print(f"     -> challenge branch selects y = t = {mp.nstr(t,14)}")
print(f"  at x=beta ={mp.nstr(beta,12)}: {[mp.nstr(z,14) for z in real_positive_roots(beta)]}")
print(f"     -> challenge branch selects y = r = {mp.nstr(r,14)}")

# --------------------------------------------------------------------------
# Root selection: at each x, enumerate ALL 5 roots and pick the real root
# nearest to a reference value (continuity).  This is robust; single-seed
# secant iteration is not (it hops between the four nearby real roots).
# --------------------------------------------------------------------------
def nearest_real_root(x, ref):
    best = None
    bestd = None
    for z in y_roots(x):
        if abs(mp.im(z)) < mp.mpf(10)**(-20):
            zr = mp.re(z)
            d = abs(zr - ref)
            if bestd is None or d < bestd:
                bestd, best = d, zr
    return best

# Continuity-tracked table of y(x) across [alpha, beta].
N_TABLE = 800
xs = [alpha + (beta - alpha) * mp.mpf(i) / N_TABLE for i in range(N_TABLE + 1)]
ys = [t]
for i in range(1, N_TABLE + 1):
    ys.append(nearest_real_root(xs[i], ys[-1]))

# sanity: endpoint of the track lands on r
print(f"\nTracked branch endpoint y(beta) = {mp.nstr(ys[-1], 30)}")
print(f"                        expected r = {mp.nstr(r, 30)}")
print(f"   |y(beta) - r| = {mp.nstr(abs(ys[-1]-r), 3)}")

def guess_for(x):
    # nearest tabulated node (grid is uniform), then interpolate linearly
    pos = (x - alpha) / (beta - alpha) * N_TABLE
    i = int(mp.floor(pos))
    i = max(0, min(N_TABLE - 1, i))
    frac = pos - i
    return ys[i] * (1 - frac) + ys[i + 1] * frac

def y_of(x):
    # polish the interpolated guess to full precision on the tracked branch
    guess = guess_for(x)
    return nearest_real_root(x, guess)

def yprime(x, y):
    return -Ax_f(x, y) / Ay_f(x, y)

# --------------------------------------------------------------------------
# Branch character check: positive, decreasing, y' <= -2.
# --------------------------------------------------------------------------
print("\nBranch character (should be positive, decreasing, y' <= -2):")
worst_yp = mp.mpf('-inf')
mono_ok = True
for frac in [mp.mpf(k)/10 for k in range(11)]:
    x = alpha + (beta - alpha) * frac
    y = y_of(x)
    yp = yprime(x, y)
    worst_yp = max(worst_yp, yp)  # closest to 0 (least negative)
    print(f"  x={mp.nstr(x,10)}  y={mp.nstr(y,12)}  y'={mp.nstr(yp,8)}")
print(f"  max y' over samples = {mp.nstr(worst_yp,8)}  (need <= -2): "
      f"{'OK' if worst_yp <= -2 else 'VIOLATED'}")

# --------------------------------------------------------------------------
# The integrand:  eta = log x d log y - log y d log x
#   = [ log(x) * y'/y  -  log(y) * (1/x) ] dx
# --------------------------------------------------------------------------
def integrand(x):
    y = y_of(x)
    yp = yprime(x, y)
    return mp.log(x) * (yp / y) - mp.log(y) * (1 / x)

I = mp.quad(integrand, [alpha, beta])
target = 4 * mp.pi**2 / 85

print("\n" + "=" * 74)
print("RESULT")
print("=" * 74)
print(f"  I  = integral of eta      = {mp.nstr(I, 45)}")
print(f"  4*pi^2/85                 = {mp.nstr(target, 45)}")
diff = abs(I - target)
print(f"  |I - 4*pi^2/85|           = {mp.nstr(diff, 4)}")
if diff == 0:
    print("  digits of agreement       : full working precision (exact zero diff)")
else:
    print(f"  digits of agreement       : ~{int(mp.floor(-mp.log10(diff)))}")
print(f"  I / (pi^2)                = {mp.nstr(I/mp.pi**2, 30)}")
print(f"  4/85                      = {mp.nstr(mp.mpf(4)/85, 30)}")
print("=" * 74)
