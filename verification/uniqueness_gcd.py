"""
Uniqueness of the character over each peripheral point (the logical crux of the
Riley<->Ptolemy bridge), NO SnapPy.

The bridge argument is: any irreducible rep of the two-bridge group is conjugate
to the Riley form, and over the peripheral point (M,L) the linear gcd
    gcd_u( R(M,u),  L*lambda_11(M,u)-1 )
leaves EXACTLY ONE value of u, hence exactly one character.  R = numerator of the
(1,2) relator entry (degree 5 in u); the longitude numerator has degree 9 in u.
We verify the gcd is linear (a single common root) at both endpoints, and that
this common root is the u used on both the Riley and Ptolemy sides.
"""
import sympy as sp
import mpmath as mp
mp.mp.dps = 80

s, u = sp.symbols('s u')
a  = sp.Matrix([[s,1],[0,1/s]]); ai = sp.Matrix([[1/s,-1],[0,s]])
b  = sp.Matrix([[s,0],[-u,1/s]]); bi = sp.Matrix([[1/s,0],[u,s]])
def W(L):
    P=sp.eye(2)
    for X in L: P=P*X
    return P
seq=[a,b,ai,bi,a,b,ai,bi,a,b]; w=W(seq); ws=W(list(reversed(seq)))
lam=ws*w*(ai**4)
R      = sp.expand(sp.numer(sp.cancel(sp.together((w*a)[0,1]-(b*w)[0,1]))))
lam11  = sp.cancel(sp.together(lam[0,0]))

Rfun    = sp.lambdify((s,u), R, 'mpmath')
lam11f  = sp.lambdify((s,u), lam11, 'mpmath')
Rpoly_u = sp.Poly(R, u)   # coefficients are Laurent-ish in s; evaluate numerically

fa = lambda x:(x**12-3*x**11+4*x**10-5*x**9+6*x**8-7*x**7+7*x**6-7*x**5+6*x**4-5*x**3+4*x**2-3*x+1)
fb = lambda x:(x**16-7*x**15+22*x**14-48*x**13+87*x**12-133*x**11+178*x**10-211*x**9
               +223*x**8-211*x**7+178*x**6-133*x**5+87*x**4-48*x**3+22*x**2-7*x+1)
t = mp.findroot(fa, mp.mpf('0.59098942867025644049'))
r = mp.findroot(fb, mp.mpf('0.40681308133678976238'))

_ccache = [sp.lambdify(s, c, 'mpmath') for c in Rpoly_u.all_coeffs()]
def coeffs_in_u(sval):
    # numeric coefficients of R as a univariate polynomial in u at s=sval (full precision)
    return [mp.mpc(f(sval)) for f in _ccache]

for label, sval, Lval, uref in [
        # generic u = (M^2-1)(L-1)/(M^2+L)
        ("ALPHA (M=t^2, L=t)", t**2, t, ((t**2)**2-1)*(t-1)/((t**2)**2+t)),
        ("BETA  (M=r,   L=r)", r,    r, (r**2-1)*(r-1)/(r**2+r))]:
    print("="*70)
    print(label)
    cs = coeffs_in_u(sval)
    # sanity: reference u really is a root of R and of the longitude relation
    print(f"  check R(M,uref)         = {mp.nstr(Rfun(sval,uref),3)}")
    print(f"  check L*lam11(M,uref)-1 = {mp.nstr(Lval*lam11f(sval,uref)-1,3)}")
    roots = mp.polyroots(cs, maxsteps=500, extraprec=600)
    common = []
    for ru in roots:
        longrel = Lval*lam11f(sval, ru) - 1
        if abs(longrel) < mp.mpf(10)**-30:
            common.append(ru)
    print(f"  R(M,u) has degree {len(cs)-1} in u; {len(roots)} roots")
    print(f"  common roots of R and (L*lambda_11-1): {len(common)}")
    for cr in common:
        print(f"    u* = {mp.nstr(cr,40)}")
    print(f"  reference u (generic/gcd/Riley/Ptolemy) = {mp.nstr(uref,40)}")
    if len(common)==1 and abs(common[0]-uref) < mp.mpf(10)**-30:
        print("  => gcd is LINEAR: exactly one character over the peripheral point,")
        print("     equal to the u used on both Riley and Ptolemy sides.  UNIQUE.")
    else:
        print("  => NOT a single common root -- uniqueness NOT established this way.")
print("="*70)
