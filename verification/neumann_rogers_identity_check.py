"""
Independent check of the Neumann extended-Rogers identity used in solution.tex eq (dS):

    2 dR(z) = log z d log(1-z) - log(1-z) d log z.

This is a pure one-form identity in the single variable z and does NOT depend on the
7_2 triangulation.  We verify it symbolically for the Rogers dilogarithm
    R(z) = Li_2(z) + (1/2) log z log(1-z),
i.e. that  2 R'(z) = -log z/(1-z) - log(1-z)/z  as rational-log functions,
and confirm numerically at several points.
"""
import sympy as sp

z = sp.symbols('z')
R = sp.polylog(2, z) + sp.Rational(1, 2)*sp.log(z)*sp.log(1 - z)

# LHS as a function multiplying dz:  2 R'(z)
lhs = sp.simplify(2*sp.diff(R, z))

# RHS as a function multiplying dz:
#   log z d log(1-z) - log(1-z) d log z
#   = log z * (-1/(1-z)) dz  -  log(1-z) * (1/z) dz
rhs = sp.simplify(sp.log(z)*(-1/(1 - z)) - sp.log(1 - z)*(1/z))

# sympy leaves Li_1(z) unexpanded; Li_1(z) = -log(1-z), so expand it before comparing
diff = sp.simplify(sp.expand_func(lhs - rhs))
print("Neumann/Rogers identity  2 dR = log z dlog(1-z) - log(1-z) dlog z")
print("  2 R'(z)          =", lhs)
print("  RHS coefficient  =", rhs)
print("  symbolic (LHS-RHS) simplifies to:", diff)
print("  SYMBOLIC:", "PASS" if diff == 0 else "FAIL")

# numeric cross-check with complex z off the branch cuts
import mpmath as mp
mp.mp.dps = 40
def Rn(zz):
    return mp.polylog(2, zz) + mp.mpf('0.5')*mp.log(zz)*mp.log(1 - zz)
maxerr = mp.mpf(0)
for zz in [mp.mpf('0.3'), mp.mpf('0.7'), mp.mpc('0.2','0.5'), mp.mpc('0.6','-0.3')]:
    h = mp.mpf(10)**-20
    dR = (Rn(zz + h) - Rn(zz - h))/(2*h)           # 2R' numerically -> compare 2*dR
    lhsn = 2*dR
    rhsn = -mp.log(zz)/(1 - zz) - mp.log(1 - zz)/zz
    maxerr = max(maxerr, abs(lhsn - rhsn))
print("  NUMERIC max|2R'(z) - RHS| over 4 points =", mp.nstr(maxerr, 3),
      "->", "PASS" if maxerr < mp.mpf(10)**-15 else "FAIL")
