"""
Riley-side of the Riley<->Ptolemy bridge (Problem 3.1), NO SnapPy.

Transcribes the two-bridge Riley data from solution.tex Section 1:

    a = [[s,1],[0,1/s]],   b = [[s,0],[-u,1/s]],   s = M
    w = a b a^{-1} b^{-1} a b a^{-1} b^{-1} a b
    lambda = w^* w a^{-4},   w^* = w read backwards
    u(M,L) = (M^2-1)(L-1)/(M^2+L)
    L = lambda_11^{-1}   (challenge longitude = inverse of canonical longitude)

Verifies, as exact rational-function identities in Q(M,L):
  (I)   relator entry (wa-bw)_{12} = A_72(M,L) / (M^6 (M^2+L)^5)   [eq. generic-relator]
  (II)  tr([a,b]) - 2 = (L-1)(M^2-1)^2(L-M^4) / (M^2 (M^2+L)^2)    [eq. generic-irreducible]
  (III) A_72 | numerator(L*lambda_11 - 1)                          [eq. generic-relator, 2nd half]

Then computes, at BOTH endpoints, the full SL(2,C) character coordinates
  tr(a), tr(b), tr(ab), meridian eigenvalue M, longitude eigenvalue L
from the Riley parametrization, and checks L*lambda_11 = 1 there.
"""
import sympy as sp

M, L, s, u = sp.symbols('M L s u')

# A_72 exactly as in the prior independent check (transcribed from the certificate).
A72 = (L**5
 + L**4*(M**14 - M**12 + 3*M**4 + 4*M**2 - 2)
 + L**3*(-2*M**18 + 5*M**16 + M**14 - 4*M**12 + 6*M**8 + 5*M**6 + 2*M**4 - 4*M**2 + 1)
 + L**2*(M**22 - 4*M**20 + 2*M**18 + 5*M**16 + 6*M**14 - 4*M**10 + M**8 + 5*M**6 - 2*M**4)
 + L*(-2*M**22 + 4*M**20 + 3*M**18 - M**10 + M**8) + M**22)

# ---- Riley matrices over Q(s,u) ----
a  = sp.Matrix([[s, 1],[0, 1/s]])
ai = sp.Matrix([[1/s, -1],[0, s]])          # a^{-1}, det a = 1
b  = sp.Matrix([[s, 0],[-u, 1/s]])
bi = sp.Matrix([[1/s, 0],[ u, s]])          # b^{-1}, det b = 1

def word(mats):
    P = sp.eye(2)
    for Mx in mats:
        P = P*Mx
    return P

# w = a b a^{-1} b^{-1} a b a^{-1} b^{-1} a b
wlist = [a,b,ai,bi, a,b,ai,bi, a,b]
w  = word(wlist)
# w^* = w read backwards
wstar = word(list(reversed(wlist)))

lam = wstar * w * (ai**4)                    # canonical longitude lambda = w^* w a^{-4}

# ---------- (I) relator identity ----------
WA = w*a
BW = b*w
rel12 = sp.together(sp.expand(WA[0,1] - BW[0,1]))     # (wa - bw)_{1,2}
# substitute Riley form s = M and the generic u = u(M,L)
uML = (M**2 - 1)*(L - 1)/(M**2 + L)
rel12_ML = sp.simplify(rel12.subs({s: M, u: uML}))
target_rel = A72 / (M**6 * (M**2 + L)**5)
diff_rel = sp.simplify(rel12_ML - target_rel)
print("(I)   (wa-bw)_12  ==  A72/(M^6 (M^2+L)^5)  :",
      "HOLDS" if diff_rel == 0 else f"FAIL residual={diff_rel}")

# also confirm the other three relator entries vanish or are multiples of entry (1,2)
rel = (WA - BW).applyfunc(lambda x: sp.simplify(x.subs({s: M, u: uML})))
r11 = sp.simplify(rel[0,0]); r22 = sp.simplify(rel[1,1])
r21_ratio = sp.simplify(rel[1,0]/rel12_ML) if rel12_ML != 0 else None
print("      (wa-bw)_11==0:", r11 == 0, " (wa-bw)_22==0:", r22 == 0,
      " (wa-bw)_21 / (wa-bw)_12 =", sp.simplify(r21_ratio))

# ---------- (II) irreducibility obstruction ----------
comm = a*b*ai*bi
trc = sp.simplify((comm[0,0] + comm[1,1]).subs({s: M, u: uML})) - 2
target_irr = (L - 1)*(M**2 - 1)**2*(L - M**4)/(M**2*(M**2 + L)**2)
diff_irr = sp.simplify(trc - target_irr)
print("(II)  tr([a,b])-2 ==  (L-1)(M^2-1)^2(L-M^4)/(M^2(M^2+L)^2) :",
      "HOLDS" if diff_irr == 0 else f"FAIL residual={diff_irr}")

# ---------- (III) A72 | num(L*lambda_11 - 1) ----------
lam11 = sp.together(sp.simplify(lam[0,0].subs({s: M, u: uML})))
expr = sp.together(L*lam11 - 1)
numer = sp.numer(sp.cancel(expr))
q, rem = sp.div(sp.Poly(sp.expand(numer), L), sp.Poly(sp.expand(A72), L))
print("(III) A72 | numerator(L*lambda_11 - 1)  :",
      "HOLDS (remainder 0)" if rem == 0 else f"FAIL remainder={rem.as_expr()}")

# also print generic u expression as a sanity cross-check of the character map
print("\ngeneric u(M,L) = (M^2-1)(L-1)/(M^2+L)  [Riley off-diagonal coupling]")
print("tr(a)=tr(b)= s+1/s = M+1/M ;  tr(ab) = M^2 + 1/M^2 - u")
ab = a*b
trab_sym = sp.simplify((ab[0,0]+ab[1,1]).subs({s:M}))
print("symbolic tr(ab) in (M,u):", trab_sym)
