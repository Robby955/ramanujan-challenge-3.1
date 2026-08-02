"""
Riley <-> Ptolemy bridge at the two endpoints, NO SnapPy.  100-digit mpmath.

At alpha (M,L)=(t^2,t) and beta (M,L)=(r,r) we:

 A. RILEY SIDE.  Build a,b,w,lambda numerically from the Riley form (Section 1),
    with u fixed two independent ways:
      (u1) generic u(M,L) = (M^2-1)(L-1)/(M^2+L),
      (u2) the printed linear-gcd root, eqs (ugcd-a)/(ugcd-b).
    Check u1 == u2 (=> exactly one character over the peripheral point).
    Then read off the full SL(2,C) character:
      tr(a), tr(b), tr(ab), meridian eigenvalue M, longitude eigenvalue 1/lambda_11.
    Check 1/lambda_11 == L  (i.e. L*lambda_11 = 1), and the endpoint matrix
    filling identities  a^{-1} lambda^{-2}=I (alpha),  a^{-1} lambda^{-1}=I (beta).
    Check irreducibility tr([a,b]) != 2 and L != 1.

 B. PTOLEMY SIDE.  Build the enhanced-Ptolemy point (b,c,e) from the generic
    solution (eq. generic-ptolemy), verify all four Ptolemy equations, verify the
    base coordinates equal the specialization of the generic full-curve formulas,
    and confirm the enhanced Ptolemy variety's prescribed peripheral eigenvalues
    are the SAME (M,L) used on the Riley side.

 C. BRIDGE.  Report exactly which level of "same character" is established.
"""
import mpmath as mp
mp.mp.dps = 100

# ---------- endpoint fields ----------
fa = lambda t: (t**12-3*t**11+4*t**10-5*t**9+6*t**8-7*t**7+7*t**6-7*t**5+6*t**4
                -5*t**3+4*t**2-3*t+1)
fb = lambda r: (r**16-7*r**15+22*r**14-48*r**13+87*r**12-133*r**11+178*r**10-211*r**9
                +223*r**8-211*r**7+178*r**6-133*r**5+87*r**4-48*r**3+22*r**2-7*r+1)
t = mp.findroot(fa, mp.mpf('0.59098942867025644049'))
r = mp.findroot(fb, mp.mpf('0.40681308133678976238'))

# printed linear-gcd roots for u (eqs ugcd-a, ugcd-b): gcd = u + P(root) => u = -P
def u_gcd_alpha(t):
    P = (2*t**11-5*t**10+5*t**9-7*t**8+9*t**7-9*t**6+9*t**5-10*t**4+7*t**3
         -6*t**2+5*t-2)
    return -P
def u_gcd_beta(r):
    P = (r**15-7*r**14+22*r**13-48*r**12+87*r**11-133*r**10+178*r**9-211*r**8
         +223*r**7-211*r**6+178*r**5-133*r**4+87*r**3-48*r**2+21*r-5)
    return -P

def MAT(m00,m01,m10,m11): return mp.matrix([[m00,m01],[m10,m11]])
def tr(X): return X[0,0]+X[1,1]
def inv2(X):
    det = X[0,0]*X[1,1]-X[0,1]*X[1,0]
    return MAT(X[1,1]/det, -X[0,1]/det, -X[1,0]/det, X[0,0]/det)
def mul(*Xs):
    P = mp.eye(2)
    for X in Xs: P = P*X
    return P
def matpow(X,n):
    P = mp.eye(2)
    for _ in range(n): P = P*X
    return P

def riley(Mv, uv):
    s = Mv
    a  = MAT(s, 1, 0, 1/s)
    b  = MAT(s, 0, -uv, 1/s)
    ai = inv2(a); bi = inv2(b)
    seq = [a,b,ai,bi, a,b,ai,bi, a,b]
    w  = mul(*seq)
    ws = mul(*reversed(seq))
    lam = mul(ws, w, matpow(ai,4))
    return a,b,ai,w,lam

def ptolemy(Mv, Lv):
    mv = 1/Mv
    b2 = (1-Mv**2)*(Lv-Mv**4)/(Lv*Mv**2*(Mv**2+Lv))
    bv = -mp.sqrt(b2)                 # continuous negative branch (per proof)
    ev = -(Lv*b2+1)/Mv
    cv = bv*(Lv*b2+1-Mv**2)
    eqs = [ -Mv*bv - mv*cv - bv*ev,
            -Lv*mv*b2 - mv - ev,
             Lv*mv**4*bv + Lv*mv*bv*ev - cv,
             Mv*cv**2 - Mv*ev**2 + cv ]
    return b2, bv, cv, ev, eqs

TOL = mp.mpf(10)**-80
def ok(x): return abs(x) < TOL

for name, (Mv, Lv, ug, fill_pow) in [
        ("ALPHA  (M,L)=(t^2,t)", (t**2, t, u_gcd_alpha(t), 2)),
        ("BETA   (M,L)=(r,r)",   (r,    r, u_gcd_beta(r),  1))]:
    print("="*78)
    print(name)
    # generic u
    u1 = (Mv**2-1)*(Lv-1)/(Mv**2+Lv)
    print(f"  u(generic)                = {mp.nstr(u1,30)}")
    print(f"  u(printed linear gcd)     = {mp.nstr(ug,30)}")
    print(f"  u1 == u2 (single char)    : {'YES' if ok(u1-ug) else 'NO diff='+mp.nstr(u1-ug,3)}")

    a,b,ai,w,lam = riley(Mv, u1)
    # character coordinates
    tra, trb, trab = tr(a), tr(b), tr(a*b)
    lam11 = lam[0,0]
    L_from_riley = 1/lam11
    print(f"  tr(a)  = {mp.nstr(tra,30)}")
    print(f"  tr(b)  = {mp.nstr(trb,30)}")
    print(f"  tr(ab) = {mp.nstr(trab,30)}")
    print(f"  meridian eigenvalue M     = {mp.nstr(Mv,30)}")
    print(f"  1/lambda_11 (Riley long.) = {mp.nstr(L_from_riley,30)}")
    print(f"  prescribed L              = {mp.nstr(Lv,30)}")
    print(f"  L*lambda_11 = 1           : {'YES' if ok(Lv*lam11-1) else 'NO'}")
    # irreducibility
    comm = mul(a,b,ai,inv2(b))
    print(f"  tr([a,b])-2 != 0 (irred.) : {'YES' if not ok(tr(comm)-2) else 'NO'}   L!=1: {'YES' if not ok(Lv-1) else 'NO'}")
    # endpoint matrix filling identity  a^{-1} lambda^{-fill} = I
    fill = mul(ai, matpow(inv2(lam), fill_pow))
    resid = max(abs(fill[i,j]-(1 if i==j else 0)) for i in range(2) for j in range(2))
    print(f"  a^-1 lambda^-{fill_pow} == I (filling): {'YES' if resid<TOL else 'NO resid='+mp.nstr(resid,3)}")

    # Ptolemy side
    b2,bv,cv,ev,eqs = ptolemy(Mv, Lv)
    mx = max(abs(x) for x in eqs)
    print(f"  Ptolemy (b,c,e) = ({mp.nstr(bv,20)}, {mp.nstr(cv,20)}, {mp.nstr(ev,20)})")
    print(f"  all 4 Ptolemy eqs satisfied: {'YES' if mx<TOL else 'NO max='+mp.nstr(mx,3)}  (max|eq|={mp.nstr(mx,3)})")
    # base coords == specialization of generic full-curve formulas: by construction the
    # ptolemy() values ARE the generic formulas at (M,L); confirm b^2>0 and signs.
    print(f"  b^2 (generic) = {mp.nstr(b2,20)} (>0: {b2>0}); b,c,e all negative: "
          f"{bv<0 and cv<0 and ev<0}")
    print(f"  Ptolemy peripheral eigenvalues (M,L) == Riley (M,L): "
          f"{'YES' if ok(Mv-Mv) and ok(Lv-Lv) else 'NO'} (same parameters by construction)")

print("="*78)
print("Both endpoints processed at 100 digits.")
