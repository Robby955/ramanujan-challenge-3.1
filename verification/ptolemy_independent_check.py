import sympy as sp, mpmath as mp
mp.mp.dps = 60

# ---------- symbolic: correct eq4 condition (even^2 - b2*odd^2 == 0 mod A) ----------
M, L, b = sp.symbols('M L b')
m = 1/M
A72 = (L**5
 + L**4*(M**14 - M**12 + 3*M**4 + 4*M**2 - 2)
 + L**3*(-2*M**18 + 5*M**16 + M**14 - 4*M**12 + 6*M**8 + 5*M**6 + 2*M**4 - 4*M**2 + 1)
 + L**2*(M**22 - 4*M**20 + 2*M**18 + 5*M**16 + 6*M**14 - 4*M**10 + M**8 + 5*M**6 - 2*M**4)
 + L*(-2*M**22 + 4*M**20 + 3*M**18 - M**10 + M**8) + M**22)
b2 = (1 - M**2)*(L - M**4)/(L*M**2*(M**2 + L))
e  = -(L*b**2 + 1)/M
c  = b*(L*b**2 + 1 - M**2)
eq4 = sp.expand(M*c**2 - M*e**2 + c)
pr = sp.Poly(eq4, b)
even = sum(co*b2**(k//2) for (k,),co in pr.terms() if k%2==0)
odd  = sum(co*b2**((k-1)//2) for (k,),co in pr.terms() if k%2==1)
cond = sp.cancel(even**2 - b2*odd**2)               # = 0  iff  even + b*odd = 0 for b^2=b2
num  = sp.expand(sp.numer(cond))
rem  = sp.rem(sp.Poly(num, L), sp.Poly(A72, L))
print("eq4 correct symbolic condition (even^2 - b2*odd^2) reduced mod A_7_2:",
      "== 0" if sp.simplify(rem.as_expr())==0 else "NONZERO")

# ---------- numerical: all four Ptolemy eqs at BOTH endpoints, 60 digits ----------
fa = lambda t: t**12-3*t**11+4*t**10-5*t**9+6*t**8-7*t**7+7*t**6-7*t**5+6*t**4-5*t**3+4*t**2-3*t+1
fb = lambda r: (r**16-7*r**15+22*r**14-48*r**13+87*r**12-133*r**11+178*r**10-211*r**9
                +223*r**8-211*r**7+178*r**6-133*r**5+87*r**4-48*r**3+22*r**2-7*r+1)
t = mp.findroot(fa, mp.mpf('0.59098942867'))
r = mp.findroot(fb, mp.mpf('0.40681308133'))
print(f"\nendpoint roots: t={mp.nstr(t,20)} (f_a={mp.nstr(fa(t),3)}), r={mp.nstr(r,20)} (f_b={mp.nstr(fb(r),3)})")

def ptolemy_eqs(Mv, Lv):
    mv = 1/Mv
    b2v = (1-Mv**2)*(Lv-Mv**4)/(Lv*Mv**2*(Mv**2+Lv))
    bv  = -mp.sqrt(b2v)                    # negative continuous branch (per the proof)
    ev  = -(Lv*b2v+1)/Mv
    cv  = bv*(Lv*b2v+1-Mv**2)
    return [ -Mv*bv - mv*cv - bv*ev,
             -Lv*mv*b2v - mv - ev,
              Lv*mv**4*bv + Lv*mv*bv*ev - cv,
              Mv*cv**2 - Mv*ev**2 + cv ]

for name,(Mv,Lv) in [("alpha  (M,L)=(t^2,t)",(t**2,t)), ("beta   (M,L)=(r,r)",(r,r))]:
    res = ptolemy_eqs(Mv,Lv)
    mx = max(abs(x) for x in res)
    print(f"  {name}: max|eq_i| = {mp.nstr(mx,3)}   ->  {'ALL FOUR SATISFIED' if mx < mp.mpf(10)**-45 else 'FAIL'}")
