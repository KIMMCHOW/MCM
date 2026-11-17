"""Prototype for the parsimony/stochastic design spec. Numbers here feed the spec's 'expected values'."""
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

Rg = 8.314; Tref = 298.15
# ---------- cell ----------
Q = 4.70  # Ah
b0, b1, a0, a1, a2, a3 = 0.60117, 36.426, 3.6193, 0.44311, -0.48605, 0.63127
dV, pst = 0.242, 3.86
def U(z):
    z = np.clip(z, 0, 1)
    return -b0*np.exp(-b1*z) + a0 + a1*z + a2*z**2 + a3*z**3 + dV*z**pst
def dUdz(z):
    return b0*b1*np.exp(-b1*z) + a1 + 2*a2*z + 3*a3*z**2 + dV*pst*z**(pst-1)
R0ref, A_R, z0R = 0.060, 0.80, 0.139
R1ref, tau1ref = 0.010, 30.0
Ea0, Ea1, Eatau = 18e3, 30e3, 21e3
Vcut = 3.40
norm = 1 + A_R*np.exp(-0.5/z0R)
def arr(Ea, T): return np.exp(Ea/Rg*(1/T - 1/Tref))
def R0(z, T=Tref, rg=0.0): return R0ref*(1 + A_R*np.exp(-z/z0R))/norm*arr(Ea0, T)*(1+rg)
def R1(T=Tref, rg=0.0): return R1ref*arr(Ea1, T)*(1+rg)

print("U(0),U(0.5),U(1):", U(0), U(0.5), U(1), " mean U:", quad(U, 0, 1)[0])
En = Q*quad(U, 0, 1)[0]
print("E_n (0->1) Wh = %.3f" % En)
zc0 = brentq(lambda z: U(z) - Vcut, 1e-4, 0.2)
print("z_c(P->0) = %.4f" % zc0)
Eu0 = Q*quad(U, zc0, 1)[0]
print("E_u(P->0) Wh = %.3f  (fraction %.4f)" % (Eu0, Eu0/En))
for z in [0.03, 0.1, 0.2, 0.5, 0.8, 0.95]:
    print(" z=%.2f U=%.4f dU/dz=%.3f R0=%.4f" % (z, U(z), dUdz(z), R0(z)))

def M1_TTE(P, T=Tref, rg=0.0, qloss=0.0, z_start=1.0, Rextra=True):
    """exact quadrature TTE (h) for constant P with series DC resistance R0+R1 (quasi-steady)"""
    Qe = Q*(1-qloss)
    Rs = lambda z: R0(z, T, rg) + (R1(T, rg) if Rextra else 0)
    f = lambda z: U(z) - Vcut - P*Rs(z)/Vcut
    if f(z_start) <= 0: return 0.0, z_start
    zc = brentq(f, 1e-5, z_start) if f(1e-5) < 0 else 1e-5
    Vt = lambda z: 0.5*(U(z) + np.sqrt(max(U(z)**2 - 4*Rs(z)*P, 0)))
    Eu = Qe*quad(Vt, zc, z_start, limit=200)[0]
    return Eu/P, zc

# ---------- components ----------
A_disp = (6.3*2.54)**2 * (19.5*9)/(19.5**2+9**2)  # cm^2
pA = 4.8e-3
print("display area cm2 = %.1f, pA*A = %.4f W" % (A_disp, pA*A_disp))
def g(b): return np.where(b <= 0.5, 0.16 + 1.68*b, 1.0 + 2.0*(b-0.5))
def Pdisp(b, apl, fr=60, on=1):
    return on*(0.12 + pA*A_disp*g(b)*apl + 0.0028*max(fr-60, 0))
P_sleep = 0.015; dP_awake = 0.080; floor_on = P_sleep + dP_awake
P_page = {'air': 0.0, 'wifi': 0.005, '4g': 0.020, '5g': 0.030}
tail = {'4g': (0.120, 10.2), '5g': (0.255, 10.3), 'wifi': (0.12, 0.24)}
def P_bg(rat, lam=1/60, Ewake_cpu=0.1):
    Pt, Tt = tail[rat]
    pid = P_page[rat]
    return pid + (Pt - pid)*(1-np.exp(-lam*Tt)) + lam*Ewake_cpu
def msig(rsrp): return 1 + 0.26*max(0, (-85 - rsrp)/10)
net = {  # activity radio power at good signal (W)
 'social': {'wifi': 0.10, '4g': 0.35, '5g': 0.50},
 'video':  {'wifi': 0.10, '4g': 0.28, '5g': 0.39},
 'game':   {'wifi': 0.03, '4g': 0.12, '5g': 0.25},
 'nav':    {'wifi': 0.05, '4g': 0.20, '5g': 0.30},
 'vcall':  {'wifi': 0.35, '4g': 1.50, '5g': 1.80},
 'camera': {'wifi': 0.02, '4g': 0.05, '5g': 0.08},
 'audio':  {'wifi': 0.02, '4g': 0.06, '5g': 0.09},
}
def state_power(act, rat='4g', b=0.4, rsrp=-85, fr_force=None, dark=False):
    m = msig(rsrp) if rat in ('4g', '5g') else 1.0
    if act == 'idle':
        return P_sleep + P_bg(rat)*(m if rat != 'wifi' else 1.0), {}
    comp = {}
    comp['floor'] = floor_on
    if act == 'social':
        apl = 0.15 if dark else 0.5; fr = 120
        comp['disp'] = Pdisp(b, apl, fr if fr_force is None else fr_force); comp['soc'] = 0.45
    elif act == 'video':
        comp['disp'] = Pdisp(b, 0.25, 60); comp['soc'] = 0.20
    elif act == 'game':
        comp['disp'] = Pdisp(max(b, 0.5), 0.35, 60 if fr_force is None else fr_force); comp['soc'] = 3.0
    elif act == 'nav':
        comp['disp'] = Pdisp(max(b, 0.6), 0.5, 60); comp['soc'] = 0.6; comp['gnss'] = 0.27
    elif act == 'vcall':
        comp['disp'] = Pdisp(b, 0.4, 60); comp['soc'] = 0.6; comp['cam'] = 0.30
    elif act == 'camera':
        comp['disp'] = Pdisp(max(b, 0.6), 0.4, 60); comp['cam'] = 1.15
    elif act == 'audio':
        comp['disp'] = 0.0; comp['soc'] = 0.03; comp['bt'] = 0.04
    comp['net'] = net[act][rat]*m
    return sum(comp.values()), comp

acts = ['idle', 'social', 'video', 'game', 'nav', 'vcall', 'camera', 'audio']
print("\nState powers (W) at b=0.4, good signal:")
for a in acts:
    row = [state_power(a, r)[0] for r in ['wifi', '4g', '5g']]
    print(" %-7s wifi %.3f  4g %.3f  5g %.3f" % (a, *row), state_power(a, '4g')[1])

# ---------- scenario constant-power TTEs (M0, M1) ----------
print("\nM0 vs M1 TTE from z=1 (25C, new):")
for name, P in [('standby4g', state_power('idle', '4g')[0]), ('social_wifi', state_power('social', 'wifi')[0]),
                ('video_wifi', state_power('video', 'wifi')[0]), ('video_5g', state_power('video', '5g')[0]),
                ('nav_4g', state_power('nav', '4g')[0]), ('vcall_4g', state_power('vcall', '4g')[0]),
                ('game_cap', 3.4), ('social_4g_weak', state_power('social', '4g', rsrp=-110)[0]),
                ('social_4g', state_power('social', '4g')[0])]:
    t0 = Eu0/P
    t1, zc = M1_TTE(P)
    print(" %-15s P=%.3f W  M0 %.2f h  M1 %.2f h  (%.2f%%) zc=%.4f" % (name, P, t0, t1, 100*(t1/t0-1), zc))
