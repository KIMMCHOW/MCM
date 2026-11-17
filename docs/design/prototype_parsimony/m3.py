import numpy as np, sys
sys.path.insert(0, '.')
from proto import *  # noqa
zt = np.array([0, .05, .1, .15, .2, .25, .3, .35, .4, .45, .5, .55, .6, .65, .7, .75, .8, .85, .9, .95, 1.0])
dudt = 1e-3*np.array([-.4485, -.7057, -.5343, -.4982, -.4787, -.4213, -.3722, -.3439, -.3068, -.2655, -.2482,
                      -.2412, -.2482, -.2281, -.2047, -.1798, -.1528, -.1246, -.0966, -.0705, -.0459])
def dUdT(z): return np.interp(z, zt, dudt)
Cth, Rth, Rsk = 150.0, 3.0, 5.5
Tthr, wthr, Psoc_floor = 43.0, 1.0, 0.3

def simulate(schedule, Tamb_C, z0=1.0, qloss=0.0, rg=0.0, h=2.0, tmax=400*3600, T0_C=None,
             rc=True, thermal=True, throttle=True, entropic=True, record=False):
    """schedule(t)-> (P_nonsoc, P_soc_dem). returns TTE_h, z_end, info"""
    Ta = Tamb_C + 273.15
    T = Ta if T0_C is None else T0_C + 273.15
    z, v1 = z0, 0.0
    Qe = Q*(1-qloss)
    t = 0.0
    hist = []
    def rhs(t, z, v1, T):
        Pns, Pso = schedule(t)
        Tsk = Tamb_C + ((T-273.15) - Tamb_C)*Rsk/Rth
        if throttle and Pso > Psoc_floor:
            phi = 1/(1+np.exp(-(Tthr - Tsk)/wthr))
            Pso = Psoc_floor + (Pso - Psoc_floor)*phi
        P = Pns + Pso
        Tk = T if thermal else Ta
        r0 = R0(z, Tk, rg); r1 = R1(Tk, rg); tau = tau1ref*arr(Eatau, Tk)
        Ue = U(z) + ((Tk - Tref)*dUdT(z) if entropic else 0.0) - (v1 if rc else 0.0)
        Rs = r0 if rc else r0 + r1
        disc = Ue*Ue - 4*Rs*P
        if disc <= 0: return None, P, Ue, np.nan
        I = (Ue - np.sqrt(disc))/(2*Rs)
        V = Ue - I*Rs
        dz = -I/(3600*Qe)
        dv = (-v1/tau + I*r1/tau) if rc else 0.0
        qheat = P + I*I*Rs + ((v1*v1/r1) if rc else 0.0) - (I*Tk*dUdT(z) if entropic else 0.0)
        dT = (qheat - (T - Ta)/Rth)/Cth if thermal else 0.0
        return (dz, dv, dT), P, V, I
    while t < tmax:
        k1, P, V, I = rhs(t, z, v1, T)
        if k1 is None or V <= Vcut:
            return t/3600, z, dict(T=T-273.15, P=P, V=V, hist=hist)
        if record and int(t) % 60 == 0: hist.append((t, z, v1, T-273.15, P, V))
        k2, *_ = rhs(t+h/2, z+h/2*k1[0], v1+h/2*k1[1], T+h/2*k1[2])
        if k2 is None: return t/3600, z, dict(T=T-273.15, P=P, V=V, hist=hist)
        k3, *_ = rhs(t+h/2, z+h/2*k2[0], v1+h/2*k2[1], T+h/2*k2[2])
        if k3 is None: return t/3600, z, dict(T=T-273.15, P=P, V=V, hist=hist)
        k4, *_ = rhs(t+h, z+h*k3[0], v1+h*k3[1], T+h*k3[2])
        if k4 is None: return t/3600, z, dict(T=T-273.15, P=P, V=V, hist=hist)
        z += h/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]); v1 += h/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1]); T += h/6*(k1[2]+2*k2[2]+2*k3[2]+k4[2])
        t += h
    return np.inf, z, {}

def const(act, rat='4g', **kw):
    P, comp = state_power(act, rat, **kw)
    ps = comp.get('soc', 0.0) if act == 'game' else 0.0
    return lambda t: (P - ps, ps)

if __name__ == '__main__':
    res = {}
    for name, sch, Ta, extra in [
        ('social_wifi', const('social', 'wifi'), 25, {}),
        ('video_wifi', const('video', 'wifi'), 25, {}),
        ('video_5g', const('video', '5g'), 25, {}),
        ('nav_4g', const('nav', '4g'), 25, {}),
        ('vcall_4g', const('vcall', '4g'), 25, {}),
        ('game_4g_25C', const('game', '4g'), 25, {}),
        ('game_4g_35C', const('game', '4g'), 35, {}),
        ('game_4g_15C', const('game', '4g'), 15, {}),
        ('game_nothrottle', const('game', '4g'), 25, dict(throttle=False)),
        ('nav_4g_b1_-10C', const('nav', '4g', b=1.0), -10, {}),
        ('nav_4g_b1_25C', const('nav', '4g', b=1.0), 25, {}),
        ('social_4g_weak', const('social', '4g', rsrp=-110), 25, {}),
        ('social_4g', const('social', '4g'), 25, {}),
    ]:
        tte, z, info = simulate(sch, Ta, h=2.0 if 'game' in name or 'vcall' in name else 5.0, **extra)
        P0 = sum(sch(0))
        print("%-18s Ta=%3d  P_dem=%.3f  TTE=%.2f h  z_end=%.4f  T_end=%.1f C  P_end=%.3f" % (name, Ta, P0, tte, z, info.get('T', np.nan), info.get('P', np.nan)))
    # bursty cold: 9 min nav (b=1) + 1 min camera (b=1)
    Pn, _ = state_power('nav', '4g', b=1.0); Pc, _ = state_power('camera', '4g', b=1.0)
    print("nav b1 P=%.3f camera b1 P=%.3f" % (Pn, Pc))
    burst = lambda t: ((Pc if (t % 600) >= 540 else Pn), 0.0)
    for Ta, qloss, rg in [(25, 0, 0), (0, 0, 0), (-10, 0, 0), (-20, 0, 0), (25, .2, .5), (0, .2, .5), (-10, .2, .5)]:
        tte, z, info = simulate(burst, Ta, qloss=qloss, rg=rg, h=2.0)
        zc_ref = 0.0263
        print("burst Ta=%4d SOH=%.1f TTE=%.2f h z_off=%.3f disp_SOC_off=%.1f%% T=%.1f" % (Ta, 1-qloss, tte, z, 100*(z-zc_ref)/(1-zc_ref), info['T']))
    # bigger burst: 2.5 W camera+vcall style 2 min every 10 min
    for Ta, qloss, rg in [(-10, 0, 0), (-10, .2, .5), (-20, 0, 0), (-20, .2, .5)]:
        burst2 = lambda t: ((4.0 if (t % 600) >= 570 else 0.6), 0.0)
        tte, z, info = simulate(burst2, Ta, qloss=qloss, rg=rg, h=1.0)
        print("burst4W Ta=%4d SOH=%.1f TTE=%.2f h z_off=%.3f disp=%.1f%% T=%.1f" % (Ta, 1-qloss, tte, z, 100*(z-0.0263)/(1-0.0263), info['T']))
