import numpy as np
Rg = 8.314; Tref = 298.15
Ecal, EH, EL, betaL, cref = 58e3, 33e3, 55e3, 2.64, 0.5
kcal = 0.0275/np.sqrt(365)      # per sqrt(day), q_cal(1y,25C,50%SOC)=2.75%
def A(T, E): return np.exp(-E/Rg*(1/T - 1/Tref))
def S(x): return 1/(1+np.exp(-x))
def gsoc(z): return np.interp(z, [0, .5, 1.0], [0.67, 1.0, 1.83])
def B(T, c): return 0.5*(A(T, EH) + np.exp(betaL*(c - cref))/A(T, EL))
print("g(0),g(.5),g(1)=", gsoc(0), gsoc(.5), gsoc(1.0))
print("B(25C,0.5C)=%.3f B(10C,1C)=%.3f B(0C,1C)=%.3f B(35C,0.5C)=%.3f" % (B(298.15, .5), B(283.15, 1), B(273.15, 1), B(308.15, .5)))
kcyc = (0.20 - kcal*np.sqrt(1000))/np.sqrt(1000)
print("kcal=%.3e kcyc=%.3e" % (kcal, kcyc))
def day_profile(strategy, Tamb, Eday_frac, chg_c=0.7):
    """returns list of (dt_h, z, T_K, dFCE/dt (1/h), c_ch) segments for one day"""
    seg = []; Tk = Tamb + 273.15
    zmax = 0.8 if strategy == 'limit80' else 1.0
    zlow = zmax - Eday_frac
    tch = (zmax - zlow)/chg_c   # h of charging
    # charging segment at night (T + 8 K)
    n = 10
    for i in range(n):
        seg.append((tch/n, zlow + (zmax - zlow)*(i+.5)/n, Tk + 8, 0.0, chg_c))
    hold = 8.0 - tch
    if strategy == 'obc':
        seg.append((hold-1, 0.8, Tk, 0, 0)); seg.append((1, 1.0, Tk, 0, 0))
    else:
        seg.append((hold, zmax, Tk, 0, 0))
    for i in range(16):
        seg.append((1.0, zmax - Eday_frac*(i+.5)/16, Tk + 3, Eday_frac/16, 0))
    return seg
def run(strategy, Tamb, years=3, Pbar=0.45, En=18.19):
    qc, qy = 1e-6, 1e-6; out = []
    Eday = Pbar*24/En
    if strategy == 'obc': pass
    for d in range(int(365*years)):
        for dt, z, T, fce_rate, c in day_profile(strategy, Tamb, Eday*(1/(1 - qc - qy))):
            # calendar: dq/dt = (kcal A g)^2/(2q), t in days
            qc += (kcal*A(T, Ecal)*gsoc(z))**2/(2*qc)*dt/24
            # cycle (discharge FCE); charge plating multiplier via B(T, c) on charge segments too
            if fce_rate > 0:
                qy += (kcyc*B(T, 0.5))**2*fce_rate*dt/(2*qy)
            if c > 0:
                # plating excess only: (B(T,c) - B(T,cref)) contributes during charge, applied to charge throughput
                ex = max(B(T, c) - B(T, cref), 0.0)
                qy += (kcyc*ex)**2*c*dt/(2*qy)
        if d % 182 == 181: out.append((d+1, 1 - qc - qy, qc, qy))
    return out
for strat in ['full100', 'obc', 'limit80']:
    for Ta in [25, 32]:
        r = run(strat, Ta)
        print("%-8s Ta=%d " % (strat, Ta) + " ".join("d%d:SOH=%.3f(cal %.3f cyc %.3f)" % x for x in r[1::2]))
