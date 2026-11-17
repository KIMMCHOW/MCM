import numpy as np, sys, io, contextlib
sys.path.insert(0, '.')
with contextlib.redirect_stdout(io.StringIO()):
    import proto as pr
    from ctmc import acts, gen_activity, asym, product
from ctmc2 import prof
rats = ['wifi', '4g', '5g']
def chain(kind, cell='4g', wifi_mean=2.5, cell_mean=2.5):
    lam, p, d = prof(kind)
    Lc = np.array([[-1/wifi_mean, 1/wifi_mean], [1/cell_mean, -1/cell_mean]])
    return product(gen_activity(lam, p, d), Lc), [('wifi'), cell]
def Pvec(rr, cap=True, **kw):
    out = []; comps = []
    for a in acts:
        for r in rr:
            P, c = pr.state_power(a, r, **kw)
            if a == 'game' and cap:   # quasi-static thermal cap at 25C: scale soc to reach ~2.9 W
                ex = P - 2.9
                if ex > 0: c = dict(c); c['soc'] = c['soc'] - ex; P = 2.9
            if a == 'idle': c = {'sleep+bg': P}
            out.append(P); comps.append(c)
    return np.array(out), comps
def summary(kind, label, cell='4g', **kw):
    L, rr = chain(kind, cell)
    Pv, comps = Pvec(rr, **kw)
    pi, Pbar, s2, h = asym(L, Pv)
    sh = {}
    for w, c in zip(pi, comps):
        for k, v in c.items(): sh[k] = sh.get(k, 0) + w*v
    return Pbar, pr.Eu0/Pbar, {k: v/Pbar for k, v in sh.items()}, s2
base = summary('typical', 'base')
print("typical base: Pbar=%.4f TTE=%.2f h" % base[:2]); print(" component shares (=-elasticities):", {k: round(v, 3) for k, v in sorted(base[2].items(), key=lambda x: -x[1])})
hv = summary('heavy', 'base'); print("heavy base: Pbar=%.4f TTE=%.2f h" % hv[:2]); print(" shares:", {k: round(v, 3) for k, v in sorted(hv[2].items(), key=lambda x: -x[1])})
# behaviour changes (M0 closed form)
import copy
def with_patch(fn):
    saved = {k: getattr(pr, k) for k in ['P_bg', 'net', 'P_page', 'tail', 'pA']}
    try:
        fn(); return summary('typical', 'x')
    finally:
        for k, v in saved.items(): setattr(pr, k, v)
cases = {}
cases['brightness 0.4->0.3'] = summary('typical', '', b=0.3)
cases['brightness 0.4->0.6 (penalty)'] = summary('typical', '', b=0.6)
cases['dark mode (social UI)'] = summary('typical', '', dark=True)
cases['cap refresh 60 Hz'] = summary('typical', '', fr_force=60)
cases['5G instead of 4G (penalty)'] = summary('typical', '', cell='5g')
cases['weak signal -110 dBm (penalty)'] = summary('typical', '', rsrp=-110)
def halve_bg():
    f = pr.P_bg; pr.P_bg = lambda rat, lam=1/120, Ewake_cpu=0.1: f(rat, lam, Ewake_cpu)
cases['background sync 1/60s -> 1/120s'] = with_patch(halve_bg)
for k, v in cases.items():
    print(" %-34s Pbar=%.4f TTE=%.2f h  dTTE=%+.2f h (%+.1f%%)" % (k, v[0], v[1], v[1]-base[1], 100*(v[1]/base[1]-1)))
# Wi-Fi share: wifi mean 4h vs cell 1h
L, rr = chain('typical', '4g', 4.0, 1.0); Pv, _ = Pvec(rr); pi, Pbar, *_ = asym(L, Pv)
print(" %-34s Pbar=%.4f TTE=%.2f h dTTE=%+.2f h" % ('Wi-Fi share 50%->80%', Pbar, pr.Eu0/Pbar, pr.Eu0/Pbar - base[1]))
# no gaming: remove game sessions (renormalise p)
from ctmc import gen_activity as ga
lam, p, d = prof('typical'); p2 = dict(p); p2['game'] = 0; s = sum(p2.values()); p2 = {k: v/s for k, v in p2.items()}
Lc = np.array([[-1/2.5, 1/2.5], [1/2.5, -1/2.5]]); L = product(ga(lam, p2, d), Lc); Pv, _ = Pvec(['wifi', '4g'])
pi, Pbar, *_ = asym(L, Pv); print(" %-34s Pbar=%.4f TTE=%.2f h dTTE=%+.2f h" % ('no gaming sessions', Pbar, pr.Eu0/Pbar, pr.Eu0/Pbar-base[1]))
