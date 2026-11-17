import numpy as np, sys, io, contextlib
sys.path.insert(0, '.')
with contextlib.redirect_stdout(io.StringIO()):
    from proto import state_power, Eu0, U, R0, R1, arr, Eatau, tau1ref, Q, Vcut, Tref
    from m3 import dUdT, Cth, Rth, Rsk, Tthr, wthr, Psoc_floor
from ctmc import acts, gen_activity, stationary, asym, product
from scipy.stats import invgauss
rats = ['wifi', '4g']
Lc = np.array([[-1/2.5, 1/2.5], [1/2.5, -1/2.5]])
def prof(kind):
    if kind == 'typical':
        lam = 4.0; p = dict(social=.55, video=.12, game=.05, nav=.04, vcall=.04, camera=.08, audio=.12)
        d = dict(social=5, video=15, game=20, nav=20, vcall=10, camera=2, audio=30)
    else:
        lam = 4.0; p = dict(social=.45, video=.20, game=.15, nav=.05, vcall=.05, camera=.05, audio=.05)
        d = dict(social=8, video=25, game=30, nav=25, vcall=15, camera=3, audio=30)
    return lam, p, {k: v/60 for k, v in d.items()}
Pfull = {(a, r): state_power(a, r) for a in acts for r in rats}
Pv_dem = np.array([Pfull[(a, r)][0] for a in acts for r in rats])
Psoc_dem = np.array([(Pfull[(a, r)][1].get('soc', 0.0) if a == 'game' else 0.0) for a in acts for r in rats])
Pv_cap = Pv_dem.copy(); gi = [acts.index('game')*2, acts.index('game')*2+1]; Pv_cap[gi] = [2.85, 2.95]
if __name__ == "__main__":
  for kind in ["typical", "heavy"]:
      lam, p, d = prof(kind)
      L = product(gen_activity(lam, p, d), Lc)
      pi, Pbar, s2, h = asym(L, Pv_cap)
      Eu = Eu0
      hP = np.sum(pi*Pv_cap*h)/Pbar
      mu0 = Eu/Pbar; mu_b = (Eu - h[0] + hP)/Pbar
      sd = np.sqrt(s2*Eu/Pbar**3)
      print("\n== %s == Pbar=%.4f  s2=%.4f  mu0=%.2f  mu(Dynkin)=%.2f  sd=%.2f  CV=%.3f  screen-on frac=%.3f" % (
          kind, Pbar, s2, mu0, mu_b, sd, sd/mu0, 1 - pi.reshape(-1, 2)[[0, 7]].sum()))
      lamIG = Eu**2/s2; ig = invgauss(mu=mu_b/lamIG, scale=lamIG)
      print(" IG q5/50/95:", np.round(ig.ppf([.05, .5, .95]), 2))
      # ---- path sampler on grid
      rng = np.random.default_rng(7)
      N = 2000; h_s = 5.0; nmax = int(120*3600/h_s)
      rates = -np.diag(L); J = (L - np.diag(np.diag(L)))/rates[:, None]; cum = np.cumsum(J, 1)
      S = np.zeros((nmax, N), dtype=np.int8)
      for k in range(N):
          s = 0; t = 0.0; i0 = 0
          while i0 < nmax:
              tau = rng.exponential(1/rates[s])*3600; t += tau
              i1 = min(int(t/h_s), nmax); S[i0:i1, k] = s; i0 = i1
              s = np.searchsorted(cum[s], rng.random()*cum[s, -1], side='right')
      # M0 on grid
      E = np.cumsum(Pv_cap[S]*h_s/3600, 0)
      tte0 = (np.argmax(E >= Eu, 0))*h_s/3600
      # M3 vectorized RK4
      Ta = 25.0; TaK = Ta + 273.15
      z = np.ones(N); v1 = np.zeros(N); T = np.full(N, TaK); alive = np.ones(N, bool); tte3 = np.full(N, np.nan); zoff = np.full(N, np.nan)
      def rhs(z, v1, T, s):
          Pns = Pv_dem[s] - Psoc_dem[s]; Pso = Psoc_dem[s]
          Tsk = Ta + (T - 273.15 - Ta)*Rsk/Rth
          phi = 1/(1+np.exp(-(Tthr - Tsk)/wthr))
          Pso = np.where(Pso > Psoc_floor, Psoc_floor + (Pso - Psoc_floor)*phi, Pso)
          P = Pns + Pso
          r0 = R0(z, T); r1 = R1(T); tau = tau1ref*arr(Eatau, T)
          Ue = U(z) + (T - Tref)*dUdT(z) - v1
          disc = Ue*Ue - 4*r0*P
          I = (Ue - np.sqrt(np.maximum(disc, 0)))/(2*r0)
          V = np.where(disc > 0, Ue - I*r0, 0.0)
          q = P + I*I*r0 + v1*v1/r1 - I*T*dUdT(z)
          return -I/(3600*Q), -v1/tau + I*r1/tau, (q - (T - TaK)/Rth)/Cth, V
      for n in range(nmax):
          s = S[n]
          a1, b1_, c1, V = rhs(z, v1, T, s)
          dead = alive & (V <= Vcut)
          tte3[dead] = n*h_s/3600; zoff[dead] = z[dead]; alive &= ~dead
          if not alive.any(): break
          a2, b2, c2, _ = rhs(z+h_s/2*a1, v1+h_s/2*b1_, T+h_s/2*c1, s)
          a3, b3, c3, _ = rhs(z+h_s/2*a2, v1+h_s/2*b2, T+h_s/2*c2, s)
          a4, b4, c4, _ = rhs(z+h_s*a3, v1+h_s*b3, T+h_s*c3, s)
          z = np.where(alive, z + h_s/6*(a1+2*a2+2*a3+a4), z); v1 = np.where(alive, v1 + h_s/6*(b1_+2*b2+2*b3+b4), v1)
          T = np.where(alive, T + h_s/6*(c1+2*c2+2*c3+c4), T)
      print(" MC M0 grid: mean=%.2f sd=%.2f q=%s" % (tte0.mean(), tte0.std(), np.round(np.percentile(tte0, [5, 50, 95]), 2)))
      print(" MC M3      : mean=%.2f sd=%.2f q=%s ; mean ratio M3/M0=%.4f ; z_off mean=%.4f max=%.4f" % (
          np.nanmean(tte3), np.nanstd(tte3), np.round(np.nanpercentile(tte3, [5, 50, 95]), 2), np.nanmean(tte3/tte0), np.nanmean(zoff), np.nanmax(zoff)))
