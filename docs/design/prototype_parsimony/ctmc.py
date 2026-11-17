import numpy as np, sys
sys.path.insert(0, '.')
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from proto import state_power, Eu0
acts = ['idle', 'social', 'video', 'game', 'nav', 'vcall', 'camera', 'audio']
def profile(kind):
    if kind == 'typical':
        lam = 3.0; p = dict(social=.55, video=.12, game=.05, nav=.04, vcall=.04, camera=.08, audio=.12)
        d = dict(social=5, video=15, game=20, nav=20, vcall=10, camera=2, audio=30)
    elif kind == 'heavy':
        lam = 4.0; p = dict(social=.45, video=.20, game=.15, nav=.05, vcall=.05, camera=.05, audio=.05)
        d = dict(social=8, video=25, game=30, nav=25, vcall=15, camera=3, audio=30)
    return lam, p, {k: v/60 for k, v in d.items()}  # hours
def gen_activity(lam, p, d):
    n = len(acts); L = np.zeros((n, n))
    for j, a in enumerate(acts[1:], 1):
        L[0, j] = lam*p[a]; L[j, 0] = 1/d[a]
    np.fill_diagonal(L, -L.sum(1)); return L
def stationary(L):
    n = L.shape[0]; A = np.vstack([L.T, np.ones(n)]); b = np.zeros(n+1); b[-1] = 1
    return np.linalg.lstsq(A, b, rcond=None)[0]
def asym(L, Pv):
    pi = stationary(L); n = len(pi); Pi = np.outer(np.ones(n), pi)
    D = np.linalg.inv(Pi - L) - Pi
    Pbar = pi @ Pv; pt = Pv - Pbar
    s2 = 2*np.sum(pi*pt*(D @ pt)); h = D @ Pv  # bias: E[int_0^t P] = Pbar t + h_i + o(1)
    return pi, Pbar, s2, h
# connectivity chain: wifi (mean 3h) <-> 4g (mean 2h)
Lc = np.array([[-1/3, 1/3], [1/2, -1/2]]); rats = ['wifi', '4g']
def product(La, Lc):
    na, nc = La.shape[0], Lc.shape[0]
    return np.kron(La, np.eye(nc)) + np.kron(np.eye(na), Lc)
Pw = {(a, r): state_power(a, r)[0] for a in acts for r in rats}
Pw[('game', 'wifi')] = min(Pw[('game', 'wifi')], 2.9); Pw[('game', '4g')] = min(Pw[('game', '4g')], 2.95)  # thermal-cap steady approx
if __name__ == "__main__":
  for kind in ["typical", "heavy"]:
      lam, p, d = profile(kind)
      La = gen_activity(lam, p, d)
      # hub closed form, connectivity frozen at mixture powers
      L = product(La, Lc)
      Pv = np.array([Pw[(a, r)] for a in acts for r in rats])
      pi, Pbar, s2, h = asym(L, Pv)
      Eu = Eu0  # Wh
      mu = Eu/Pbar; var = s2*Eu/Pbar**3; cv = np.sqrt(var)/mu
      pia = pi.reshape(len(acts), 2).sum(1)
      share = {a: (pi.reshape(len(acts), 2)[i] @ np.array([Pw[(a, r)] for r in rats]))/Pbar for i, a in enumerate(acts)}
      print("\n== %s ==  Pbar=%.4f W  sigma2=%.4f W^2 h  mu_TTE=%.2f h  sd=%.2f h  CV=%.3f" % (kind, Pbar, s2, mu, np.sqrt(var), cv))
      print(" pi_activity:", {a: round(x, 4) for a, x in zip(acts, pia)})
      print(" energy share:", {a: round(x, 4) for a, x in share.items()})
      print(" tau_c = s2/(2 Var_pi P) = %.3f h; CV_P=%.3f" % (s2/(2*np.sum(pi*(Pv-Pbar)**2)), np.sqrt(np.sum(pi*(Pv-Pbar)**2))/Pbar))
      print(" bias h (idle,wifi)=%.4f Wh ; (game,4g)=%.4f Wh" % (h[0], h[acts.index('game')*2+1]))
      # hub closed form for activity only at fixed rat = 4g
      Pa = np.array([Pw[(a, '4g')] for a in acts])
      pi_a, Pbar_a, s2_a, _ = asym(La, Pa)
      EL = 1/lam + sum(p[a]*d[a] for a in acts[1:])
      s2_hub = 2*sum(p[a]*d[a]**2*(Pa[i+1]-Pbar_a)**2 for i, a in enumerate(acts[1:]))/EL
      print(" activity-only 4g: Pbar=%.4f  s2(deviation matrix)=%.5f  s2(hub closed form)=%.5f" % (Pbar_a, s2_a, s2_hub))
      # Monte Carlo M0 on product chain
      rng = np.random.default_rng(1)
      N = 20000; ttes = np.empty(N)
      rates = -np.diag(L); Pjump = (L - np.diag(np.diag(L)))/rates[:, None]
      cum = np.cumsum(Pjump, 1)
      start = 0  # idle, wifi
      for k in range(N):
          s = start; E = Eu; t = 0.0
          while True:
              tau = rng.exponential(1/rates[s]); e = Pv[s]*tau
              if e >= E: t += E/Pv[s]; break
              E -= e; t += tau
              s = np.searchsorted(cum[s], rng.random()*cum[s, -1], side='right')
          ttes[k] = t
      q = np.percentile(ttes, [5, 50, 95])
      # inverse Gaussian quantiles
      from scipy.stats import invgauss
      mu_b = (Eu - h[start])/Pbar; lamIG = Eu**2/s2
      ig = invgauss(mu=mu_b/lamIG, scale=lamIG)
      print(" MC M0: mean=%.2f sd=%.2f CV=%.3f  q5/50/95=%s ; IG(bias-corr) mean=%.2f q=%s" % (
          ttes.mean(), ttes.std(), ttes.std()/ttes.mean(), np.round(q, 2), mu_b, np.round(ig.ppf([.05, .5, .95]), 2)))
