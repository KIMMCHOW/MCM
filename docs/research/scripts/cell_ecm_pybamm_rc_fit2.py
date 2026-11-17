import os
import pybamm, numpy as np, json, warnings
from scipy.optimize import least_squares
warnings.filterwarnings("ignore")
S = os.environ.get("CELL_ECM_OUT", ".") + "/"
ocvfit = json.load(open(S+"pybamm_ocv_fit.json"))
Qc50 = {"Ai2020":2.4651, "Marquis2019":0.8724}
out = {}
def run(name, soc0, T_C, crate=1.0, tp=60.0, trest=3600.0):
    pv = pybamm.ParameterValues(name)
    pv["Ambient temperature [K]"] = 273.15+T_C; pv["Initial temperature [K]"] = 273.15+T_C
    I = crate*pv["Nominal cell capacity [A.h]"]
    exp = pybamm.Experiment([pybamm.step.rest(60, period=1), pybamm.step.current(I, duration=tp, period=0.02), pybamm.step.rest(trest, period=1)])
    sim = pybamm.Simulation(pybamm.lithium_ion.DFN(), parameter_values=pv, experiment=exp)
    sol = sim.solve(initial_soc=soc0)
    return I, sol["Time [s]"].entries, sol["Voltage [V]"].entries
def fit(name, soc0, I, t, V, tp=60.0):
    tab = ocvfit[name]["OCV_table"]; zz = np.array([float(k) for k in tab]); vv = np.array(list(tab.values()))
    # local OCV slope (V per unit SOC) around soc0
    slope = (np.interp(soc0+0.05, zz, vv) - np.interp(soc0-0.05, zz, vv))/0.10
    Q = Qc50[name]*3600.0
    t0 = 60.0; t1 = t0+tp
    Vocv0 = V[np.searchsorted(t, t0)-1]
    m = t > t0+1e-6; tt = t[m]; VV = V[m]
    q = I*np.clip(tt-t0, 0, tp)  # charge removed [C]
    Vocv = Vocv0 - slope*q/Q
    def model(p, n_rc):
        R0 = p[0]; out = Vocv - 0.0
        on = tt <= t1
        vrc = np.zeros_like(tt)
        for j in range(n_rc):
            R, tau = p[1+2*j], p[2+2*j]
            v = np.empty_like(tt)
            v[on] = I*R*(1-np.exp(-(tt[on]-t0)/tau))
            vend = I*R*(1-np.exp(-tp/tau))
            v[~on] = vend*np.exp(-(tt[~on]-t1)/tau)
            vrc += v
        return Vocv - I*R0*on - vrc
    res = {}
    for n_rc, p0, lb, ub in [(1,[0.04,0.01,10],[0,0,0.01],[2,2,3000]), (2,[0.04,0.01,1,0.01,100],[0,0,0.01,0,5],[2,2,30,2,3000])]:
        r = least_squares(lambda p: model(p,n_rc)-VV, p0, bounds=(lb,ub))
        e = model(r.x,n_rc)-VV
        d = {"R0_mOhm": round(1000*r.x[0],2)}
        for j in range(n_rc):
            d[f"R{j+1}_mOhm"] = round(1000*r.x[1+2*j],2); d[f"tau{j+1}_s"] = round(float(r.x[2+2*j]),2)
        d["rmse_mV"] = round(1000*float(np.sqrt(np.mean(e**2))),3)
        res[f"{n_rc}RC"] = d
    res["I_A"] = round(I,4); res["dOCVdz_V"] = round(float(slope),3)
    return res
for name, temps in [("Ai2020",[25,0]),("Marquis2019",[25,0,-10,-20])]:
    out[name] = {}
    for T in temps:
        for s in [0.8,0.5,0.2]:
            try:
                I,t,V = run(name, s, T)
                out[name][f"T{T}_soc{s}"] = fit(name, s, I, t, V)
            except Exception as ex:
                out[name][f"T{T}_soc{s}"] = "fail "+str(ex)[:80]
            print(name, T, s, out[name][f"T{T}_soc{s}"], flush=True)
json.dump(out, open(S+"pybamm_rc_fit2.json","w"), indent=1)
