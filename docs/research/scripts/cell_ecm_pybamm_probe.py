import os
import pybamm, numpy as np, json, warnings
warnings.filterwarnings("ignore")
out = {}
def ocv_table(name):
    pv = pybamm.ParameterValues(name)
    # Use electrode OCPs at stoichiometries given by esoh-like linear interpolation between 0% and 100% SOC
    model = pybamm.lithium_ion.SPM()
    sim = pybamm.Simulation(model, parameter_values=pv, experiment=pybamm.Experiment(["Discharge at C/50 until %s V" % pv["Lower voltage cut-off [V]"]]))
    sol = sim.solve(initial_soc=1.0)
    Q = sol["Discharge capacity [A.h]"].entries
    V = sol["Voltage [V]"].entries
    Qmax = Q[-1]
    soc = 1 - Q/Qmax
    grid = np.round(np.arange(0,1.0001,0.05),2)
    Vg = np.interp(grid, soc[::-1], V[::-1])
    return Qmax, dict(zip([float(g) for g in grid],[round(float(v),4) for v in Vg]))

def cap_vs_rate(name, crates, T_C=25.0, model_cls=pybamm.lithium_ion.DFN):
    res = {}
    for c in crates:
        pv = pybamm.ParameterValues(name)
        pv["Ambient temperature [K]"] = 273.15+T_C
        pv["Initial temperature [K]"] = 273.15+T_C
        try:
            exp = pybamm.Experiment([f"Discharge at {c}C until {pv['Lower voltage cut-off [V]']} V"])
            sim = pybamm.Simulation(model_cls(), parameter_values=pv, experiment=exp)
            sol = sim.solve(initial_soc=1.0)
            res[c] = round(float(sol["Discharge capacity [A.h]"].entries[-1]),4)
        except Exception as e:
            res[c] = "fail:"+str(e)[:80]
    return res

def pulse_R(name, soc0, T_C=25.0, crate=1.0, tpulse=10.0):
    pv = pybamm.ParameterValues(name)
    pv["Ambient temperature [K]"] = 273.15+T_C
    pv["Initial temperature [K]"] = 273.15+T_C
    Qn = pv["Nominal cell capacity [A.h]"]
    I = crate*Qn
    exp = pybamm.Experiment([f"Rest for 60 seconds", f"Discharge at {I} A for {tpulse} seconds or until 2.5 V", "Rest for 600 seconds"], period="0.1 seconds")
    sim = pybamm.Simulation(pybamm.lithium_ion.DFN(), parameter_values=pv, experiment=exp)
    sol = sim.solve(initial_soc=soc0)
    t = sol["Time [s]"].entries; V = sol["Voltage [V]"].entries
    V0 = V[np.searchsorted(t,59.9)]
    # instantaneous ~0.1 s after start
    i1 = np.searchsorted(t, 60.1); i10 = np.searchsorted(t, 60+tpulse-0.01)
    R_01 = (V0 - V[i1])/I; R_10 = (V0 - V[i10])/I
    # relaxation: voltage at end of rest
    Vend = V[-1]
    return dict(I=round(I,4), R_0p1s_mOhm=round(1000*R_01,2), R_10s_mOhm=round(1000*R_10,2), V_rest_before=round(float(V0),4), V_after_600s_rest=round(float(Vend),4))

for name in ["Ai2020","Marquis2019"]:
    d = {}
    try:
        Q, tab = ocv_table(name); d["C50_capacity_Ah"]=round(float(Q),4); d["OCV_C50_SPM"]=tab
    except Exception as e: d["ocv_err"]=str(e)[:200]
    d["cap_vs_rate_25C"] = cap_vs_rate(name,[0.1,0.5,1,2])
    for T in [0,-10,-20]:
        d[f"cap_1C_{T}C"] = cap_vs_rate(name,[0.5,1],T_C=T)
    d["pulse"] = {}
    for T in [25,0,-10,-20]:
        for s in [0.8,0.5,0.2,0.1]:
            try:
                d["pulse"][f"T{T}_soc{s}"] = pulse_R(name,s,T_C=T)
            except Exception as e:
                d["pulse"][f"T{T}_soc{s}"] = "fail:"+str(e)[:100]
    out[name]=d
    print(name, json.dumps(d, indent=1), flush=True)
json.dump(out, open(os.path.join(os.environ.get("CELL_ECM_OUT", "."), "pybamm_probe.json"),"w"), indent=1)
