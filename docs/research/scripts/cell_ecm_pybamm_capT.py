import os
import pybamm, numpy as np, json, warnings
warnings.filterwarnings("ignore")
S = os.environ.get("CELL_ECM_OUT", ".") + "/"
out = {}
for name in ["Ai2020","Marquis2019"]:
    out[name] = {}
    for T in [25,10,0,-10,-20]:
        for c in [0.1,0.2,0.5]:
            pv = pybamm.ParameterValues(name)
            pv["Ambient temperature [K]"] = 273.15+T; pv["Initial temperature [K]"] = 273.15+T
            try:
                exp = pybamm.Experiment([f"Discharge at {c}C until {pv['Lower voltage cut-off [V]']} V"])
                sol = pybamm.Simulation(pybamm.lithium_ion.DFN(), parameter_values=pv, experiment=exp).solve(initial_soc=1.0)
                q = float(sol["Discharge capacity [A.h]"].entries[-1]); e = float(sol["Discharge energy [W.h]"].entries[-1]) if "Discharge energy [W.h]" in sol.all_models[0].variables else None
                out[name][f"T{T}_C{c}"] = {"Q_Ah": round(q,4), "E_Wh": (round(e,4) if e else None)}
            except Exception as ex:
                out[name][f"T{T}_C{c}"] = "fail "+str(ex)[:60]
            print(name,T,c,out[name][f"T{T}_C{c}"],flush=True)
json.dump(out, open(S+"pybamm_capT.json","w"), indent=1)
