import os
import pybamm, numpy as np, json, warnings
from scipy.optimize import least_squares
warnings.filterwarnings("ignore")
res = {}
for name in ["Ai2020", "Marquis2019"]:
    pv = pybamm.ParameterValues(name)
    param = pybamm.LithiumIonParameters()
    Vmin = pv["Lower voltage cut-off [V]"]; Vmax = pv["Upper voltage cut-off [V]"]
    x0, x100, y100, y0 = pybamm.lithium_ion.get_min_max_stoichiometries(pv)
    z = np.linspace(0, 1, 201)
    x = x0 + z*(x100 - x0); y = y0 - z*(y0 - y100)
    T = 298.15
    Un = np.array([pv.evaluate(param.n.prim.U(pybamm.Scalar(xi), pybamm.Scalar(T))) for xi in x]).astype(float).ravel()
    Up = np.array([pv.evaluate(param.p.prim.U(pybamm.Scalar(yi), pybamm.Scalar(T))) for yi in y]).astype(float).ravel()
    ocv = Up - Un
    d = {"x0":float(x0),"x100":float(x100),"y100":float(y100),"y0":float(y0),"Vmin":Vmin,"Vmax":Vmax}
    grid = np.round(np.arange(0,1.0001,0.05),2)
    d["OCV_table"] = {f"{g:.2f}": round(float(np.interp(g, z, ocv)),4) for g in grid}
    # entropic (only Ai2020 physically measured)
    if name == "Ai2020":
        dUn = np.array([pv.evaluate(pybamm.Scalar(0)+pv["Negative electrode OCP entropic change [V.K-1]"](pybamm.Scalar(xi), pybamm.Scalar(T)) if pv["Negative electrode OCP entropic change [V.K-1]"].__code__.co_argcount==2 else pv["Negative electrode OCP entropic change [V.K-1]"](pybamm.Scalar(xi))) for xi in x]).astype(float).ravel()
        dUp = np.array([pv.evaluate(pybamm.Scalar(0)+pv["Positive electrode OCP entropic change [V.K-1]"](pybamm.Scalar(yi), pybamm.Scalar(T)) if pv["Positive electrode OCP entropic change [V.K-1]"].__code__.co_argcount==2 else pv["Positive electrode OCP entropic change [V.K-1]"](pybamm.Scalar(yi))) for yi in y]).astype(float).ravel()
        dU = dUp - dUn
        d["dUdT_cell_mV_per_K"] = {f"{g:.2f}": round(1000*float(np.interp(g, z, dU)),4) for g in grid}
        d["dUdT_validity_note"] = "Ai2020 fits valid 0.43<sto<0.9936 per PyBaMM docstring"
    # Fit Plett combined model on z in [eps, 1-eps] with scaling s'=(1-2e)z+e, e=0.175 (Ahmed et al.) and also raw z in [0.02,0.98]
    def design(s):
        return np.vstack([np.ones_like(s), 1/s, s, np.log(s), np.log(1-s)]).T
    for tag, eps in [("combined_scaled_eps0.175", 0.175), ("combined_raw_z", None)]:
        if eps is None:
            m = (z>=0.02)&(z<=0.98); s = z[m]; V = ocv[m]
        else:
            s = (1-2*eps)*z + eps; V = ocv; m = np.ones_like(z,bool)
        A = design(s); k, *_ = np.linalg.lstsq(A, V, rcond=None)
        rms = float(np.sqrt(np.mean((A@k - V)**2)))*1000; mx = float(np.max(np.abs(A@k - V)))*1000
        d[tag] = {"K0..K4": [float(f"{v:.6g}") for v in k], "rmse_mV": round(rms,2), "max_err_mV": round(mx,2), "fit_range": "z in [0,1] scaled" if eps else "z in [0.02,0.98]"}
    # Fit Chen & Rincon-Mora form: V = -b0 exp(-b1 z) + a0 + a1 z + a2 z^2 + a3 z^3
    def f(p, zz): return -p[0]*np.exp(-p[1]*zz) + p[2] + p[3]*zz + p[4]*zz**2 + p[5]*zz**3
    p0 = [1.0, 35, 3.6, 0.3, -0.1, 0.3]
    r = least_squares(lambda p: f(p, z) - ocv, p0, bounds=([0,1,2,-5,-5,-5],[5,500,5,5,5,5]))
    e = f(r.x, z) - ocv
    d["chen_rm_form"] = {"b0,b1,a0,a1,a2,a3": [float(f"{v:.5g}") for v in r.x], "rmse_mV": round(float(np.sqrt(np.mean(e**2)))*1000,2), "max_err_mV": round(float(np.max(np.abs(e)))*1000,2)}
    # dOCV/dz at mid
    res[name] = d
print(json.dumps(res, indent=1))
json.dump(res, open(os.path.join(os.environ.get("CELL_ECM_OUT", "."), "pybamm_ocv_fit.json"),"w"), indent=1)
