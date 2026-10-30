"""Blowdown of a nitrogen pressurant tank through a choked orifice, a simple launch site pneumatic line model.
Isothermal and adiabatic limits are compared. All inputs are ASSUMED (self initiated study)."""
import numpy as np, pandas as pd
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = 296.8           # J/kgK nitrogen
GAM = 1.4
V = 0.050           # m3 tank (assumed)
T0 = 293.15         # K
P0 = 200e5          # Pa start (assumed)
PEND = 50e5         # Pa end of useful pressure (assumed)
D = 6e-3            # m orifice diameter (assumed)
CD = 0.8            # discharge coefficient (assumed)
AO = CD * np.pi * D ** 2 / 4
PAMB = 1.013e5
CH = np.sqrt(GAM / R) * (2 / (GAM + 1)) ** ((GAM + 1) / (2 * (GAM - 1)))   # choked flow factor

def mdot(p, T):
    return AO * p * CH / np.sqrt(T)

def run(mode):
    m0 = P0 * V / (R * T0)
    def f(t, y):
        m, T = y
        p = m * R * T / V
        md = mdot(p, T)
        if mode == "isothermal":
            return [-md, 0.0]
        # adiabatic: energy balance of the tank contents, T falls as gas leaves
        dT = -(GAM - 1) * T * md / m
        return [-md, dT]
    ev = lambda t, y: y[0] * R * y[1] / V - PEND
    ev.terminal = True
    sol = solve_ivp(f, [0, 600], [m0, T0], events=ev, max_step=0.05, rtol=1e-8)
    t = sol.t; m, T = sol.y
    return t, m * R * T / V, T

rows = []
out = {}
for mode in ("isothermal", "adiabatic"):
    t, p, T = run(mode)
    out[mode] = (t, p)
    rows.append(dict(model=mode, time_to_50bar_s=t[-1], gas_mass_used_kg=(P0 * V / (R * T0)) - p[-1] * V / (R * T[-1]), end_temperature_K=T[-1]))
df = pd.DataFrame(rows)
df.to_csv("results.csv", index=False)
print(df.round(3).to_string())
plt.figure(figsize=(5, 3.4))
for k, (t, p) in out.items():
    plt.plot(t, p / 1e5, label=k)
plt.xlabel("time (s)"); plt.ylabel("tank pressure (bar)"); plt.legend(); plt.title("Tank blowdown through a 6 mm orifice (assumed data)", fontsize=9)
plt.tight_layout(); plt.savefig("blowdown.png", dpi=130)
