# Pressurant tank blowdown through a choked orifice

Self initiated study. Every input is ASSUMED, none comes from a measurement or an employer.

## What it does
`blowdown.py` models a 50 L nitrogen tank at 200 bar and 293.15 K discharging through a 6 mm orifice (discharge coefficient 0.8)
to ambient, as a simple launch site pneumatic line. The flow is choked throughout. Two limits are solved as ODEs with SciPy:
- isothermal: the gas stays at 293.15 K (heat input from the tank wall keeps up)
- adiabatic: no heat exchange, the gas in the tank cools as it expands

## Results (results.csv)
| model | time to 50 bar s | gas mass used kg | end temperature K |
|---|---|---|---|
| isothermal | 15.17 | 8.62 | 293.15 |
| adiabatic | 11.99 | 7.22 | 197.28 |

## Assumptions and limits
- Ideal gas nitrogen, constant specific heat ratio 1.4. Real gas effects at 200 bar are ignored.
- Choked flow at every instant, constant discharge coefficient, no line losses, no valve dynamics, no regulator.
- Real tanks sit between the two limits. The study does not model wall heat transfer.
- Tank volume, pressures, orifice size and discharge coefficient are typed in as assumptions.
- Not validated against a measurement.

## Run
`pip install numpy scipy pandas matplotlib` then `python blowdown.py`. It writes results.csv and blowdown.png.
