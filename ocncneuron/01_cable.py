from neuron import h
from neuron.units import ms, mV
import matplotlib.pyplot as plt
import os
import numpy as np

# Create two section
passive_cable = h.Section()

passive_cable.L = 1000  # Length in um
passive_cable.diam = 1  # Diameter in um
passive_cable.Ra = 100  # 100 MOhm axial resistance
passive_cable.insert(h.pas)
passive_cable.nseg = 100

active_cable = h.Section()

active_cable.L = 1000  # Length in um
active_cable.diam = 1  # Diameter in um
active_cable.Ra = 100  # 100 MOhm axial resistance
active_cable.insert(h.pas)
active_cable.insert(h.hh)  # Insert Hodgkin-Huxley mechanism into 
active_cable.nseg = 100

# Insert a current clamp process on one end of each cable
passive_iclamp = h.IClamp(passive_cable(0.0))
active_iclamp = h.IClamp(active_cable(0.0))

for ic in [passive_iclamp, active_iclamp]:
    ic.delay = 10 * ms  # Inject after 200 ms
    ic.dur = 0.1 * ms  # Inject for 10 ms
    ic.amp = 1  # Inject 500 nA peak current

# Measure the voltage at multiple places along each cable
mloc = np.arange(0, 1.1, 0.2)
passive_vectors = [h.Vector().record(passive_cable(loc)._ref_v) for loc in mloc]
active_vectors = [h.Vector().record(active_cable(loc)._ref_v) for loc in mloc]

# Set initial values and run the simulation
h.finitialize(-65)  # Initialize voltage to -65 mV
tstop = 30
h.dt = 0.01

while h.t < tstop:
    h.fadvance()

# Plot the results    
passive_arrays = np.array(passive_vectors)
active_arrays = np.array(active_vectors)

fig, ax = plt.subplots(2,1)
t = np.arange(0, tstop+0.02, h.dt)
for idx, loc in enumerate(mloc):
    ax[0].plot(t, passive_arrays[idx, :])
    ax[1].plot(t, active_arrays[idx, :])

ax[1].set_xlabel("Time (ms)")
for a in ax:
    a.set_ylabel("Voltage (mV)")
    a.legend(mloc * 1000)

ax[0].set_title("Passive Cable")
ax[1].set_title("Active Cable")

