from neuron import h
from neuron.units import ms, mV
import matplotlib.pyplot as plt
import os
import numpy as np
import sys


sheath_length = 15
node_length = 2
n_sheaths = 100
node_conductance = 0.001  # S/cm2
sheath_conductance = 0.000000001  # S/cm2

first_node = h.Section()
nodes = [first_node]
sheaths = []
for i in range(n_sheaths):
    current_sheath = h.Section()
    nodes[i].connect(current_sheath)
    next_node = h.Section()
    current_sheath.connect(next_node)
    nodes.append(next_node)
    sheaths.append(current_sheath)

for n in nodes:
    n.insert('pas')
    n.insert('hh')
    n.L = node_length
    n.diam = 1
    n.nseg = 10
    for seg in n:
        seg.pas.g = node_conductance

for s in sheaths:
    s.insert('pas')  
    s.nseg = 10
    s.L = sheath_length
    s.diam = 1
    s.cm = 0.000001
    for seg in s:
        seg.pas.g = sheath_conductance

control_node = h.Section()
control_length = len(nodes) * node_length + len(sheaths) * sheath_length
control_segs = len(nodes) * 10 + len(sheaths) * 10
control_node.L = control_length
control_node.diam = 1
control_node.insert('pas')
control_node.insert('hh')
control_node.nseg = control_segs
for seg in control_node:
    seg.pas.g = node_conductance
    
clamp = h.IClamp(nodes[0](0.0))
clamp.delay = 20 * ms
clamp.dur = 0.1 * ms
clamp.amp = 2

control_clamp = h.IClamp(control_node(0.0))
control_clamp.delay = 20 * ms
control_clamp.dur = 0.1 * ms
control_clamp.amp = 2

voltage_end = h.Vector()
voltage_end.record(nodes[-1](1.0)._ref_v)

voltage_start = h.Vector()
voltage_start.record(nodes[0](0.0)._ref_v)

control_voltage_end = h.Vector()
control_voltage_end.record(control_node(1.0)._ref_v)

control_voltage_start = h.Vector()
control_voltage_start.record(control_node(0.0)._ref_v)

h.finitialize(-65)  # Initialize voltage to -65 mV
tstop = 100 * ms
h.dt = 0.01

while h.t < tstop:
    h.fadvance()

plt.figure()
plt.plot(voltage_start, color='k')
plt.plot(control_voltage_start, color='r')
plt.plot(voltage_end, color='k', linestyle='--')
plt.plot(control_voltage_end, color='r', linestyle='--')
plt.legend(["Sheath Start", "Control Start", "Sheath End", "Control End"])


