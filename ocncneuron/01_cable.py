from neuron import h, gui
import matplotlib.pyplot as plt
import os

cable = h.Section()

cable.L = 100  # Length
cable.diam = 1  # Diameter
# cable