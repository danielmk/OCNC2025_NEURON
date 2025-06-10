from neuron import h
from neuron.units import ms, mV
import matplotlib.pyplot as plt
import os
import numpy as np

magee_file = "main.hoc"

h.load_file(magee_file)

ps = h.PlotShape(False)
ps.plot(plt)

