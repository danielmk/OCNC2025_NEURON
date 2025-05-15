from setuptools import setup

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='ocncneuron',
    version='0.1.0',
    description='Code for the NEURON tutorial at OCNC',
    long_description=readme,
    author='Daniel Müller-Komorowska',
    author_email='danielmuellermsc@gmail.com',
    url='https://github.com/danielmk/OCNC2025_NEURON',
    license=license,
    packages=['ocncneuron'],
    install_requires=[
          'neuron',
          'numpy',
          'pandas',
          'matplotlib',
          'spyder',],)
