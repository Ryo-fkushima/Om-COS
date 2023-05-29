# **Om-COS: Cation ordering simulator for omphacite**

**Om-COS** is a set of Python routines to simulate antiphase domain (APD) growth in metastable disordered omphacite. Based on a phenomenological phase-field approach with a Landau expansion of the excess free energy, these scripts allow to calculate ordered domain growth in a 2D/3D space. We can also simulate corresponding transmission electron microscope images for the calculated phase fields by setting a brightness curve against the degree of long range order. Om-COS should be a useful tool to properly interpret both the cation ordering process and temporal changes of APD size/morphology in natural omphacites from low-temperature eclogites.

## Requirement
Python 3.7

Modules required: matplotlib, numba, numpy, pandas, and scipy

## How to use
0. Install python3 to your computer. (The scripts were tested on: Python 3.8.2 (pip3), macOS Big Sur (11.6.4), MacBook Air (2020, 1.1 GHz Quad-Core Intel Core i5, 16 GB RAM)). 
Note: This program requires several modules: matplotlib, numba, numpy, pandas, and scipy. (You would be able to run the scripts in the default Anaconda distribution.)

1. Make a new directory for the calculation, and copy & paste all the scripts there.

2. Edit the contents of params.py and initials.py (optional) to set up calculation parameters.

3. Run pfcalc.py with some command-line keys (see below).

4. Run the other programs (pfvisual.py, genplot.py, pfsnap.py) according to your preference.

## Author
Ryo Fukushima (Tohoku University, Sendai, Japan)

## References



