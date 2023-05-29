# **Om-COS: Cation ordering simulator for omphacite**

**Om-COS** is a set of Python routines to simulate antiphase domain (APD) growth in metastable disordered omphacite. 
Based on a phenomenological phase-field approach with a Landau expansion of the excess free energy, these scripts allow to calculate ordered domain growth in a 2D/3D space. 
We can also simulate corresponding transmission electron microscope images for the calculated phase fields by setting a brightness curve against the degree of long range order. 
This should be a useful tool to properly interpret both the cation ordering process and temporal changes of APD size/morphology in natural omphacites from low-temperature eclogites.

## Requirement

The scripts were tested on: Python 3.8.2 (pip3), macOS Big Sur (11.6.4), MacBook Air (2020, 1.1 GHz Quad-Core Intel Core i5, 16 GB RAM). 

Required modules: matplotlib, numba, numpy, pandas, and scipy

## How to use

1. Make a new directory for the calculation, and copy & paste all the scripts there.

2. Edit the contents of **params.py** and **initials.py** to set up calculation parameters.

3. Run **pfcalc.py**.

4. Run the other programs (**pfvisual.py**, **genplot.py**, **pfsnap.py**) according to your preference.

For more details, please refer to **README.txt**.

## Author

Ryo Fukushima (Tohoku University, Sendai, Japan)

## References

Allen, S. M., & Cahn, J. W. (1979). A microscopic theory for antiphase boundary motion and its application to antiphase domain coarsening. Acta Metallurgica, 27, 1085–1095.

Buscombe, D., Rubin, D. M., & Warrick, J. A. (2010). A universal approximation of grain size from images of noncohesive sediment. Journal of Geophysical Research, 115, F02015.

Cahn, J. W., & Hilliard, J. E. (1958). Free energy of a nonuniform system. I. Interfacial free energy. The Journal of Chemical Physics, 28(2), 258–267.

Carpenter, M. A., Domeneghetti, M. -C., & Tazzoli, V. (1990). Application of Landau theory to cation ordering in omphacite I: Equilibrium behaviour. European Journal of Mineralogy, 2, 7–18.

Holland, T. J. B., & Powell, R. (2011). An improved and extended internally consistent thermodynamic dataset for phases of petrological interest, involving a new equation of state for solids. Journal of Metamorphic Geology, 29, 333–383.
