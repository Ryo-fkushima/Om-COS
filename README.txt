Om-COS: Cation ordering simulator for omphacite
version: 1.0.1 (Apr 22, 2024)
author: Ryo Fukushima (Tohoku University) rpifukushima@gmail.com

===== List of the scripts =====

[What you can basically edit]

params.py: parameters for calculations
initials.py: used to set an initial phase-field

[Other scripts (just to run from your console)]

pfcalc.py: script for phase-field calculations
pfvisual.py: script for phase-field visualization
pfhist.py: script for histogram generation (power spectrum of the order parameter)
genplot.py: script for plotting characteristic values of the calculated phase-field
pfsnap.py: script for taking a snapshot of the calculated phase-field

[Function list]

OmpApd.py

===== How to use =====

0. Install python3 to your computer. (The scripts were tested on: Python 3.8.2 (pip3), macOS Big Sur (11.6.4), MacBook Air (2020, 1.1 GHz Quad-Core Intel Core i5, 16 GB RAM)). 
Note: This program requires several modules: matplotlib, numba, numpy, pandas, and scipy. (You would be able to run the scripts in the default Anaconda distribution.)

1. Make a new directory for the calculation, and copy & paste all the scripts there.
2. Edit the contents of params.py and initials.py (optional) to set up calculation parameters.
3. Run pfcalc.py with some command-line keys (see below).
4. Run the other programs (pfvisual.py, pfhist.py, genplot.py, pfsnap.py) according to your preference.


===== Descriptions of the calculation scripts =====

[pfcalc.py]

This calculates the phase field and makes outputs as a summary.
"2d" option enables to calculate in a 2d space (fastest).
"3d" option (= normal 3d mode) enables to calculate in a 3d space, saving 3d matrices for all time steps (slowest and much memory is required). 
"3d_cs" option (= cross-section mode) enables to calculate in a 3d space (making a 2d slice for each time step), without saving any 3d matrices.
"omit" option enables to omit summary calculations.
"omit_wl" option enables to omit wavelength calculation, but generate a summary for the other values.
A "summary" folder will be automatically generated. (If already present, this function overwrite the summary folder.)
Calculated phase-field data will be stored in "PhaseField.npy" or "PhaseField3d.npy".

ex. 
python3 pfcalc.py 2d
python3 pfcalc.py 3d
python3 pfcalc.py 3d_cs
python3 pfcalc.py 2d omit
python3 pfcalc.py 2d omit_wl

[pfvisual.py]

This generates a gif movie of the phase-field time evolution.
First command-line key designates the plot type.
(-1/0: order parameter; 1: DF-like image; 2: Binary image)
Second command-line key designates the animation speed.
Third command-line key will be the name of the exported gif file.
You can visualize the phase field even if you did not calculate the summary. 
(i.e. "omit" options in pfcalc.py)

ex.
python3 pfvisual.py -1 100 animation
python3 pfvisual.py -1 400 animation
python3 pfvisual.py 0 100 animation
python3 pfvisual.py 1 100 animation
python3 pfvisual.py 2 100 animation

[pfhist.py]

This generates a gif movie of the time-evolution of histograms of order-parameter values.
First command-line key designates the number of the bins.
Second command-line key designates the animation speed.
Third command-line key will be the name of the exported gif file.
You can get raw data for the histogram as .txt file by setting the fourth command-line key as "Y".
You can get histograms even if you did not calculate the summary. 
(i.e. "omit" options in pfcalc.py)

ex.
python3 pfhist.py 10 100 histogram Y
python3 pfhist.py 40 100 histogram N



[genplot.py]

This generates several plots of characteristic values of the phase-field based on the calculated summary.
The command-line keys designates parameters for the plot.
("time": real time; "order": mean order parameter; "fraction": P2/n fraction; "apdsize": mean APD size;
 "apdno": numbers of the measured APD length; "sd": standard deviation of the measured APD lengths; 
 "wavelength": spatial wavelength of the phase field)
Note: You can also see a parameter list above with a single command-line key "help".
If you set no command-line keys, some default plots would be generated automatically.
If you set only the first command-line key, the designated value would be plotted to the time step.
(With a default setting, you can see a log-log plot when you choose "wavelength" as a command-line key.)
If you set both the first and second command-line key, the second value (Y) would be plotted to the first value (X). 
When you set the "Repeat" value as >1 for pfcalc.py, all the calculated data will be simulanaously plotted.
This does not work when you chose the "omit" option in pfcalc.py.


ex.
python3 genplot.py
python3 genplot.py apdsize
python3 genplot.py wavelength
python3 genplot.py apdsize fraction


[pfsnap.py]

This is useful to save the calculated phase-field at a specific time step.
The command-line keys designates the dimension type and which time step you want to save the phase-field with.
("2d": 2d phase-field will be saved in a .npy file; "2dtxt": 2d phase-field will be saved in a .txt file; 
"3d": 3d phase-field will be saved in a .npy file)
You may use the output file as the next initial phase-field by editing initials.py properly.
This script works even if you did not calculate the summary. (i.e. "omit" options in pfcalc.py)

ex.
python3 pfsnap.py 2d 100
python3 pfsnap.py 2dtxt 100
python3 pfsnap.py 3d 1000


===== General comments =====

> For the present, all the summary calculations, phase-field visualization, and histogram generations will be done
for selected 2d phase-fields, even if you choose the "3d" option in pfcalc.py.
> "mob", "premob", and "actEmob" values are important only when you consider how long each time step is. 
These values do not affect APD patterns.
> The "repeat" value may be larger than 1 when you want to check effects of the initial randomness with genplot.py. 
(note that only the final phase-field data will be saved.)
> The file size of PhaseField.npy and PhaseField3d.npy may be so large (> a few GB).
> initials.py would be read only when you set "ini" as -1.
> "ApdMeasureWidth", "MaxOffset", and "slicefor3d/slicefor3d_cs" values should be smaller than the "nx", "ny", and "nz" values, respectively.
> When you calculate "wavelength" with a too small area for the APD size, it sometimes fails. Then choose the "omit_wl" option instead.


===== Examples for parameter settings =====

[When you calculate 2d phase-field with initial fluctuation of Gaussian or uniform distributions]

> Make sure that "nx", "ny", "dx", "dy", "periodi", and "periodj" values are correct, and that "ini" value is 0/1 according to your preferrence.
"nz", "dz", "periodk", and "slicefor3d/slicefor3d_cs" values are irrelevant to such 2d calculations.


[When you calculate 2d phase-field with specific initial fluctuation patterns]

> Make sure that "nx", "ny", "dx", "dy", "periodi", and "periodj" values are correct, and that "ini" and "inifield" values are -1 and "initials.initial_PF_2d".
"nz", "dz", "periodk", and "slicefor3d/slicefor3d_cs" values are irrelevant to such 2d calculations.
In initials.py, you should confirm that the matrix size of "initial_PF_2d" (and temperature when you use apd.phi_m(temp), which is the positive phi value with the minimum Gibbs energy) is consistent with the input parameters written in params.py.


[When you calculate 3d phase-field with initial fluctuation of Gaussian or uniform distributions]

> Make sure that "nx", "ny", "nz", "dx", "dy", "dz", "periodi", "periodj", and "periodk" values are correct, and that "ini" value is 0/1 according to your preferrence.


[When you calculate 3d phase-field with specific initial fluctuation patterns:]

> Make sure that "nx", "ny", "nz", "dx", "dy", "dz", "periodi", "periodj", and "periodk" values are correct, and that "ini" and "inifield" values are -1 and "initials.initial_PF_3d", respectively.
In initials.py, you should confirm that the matrix size of "initial_PF_3d" (and temperature when you use apd.phi_m(temp)) is consistent with the input parameters written in params.py.











