# pfcalc.py
# version: 1.0.0 (May 29, 2023)
# author: Ryo Fukushima (Tohoku University) rpifukushima@gmail.com
#
import OmpApd as apd
import numpy as np
import pandas as pd
import params
import importlib
import sys
import os
import shutil
import datetime

if len(sys.argv) != 2 and len(sys.argv) != 3:
	print("Error: Add 1 (or 2) values.\n ex.\n pfcalc.py 2d\n pfcalc.py 3d\n pfcalc.py 3d_cs\n pfcalc.py 2d omit\n pfcalc.py 2d omit_wl")
	sys.exit()


print("\n", "Phase-field calculation (", datetime.datetime.now(), ")\n")

importlib.reload(params)

Repeat = params.CalcParams["repeat"] 
print("loop for %s time(s)" % str(Repeat))

if os.path.exists("summary") == True:
	
	shutil.rmtree("summary")
		
	os.mkdir("summary")
		
else:

	os.mkdir("summary")
		
print("Summary updated")

Summary_element = []

if (len(sys.argv) == 2) or (len(sys.argv) == 3 and sys.argv[2] != "omit") or (len(sys.argv) == 3 and sys.argv[2] != "omit_wl"):

	Summary_element = pd.DataFrame(index=np.arange(0,(params.CalcParams["nsteps"] + 1),1), columns=[])
	

for i in range(Repeat):
	
	print("\n",str(i + 1))
	
	Name = str(i)
    
	if sys.argv[1] == "2d":
    	
		PhaseField = apd.CalcApd2d(**params.CalcParams)
	
	if sys.argv[1] == "3d":
	
		PhaseField3d = apd.CalcApd3d(**params.CalcParams)
		PhaseField = PhaseField3d[:,:,params.VisualParams["slicefor3d"],:]
		print("slicefor3d = %d" % params.VisualParams["slicefor3d"])
		
	if sys.argv[1] == "3d_cs":
	
		PhaseField = apd.CalcApd3d_CS(**params.CalcParams)
		print("z = %d cross sections" % params.CalcParams["slicefor3d_cs"])
		
	if len(sys.argv) == 2:
    
		Summary_element = apd.QuantEvalApd2d(PhaseField, Name = Name, WLswitch = 1, **params.ParamsSummary)

		
	if len(sys.argv) == 3:
	
		if sys.argv[2] == "omit_wl":
    
			Summary_element = apd.QuantEvalApd2d(PhaseField, Name = Name, WLswitch = 0, **params.ParamsSummary)
			
			
		if (sys.argv[2] != "omit") & (sys.argv[2] != "omit_wl"):
		
			Summary_element = apd.QuantEvalApd2d(PhaseField, Name = Name, WLswitch = 1, **params.ParamsSummary)
			

	Summary_element.to_csv("summary/summary_%s.csv" % str(i))



if sys.argv[1] == "2d":

	np.save("PhaseField.npy", PhaseField)
	print("2d phase-field generated")
	
	if os.path.isfile("PhaseField3d.npy") == True:

		os.remove("PhaseField3d.npy")
	
	
if sys.argv[1] == "3d":

	np.save("PhaseField3d.npy", PhaseField3d)
	print("3d phase-field generated")
	
	if os.path.isfile("PhaseField.npy") == True:

		os.remove("PhaseField.npy")
	

if sys.argv[1] == "3d_cs":

	np.save("PhaseField.npy", PhaseField)
	print("2d phase-field (cross section) generated")
	
	if os.path.isfile("PhaseField3d.npy") == True:

		os.remove("PhaseField3d.npy")
    


