# pfhist.py
# version: 1.0.1 (Apr 22, 2024)
# author: Ryo Fukushima (Tohoku University) rpifukushima@gmail.com
#
import OmpApd as apd
import params
import importlib
import numpy as np
import sys
import os
import datetime

if len(sys.argv) != 5:
	print("Error: Add 4 values for the histogram bin numbers, gif speed, output name, and Y/N for the raw-data output. \n ex.\n pfhist.py 10 100 hist N\n pfhist.py 20 100 hist N\n pfhist.py 40 100 hist Y")
	sys.exit()

importlib.reload(params)


if os.path.isfile("PhaseField.npy") == True:

	PhaseField = np.load("PhaseField.npy")
	

elif os.path.isfile("PhaseField3d.npy") == True:

	PhaseField3d = np.load("PhaseField3d.npy")
	PhaseField = PhaseField3d[:,:,params.VisualParams["slicefor3d"],:]
	

else:
	
	print("Error: Phase field not found")
	sys.exit()

print("============================================================")
print("\n", "Om-COS v1.0.1  pfhist.py")
print("\n", "Histogram calculation (", datetime.datetime.now(), ")\n")
print("============================================================")


apd.HistAnime2d(PhaseField, BinsInput = int(sys.argv[1]), speed = int(sys.argv[2]), filename = sys.argv[3], **params.VisualParams)

print("%s.gif generated" % sys.argv[3])

if (sys.argv[4] == "Y") or (sys.argv[4] == "y"): 

	AbsPF_allTime = apd.PFAbsRavel2d_jit(PhaseField)

	np.savetxt("histdata_%s.txt" % sys.argv[3], AbsPF_allTime[:,params.VisualParams["plotNo"]])

	print("histdata_%s.txt generated" % sys.argv[3])




