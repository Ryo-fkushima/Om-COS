# pfvisual.py
# version: 1.0.0 (May 29, 2023)
# author: Ryo Fukushima (Tohoku University) rpifukushima@gmail.com
#
import OmpApd as apd
import params
import importlib
import numpy as np
import sys
import os
import datetime
import pprint

if len(sys.argv) != 4:
	print("Error: Add 3 values for the plot type, gif speed, and output name.\n ex.\n pfvisual.py -1 100 animation\n pfvisual.py 0 100 animation\n pfvisual.py 1 100 animation\n pfvisual.py 2 100 animation")
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

print("\n", "Phase-field visualization (", datetime.datetime.now(), ")\n")


def Saturation2Order(**args):

	saturation = args["saturation"]
	x1, yfactor, A = args["phiadjust"], args["Ifactor"], args["Icoef"]
	y1 = yfactor * A * (1 - x1) * (1 - x1) / (1 - yfactor)
	
	intensity = saturation * (A * (1 - x1) * (1 - x1) + y1)
	
	if intensity > y1:
		return ((intensity - y1) / A)**(0.5) + x1
		
	else:
		return (((y1 - intensity) / y1)**0.5) * x1 * (-1) + x1


if int(sys.argv[1]) == -1:

    print("Phase-field visualization (normalized by maximum phi in each step)")
    
if int(sys.argv[1]) == 0:

    print("Phase-field visualization (standard)")

if int(sys.argv[1]) == 1:

    print("DF-image like visualization\nSaturated at phi = ", "{:.4f}".format(Saturation2Order(**params.VisualParams)))

if int(sys.argv[1]) == 2:

    print("Binary visualization\nThreshold: phi = ", "{:.4f}".format(params.VisualParams["binthres"]))



#apd.PlotApd2d(PhaseField, plotmode = int(sys.argv[1]), **params.VisualParams)
apd.ApdAnime2d(PhaseField, speed = int(sys.argv[2]), plotmode = int(sys.argv[1]), filename = sys.argv[3], **params.VisualParams)

print("%s.gif was generated with the following parameters:" % sys.argv[3])

TrueVps = params.VisualParams
del TrueVps["ApdMeasureWidth"]
del TrueVps["MaxOffset"]
del TrueVps["plotNo"]

pprint.pprint(TrueVps)




