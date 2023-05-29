# pfsnap.py
# version: 1.0.0 (May 29, 2023)
# author: Ryo Fukushima (Tohoku University) rpifukushima@gmail.com
#
import numpy as np
import sys
import os
import datetime

if len(sys.argv) != 3:
	print("Error: Add 2 values for dimensions and slice No.\n ex.\n pfsnap.py 2d 100\n pfsnap.py 2dtxt 100\n pfsnap.py 3d 1000")
	sys.exit()


if (os.path.isfile("PhaseField.npy") == False) and (os.path.isfile("PhaseField3d.npy") == False) :

	print("Error: Phase field not found")
	sys.exit()


if os.path.isfile("PhaseField.npy") == True:

	PhaseField = np.load("PhaseField.npy")
	
if os.path.isfile("PhaseField3d.npy") == True:

	PhaseField3d = np.load("PhaseField3d.npy")


SliceNo = int(sys.argv[2])

if sys.argv[1] == "2d":
	
	np.save("pfsnap2d_%s.npy" % sys.argv[2], PhaseField[:,:,SliceNo])
	
if sys.argv[1] == "2dtxt":
	
	np.savetxt("pfsnap2d_%s.txt" % sys.argv[2], PhaseField[:,:,SliceNo])

	
if sys.argv[1] == "3d":

	np.save("pfsnap3d_%s.npy" % sys.argv[2], PhaseField3d[:,:,:,SliceNo])


print("file saved")
