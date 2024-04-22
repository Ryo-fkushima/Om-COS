# genplot.py
# version: 1.0.1 (Apr 22, 2024)
# author: Ryo Fukushima (Tohoku University) rpifukushima@gmail.com 
#
import OmpApd as apd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pandas as pd
import datetime

if len(sys.argv) > 3:
	print("Error: Add 0-2 values. ex) genplot.py apdsize fraction")
	sys.exit()

if os.path.exists("summary") == False:

	print("Error: Summary not found")
	sys.exit()

print("============================================================")
print("\n", "Om-COS v1.0.1  genplot.py")
print("\n", "Summary plot (", datetime.datetime.now(), ")\n")
print("============================================================")

	
	
Summary_list = []
Name = 0

while os.path.isfile("summary/summary_%s.csv" % str(Name)) == True:

		Summary_element = pd.read_csv("summary/summary_%s.csv" % str(Name))
		
		Summary_list.append(Summary_element)
	
		Name = Name + 1
		
		
Repeat = Name


if len(sys.argv) == 1: # default plots

	ForPlot = Summary_list[0]
	MeanOrder = ForPlot.iloc[:,2]
	MeanApdSize = ForPlot.iloc[:,4]
	pFraction = ForPlot.iloc[:,3]

	

	plt.figure()
	plt.xlim(0,4.0e-7)
	plt.ylim(0,1)
	plt.ylabel("P2/n fraction")
	plt.xlabel("Mean APD size, m")
    
	plt.plot(np.ma.masked_where(MeanApdSize == 0, MeanApdSize), np.ma.masked_where(pFraction == 0, pFraction), ".", color = "k")
    



	plt.figure()
	plt.xlim(0,3000)
	plt.ylim(0,1)
	plt.ylabel("Mean order")
	plt.xlabel("Time step")
    
	plt.plot(np.ma.masked_where(MeanOrder == 0, MeanOrder), ".", color = "k")
    


	plt.figure()
	plt.xlim(0,3000)
	plt.ylim(0,1)
	plt.ylabel("P2/n fraction")
	plt.xlabel("Time step")
    
	plt.plot(np.ma.masked_where(pFraction == 0, pFraction), ".", color = "k")
  
	# showing all plots
	plt.show()

	'''
	#%%  histogram

	apd.HistPlot2d(PhaseField, 1500)

	#%%  autocorrelation curves

	for i in np.arange(0,3001,100):

    	apd.AcCalcXPlot(PhaseField,100, i, 0)
	'''



if len(sys.argv) == 2: # with 1 command-line variable

	if sys.argv[1] == "help":
		print("\n time: real time\n order: mean order parameter\n fraction: P2/n fraction\n apdsize: mean APD size\n apdno: numbers of the measured APD length\n sd: standard deviation of the measured APD lengths\n wavelength: spatial wavelength of the phase field\n")
		sys.exit()

	plt.figure()
	plt.xlabel("time step")
	
	if sys.argv[1] == "time":
		plt.ylabel("time, yr")
		fory = 1

	if sys.argv[1] == "order":
		plt.ylabel("Mean order")
		fory = 2
	
	if sys.argv[1] == "fraction":
		plt.ylabel("P2/n fraction")
		fory = 3	
	
	if sys.argv[1] == "apdsize":
		plt.ylabel("Mean APD size, m")
		fory = 4

	if sys.argv[1] == "apdno":
		plt.ylabel("APD counts")
		fory = 5
	
	if sys.argv[1] == "sd":
		plt.ylabel("sd of APD size")
		fory = 6	
	
	if sys.argv[1] == "wavelength":
		plt.ylabel("Wavelength, m")
		plt.xscale("log")
		plt.yscale("log")
		fory = 7
		
	if (fory == 7) & (Summary_list[0].shape[1] < 8):
		print("Error: No wavelength data")
		sys.exit()
	
	for i in range(Repeat):

		ForPlot = Summary_list[i]
		yvalue = ForPlot.iloc[:,fory]

		plt.plot(np.ma.masked_where(yvalue == 0, yvalue), ".", color = "k")
 
	# showing all plots
	plt.show()


if len(sys.argv) == 3: # with 2 command-line variables

	plt.figure()
	
	if sys.argv[1] == "time":
		plt.xlabel("time, yr")
		forx = 1

	if sys.argv[1] == "order":
		plt.xlabel("Mean order")
		forx = 2
	
	if sys.argv[1] == "fraction":
		plt.xlabel("P2/n fraction")
		forx = 3	
	
	if sys.argv[1] == "apdsize":
		plt.xlabel("Mean APD size, m")
		forx = 4

	if sys.argv[1] == "apdno":
		plt.xlabel("APD counts")
		forx = 5
	
	if sys.argv[1] == "sd":
		plt.xlabel("sd of APD size")
		forx = 6	
	
	if sys.argv[1] == "wavelength":
		plt.xlabel("Wavelength, m")
		forx = 7
		
	if (forx == 7) & (Summary_list[0].shape[1] < 8):
		print("Error: No wavelength data")
		sys.exit()	
	
	if sys.argv[2] == "time":
		plt.ylabel("time, yr")
		fory = 1

	if sys.argv[2] == "order":
		plt.ylabel("Mean order")
		fory = 2
	
	if sys.argv[2] == "fraction":
		plt.ylabel("P2/n fraction")
		fory = 3	
	
	if sys.argv[2] == "apdsize":
		plt.ylabel("Mean APD size, m")
		fory = 4

	if sys.argv[2] == "apdno":
		plt.ylabel("APD counts")
		fory = 5
	
	if sys.argv[2] == "sd":
		plt.ylabel("sd of APD size")
		fory = 6	
	
	if sys.argv[2] == "wavelength":
		plt.ylabel("Wavelength, m")
		fory = 7
	
	if (fory == 7) & (Summary_list[0].shape[1] < 8):
		print("Error: No wavelength data")
		sys.exit()	

	for i in range(Repeat):

		ForPlot = Summary_list[i]
		xvalue = ForPlot.iloc[:,forx]
		yvalue = ForPlot.iloc[:,fory]

		plt.plot(np.ma.masked_where(xvalue == 0, xvalue), np.ma.masked_where(yvalue == 0, yvalue), ".", color = "k")
    
    # showing all plots
	plt.show()

