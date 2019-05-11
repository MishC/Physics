import numpy

import pylab

import scipy

import scipy.interpolate

def ReadMultiColFile(filename, numArrs):

	f = open(filename, "r")

	arrCounter = 0

	outPut = [[]]* numArrs

	for i in range(numArrs):
		outPut[i] = []

	for line in f.readlines():

		if line[0] == "#":
			continue

		else:
			syms = line.split(" ")
			for char in range(len(syms)):
				flo = -1.0
				try:
					flo = float(syms[char])
				except:
					continue
				if flo != -1.0:
					outPut[arrCounter].append(flo)
					arrCounter += 1
					if arrCounter == numArrs:
						arrCounter = 0

	return numpy.array(outPut).transpose()

def GetRho(Egrid):
	"""Returns density of states"""

	N = Egrid.shape[0]

	rhoGrid = numpy.zeros((N))

	rhoGrid[0] = Egrid[1] - Egrid[0]

	rhoGrid[N-1] = Egrid[N-1] - Egrid[N - 2]

	for i in range(1, N - 1):

		rhoGrid[i] = ( Egrid[i + 1] - Egrid[i - 1] ) / 2.0


	return rhoGrid

def MakedPdE(energies, amplitudes, Egrid):

	sliceSize = N
	amps = abs(amplitudes)**2


	dPdE = numpy.zeros((numE))

	counter = 0

	#Loop over l
	for l in range(lMax + 1):

		#Extracting slice
		Eslize = energies[int(counter) : int(counter + N)]
		Cslize =     amps[int(counter) : int(counter + N)]
		#Checking energies
		for i in range(int(N)):
			#If energy negative, set Cslice to zero
			if Eslize[i] < 0:
				Cslize[i] = 0.0
				pass

		rhoGrid = GetRho(Eslize)

		#Interpolate
		intObj = scipy.interpolate.interp1d(Eslize, Cslize * 1. / rhoGrid, kind = 'cubic', fill_value = 0.0, bounds_error = False)

		#Add to energy distribution
		dPdE = dPdE + intObj(Egrid)

		counter = int(counter + N)

	return dPdE

Emin = 0.0001
Emax = 3.4
numE = 3000

lMax = 5

filename_E = "amplitudes_B200_N600_l25_om350_c15_E45_dip_nom.txt"
filename_C = "amplitudes_B200_N600_l25_om350_c15_E45_dip_nom.txt"



print("Running Emin: ", Emin)
print("Running Emax: ", Emax)
print("Running ", numE, " grid points in energy ")

print("Expecting lmax = ", lMax)


print("Looking for energies in: ", filename_E)
print("Looking for amplitudes in: ", filename_C)

innArr1 = ReadMultiColFile(filename_E, 3)

E = innArr1[:,0]

Emin = max(Emin, numpy.min(E))
Emax = min(Emax, numpy.max(E))

print("Corrected Emin, Emax: ", Emin, Emax)

Egrid = numpy.linspace(Emin, Emax, numE)

innArr = ReadMultiColFile(filename_C, 3)


numC = innArr.shape[0]
print("Found : ", numC , " amplitudes")

N = numC / (lMax  + 1)

print("Found ", N, " amplitudes per l")

print("Loaded : ", numC, " amplitudes. ")


if numC != (lMax + 1) * N:
	print("Wrong Dimension !!!!!")

E = innArr[:,0]

amps = innArr[:,1] + 1.0j * innArr[:,2]

print("Calling energy distribution module....")

dPdE = MakedPdE(E, amps, Egrid)


dE = Egrid[1] - Egrid[0]


print("Integral: ",numpy.sum(dPdE) * dE)#(Egrid[1] - Egrid[0]))

print(Egrid.shape, dPdE.shape)
#print dPdE[:10]



pylab.figure()
pylab.plot(Egrid, dPdE, label = "dP/dE")

#pylab.xlim(0.0, 35.0)
#pylab.ylim(-4, 1.)
pylab.legend()
pylab.savefig("dpdefig.pdf")
pylab.show()
