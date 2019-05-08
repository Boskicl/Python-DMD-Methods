import numpy as np
from numpy import diag, power
from scipy.linalg import expm, sinm, cosm, pinv, eig
import matplotlib.pyplot as plt
import pandas as pd

###########################. Import Data from Excel Sheet. ###################################
#DataCompanionMatrix.xlsx is the first set of power grid data used in Paper (Yoshi,Mezić)
#DataCompanionMatrixUCTE.xlsx is the second set of power grid data from Paper (Yoshi,Mezić)

df = pd.read_excel('DataCompanionMatrix.xlsx', header=None) #Insert your own excel file in between '__', if header is desired, delete header=None.
data = np.array(df,dtype=float)					    #Make sure the data file will be a numpy array.

###########################. FUNCTION DEFINE. #################################################
def Arnoldi(data):
	# Get dimensions of Data Matrix
	m = data.shape[0] #rows of matrix
	n = data.shape[1] #columns of matrix

	### Re-define data matrix into x->{1,..,m-1} and y->{last column}
	x = data[0:-1,:]
	y = data[-1,:]

	### Create A matrix and xx matrix which is used to find c_j values
	A = np.dot(x,x.T)
	xx = np.dot(x,y.T)
	Cj_values = np.dot(pinv(A),xx)

	### Building Companion Matrix
	CompanionMatrix = np.zeros((m-1,m-1))

	for ii in range(0,m-2):
		CompanionMatrix[ii+1,ii] = 1

    ### Fill last row of Companion Matrix with cj values
	CompanionMatrix[:,m-2] = Cj_values

	### Compute empirical Ritz eigenvalues/vectors --> Same as Koopman Eigenvalues/vectors
	eigV,eigW = eig(CompanionMatrix)

	### Koopman Eigs/Empirical Ritz values
	Keigs = eigV

	### Define Vandemonde Matrix
	VandemondeMatrix = np.ones((m-1,m-1),dtype=complex)
	for ii in range(1,m-1):
		VandemondeMatrix[:,ii] = np.power(Keigs,(ii))

	### Koopman Eigenvectors/Emperical Ritz Vector
	Kmodes = np.dot(x.T,pinv(VandemondeMatrix))

	########################### Rest of Computation from Arnoldi Paper ########################### 
	# args = np.zeros((m-1,m-1))
	# modes = 0 * Kmodes
	# for ii in range(0,m-2):
	# 	modes[:,ii] = abs(Kmodes[:,ii])
	# 	args[:,ii] = np.angle(Kmodes[:,ii])

	# mode1_mode5 = np.zeros((m-1,2))
	# phi1_phi5 = np.zeros((m-1,2))

	# for ii in range(0,m-1):
	# 	mode1_mode5[ii,:] = modes[ii,[0,4]]
	# 	phi1_phi5[ii,:] = args[ii,[0,4]]
	# x_pos = np.arange(len(mode1_mode5))

	# # plt.subplot(1)
	# plt.subplot(111)
	# w = 0.3 
	# plt.bar(x_pos,mode1_mode5[:,0],width=w, color='b',align='center',label='Mode 1')
	# plt.bar(x_pos + w,mode1_mode5[:,1],width=w, color='r',align='center',label='Mode 5')
	# plt.autoscale(tight=True)
	# plt.ylabel('Modulus [MW]')
	# plt.ylim((0,10))
	# plt.legend()
	# plt.show()

	return Keigs, Kmodes

###########################. Call Function . #################################################

Keigs, Kmodes = Arnoldi(data)

###########################. Plots .   #################################################

## Plot of Signals 
plt.figure(0)
plt.plot(data[:,1:]) #In the data file from the paper, the first signal is not used thus we do 1:
plt.xlabel('Time [h]')
plt.ylabel('Power Flow [MW]')
plt.show()

## Koopman Eigs Plotted on Unit Circle 
t = np.linspace(0,2*np.pi,101)
plt.figure(0)
plt.plot(np.cos(t),np.sin(t),'--')
plt.plot(Keigs.real,Keigs.imag,'ro')
plt.xlabel(r'Real $\ lambda$')
plt.ylabel(r'Imaginary $\ lambda$')
plt.axes().set_aspect('equal')
m = max(max(abs(Keigs.real)),max(abs(Keigs.imag)))
plt.xlim(-1.1*m,1.1*m)
plt.ylim(-1.1*m,1.1*m)
plt.show()

############################################################################
#
#				Author: Ljuboslav Boskic
#				UCSB Graduate Student- Mezic Group
#				Contact: lboskic@ucsb.edu
#
#
############################################################################

