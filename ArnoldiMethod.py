import numpy as np
from numpy import diag, power
from scipy.linalg import expm, sinm, cosm
import matplotlib.pyplot as plt
import pandas as pd

###########################. Import Data from Excel Sheet. ###################################

df = pd.read_excel('DataCompanionMatrix.xlsx', header=None) #Insert your own excel file in between '__', if header is desired, delete header=None.
data = np.array(df)					    #Make sure the data file will be a numpy array.

###########################. FUNCTION DEFINE. #################################################
def Arnoldi(data):
	# Get dimensions of Data Matrix
	m = data.shape[0] #rows of matrix
	n = data.shape[1] #columns of matrix
	# Re-define data matrix into x->{1,..,m-1} and y->{last column}
	x = data[0:-1,:]
	y = data[-1,:]
	# Create A matrix and xx matrix which is used to find c_j values
	A = np.dot(x,np.transpose(x))
	xx = np.dot(x,np.transpose(y))
	Cj_values = np.dot(np.linalg.pinv(A),xx)
	# Building Companion Matrix
	CompanionMatrix = np.zeros((m-1,m-1))

	for i in range(0,m-2):
		CompanionMatrix[i+1,i] = 1
    #Fill last row of Companion Matrix with cj values
	CompanionMatrix[:,m-2] = Cj_values
	# Compute empirical Ritz eigenvalues/vectors --> Same as Koopman Eigenvalues/vectors
	eigV,eigW = np.linalg.eig(CompanionMatrix)

	#Koopman Eigs
	Keigs = eigV
	#Koopman Modes
	Kmodes = np.matmul(eigW,x)
	return Keigs, Kmodes
#
###########################. Call Function . #################################################

Keigs, Kmodes = Arnoldi(data)

###########################. Plots .   #################################################

## Plot of Signals 
plt.figure(0)
plt.plot(data[:,1:]) #In the data file from the paper, the first signal is not used thus we do 2:
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



