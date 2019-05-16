import numpy as np
import matplotlib.pyplot as plt
from numpy import dot, multiply, diag, power, zeros, linspace, meshgrid
from numpy import pi, exp, sin, cos, cosh, tanh, real, imag
from numpy.linalg import inv, eig, pinv
from scipy.linalg import svd, svdvals, norm



###########################. Create Dynamics. ###################################

#define time and space domains
x = linspace(-10, 10, 100)
t = linspace(0, 6*pi, 80)
dt = t[2] - t[1]
Xm,Tm = meshgrid(x, t)

# create three spatiotemporal patterns
f1 = multiply(20-0.2*power(Xm, 2), exp((2.3j)*Tm))
f2 = multiply(Xm, exp(0.6j*Tm))
f3 = multiply(5*multiply(1/cosh(Xm/2), tanh(Xm/2)), 2*exp((0.1+2.8j)*Tm))

# combine signals and make data matrix
D = (f1 + f2 + f3).T

###########################. Define Function . ###############################################

def ExactDMD(data,relTol):
    # Define Data matrix X and Y
    #   X = [Do D1 ... Dn-1]    &    Y = [D1 D2 ... Dn]
    m = data.shape[0] #rows of data matrix
    n = data.shape[1] #columns of data matrix

    X = data[:,:-1] # take X_o to X_n-1
    Y = data[:,1:]  # take X_1 to X_n

    ################ Compute SVD of data matrix X
    U, Sdiag, Vh = svd(X, False)

    #Plot Signular values
    plt.figure()
    plt.plot(Sdiag,'ro')
    plt.title('Singular Values Plot')


    S = zeros((Sdiag.shape[0], Sdiag.shape[0]))  # Create S matrix with zeros based on Diag of S
    np.fill_diagonal(S, Sdiag)  # Fill diagonal of S matrix with the nonzero values
    V = Vh.conj().T  # Create V matrix, we are given Vh which is conjugate transpose of V (convert back to V)
    r = np.count_nonzero(np.diag(S) > S[0, 0] * relTol)  # Find the diminsion of the subspace we are using- truncation
    U = U[:, 0:r]
    S = S[0:r,0:r]
    V = V[:,0:r]

    #################### Create Atilde
    Atilde = dot(dot(dot(U.conj().T, Y), V), inv(S))

    ## Get Eigenvalues and Vectors from Atilde (eigen pairs)
    Evalue, Wvector = eig(Atilde)

    ############################Koopman Eigenvalues
    Keigs = Evalue

    #Plot Koopam Eigenvalues on unit circle
    tt = linspace(0,2*np.pi,101)
    plt.figure()
    plt.plot(np.cos(tt),np.sin(tt),'--')
    plt.plot(Keigs.real,Keigs.imag,'ro')
    plt.title('DMD Eigenvalues')
    plt.xlabel(r'Real $\ lambda$')
    plt.ylabel(r'Imaginary $\ lambda$')
    plt.axes().set_aspect('equal')
    m = max(max(abs(Keigs.real)),max(abs(Keigs.imag)))
    plt.xlim(-1.1*m,1.1*m)
    plt.ylim(-1.1*m,1.1*m)



    ##################### Koopman Eigenmodes (Dynamics Modes)
    Kmodes = dot(dot(dot(Y, V), inv(S)), Wvector)

    #Plot Koopman Modes
    plt.figure()
    plt.plot(Kmodes[:,0],'b-',label='1st Mode')
    plt.plot(Kmodes[:,1],'r-',label='2nd Mode')
    plt.plot(Kmodes[:,2],'k-',label='3rd Mode')
    plt.legend(loc='upper left')
    plt.title('Exact DMD Modes')


    ###############################Reconstruction
    b = dot(pinv(Kmodes), X[:,0])
    Xdmd = zeros([r, len(t)], dtype='complex')
    for i,_t in enumerate(t):
        Xdmd[:,i] = multiply(power(Keigs, _t/dt), b)

    #Plot time dynamics of Modes
    plt.figure()
    plt.plot(Xdmd[0,:],'r-')
    plt.title('1st Mode Time Dynamics')
 

    plt.figure()
    plt.plot(Xdmd[1,:],'b-')
    plt.title('2nd Mode Time Dynamics')

    plt.figure()
    plt.plot(Xdmd[2,:],'g-')
    plt.title('3rd Mode Time Dynamics')
  
    #Reconstruct full spectrum of data
    Xreconstruct = dot(Kmodes, Xdmd)
    #np.allclose(data.T, Xreconstruct)

    plt.figure()
    plt.plot(data[0,:],'ro',label='Data')
    plt.plot(Xreconstruct[0,:],'k-',label='Reconstruction')
    plt.title('Data vs Reconstruction')
    plt.legend()
    plt.show()

    ######### Realtive Error
    Q = dot(dot(dot(Kmodes,diag(Keigs)),pinv(Kmodes)),X)
    relativeError = norm(Q-Y,'fro')/norm(Y,'fro')

    return Keigs, Kmodes, Xdmd, relativeError

relTol = 10**-6
Keigs, Kmodes, Xdmd, relativeError = ExactDMD(D,relTol)
print('Amount of Koopman Eigenvalues :',Keigs.shape)
print('Amount of Koopman Modes :',Kmodes.shape[1])
print(relativeError)
