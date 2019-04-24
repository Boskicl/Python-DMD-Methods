# Python-DMD-Methods
The following codes will be for Dynamic Mode Decomposition(DMD)/Koopman Methods. Codes will have different methods to compute Koopman Eigenvalues and Modes. They will be verified by reproducing various papers from academia. The goal is to have these numerical methods written in Python (currently use Matlab). I currently am a researcher in Mezić Group at UCSB in the Department of Mechanical Engineering. 


# Current codes:
1) Arnoldi-Like Method Algorithm.  
  Varified Paper: Nonlinear Koopman Modes and Power System Stability Assessment Without Models- Yoshihiko Susuki and Igor Mezić
2) Hankle-DMD Algorithm. (Working on)
3) exact-DMD Algorithm. (Working on)

# Code Use:
The codes have their own data files for them specifically data from the papers. Some of the data files are available to the public while some are not. If they are allowed, I will add them into a new folder. Most of the codes I will upload will only require you to put in your own data file (excel file or mat file or whatever you use) and run the code. 

## Formal Definition (Theory)
Koopman Operator Theory is an alternative formulation of dynamical system theory which provides a versitle framework for data-driven methods of high-dimensional nonlinear systems. The theory orginated in the 1930s through the work of Koopman and Von Neuman. Work done in the previous few years has proved the spectrial decompostion, introducing the idea of Koopman Modes. This theory has led to data-driven methods to approximate the Koopman operator spectrum and mode.

In discrete time setting, if:

<img src="http://www.sciweavers.org/tex2img.php?eq=x%27%3DT%28x%29&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="x'=T(x)" width="78" height="19" />

is a discrete time dynamical system where <img src="http://www.sciweavers.org/tex2img.php?eq=x%20%5Cin%20M&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="x \in M" width="50" height="15" /> and <img src="http://www.sciweavers.org/tex2img.php?eq=T%3A%20M%20%5Cto%20M&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="T: M \to M" width="86" height="15" />, then the associated Koopman Operator U is defined as:

<img src="http://www.sciweavers.org/tex2img.php?eq=Uf%28x%29%20%3D%20f%20%20%5Ccirc%20T%28x%29&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="Uf(x) = f  \circ T(x)" width="131" height="19" />

We call <img src="http://www.sciweavers.org/tex2img.php?eq=%20%5Cphi%20%20%3A%20M%20%5Cto%20%5Cmathbb%7BC%7D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt=" \phi  : M \to \mathbb{C}" width="211" height="19" /> an eigenfunction of U associated with <img src="http://www.sciweavers.org/tex2img.php?eq=%5Clambda%20%5Cin%20%5Cmathbb%7BC%7D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="\lambda \in \mathbb{C}" width="212" height="15" /> then
<img src="http://www.sciweavers.org/tex2img.php?eq=U%5Cphi%20%3D%20%5Clambda%20%5Cphi&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="U\phi = \lambda \phi" width="72" height="19" />
And in continous time,
<img src="http://www.sciweavers.org/tex2img.php?eq=U%5Et%5Cphi%20%3D%20exp%28%5Clambda%20t%29%5Cphi&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="U^t\phi = exp(\lambda t)\phi" width="133" height="22" />

