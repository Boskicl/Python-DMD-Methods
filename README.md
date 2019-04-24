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

[first equation][https://latex.codecogs.com/gif.latex?x%27%20%3D%20T%28x%29]

  ${(x0,y0),(x1,y1),…(xn,yn)}$

in which each xi and yi is a column vector of size m. We now define two m×n matrices:

  $X=[x_o x_1 ... x_n]$,  $Y=[y_o y_1 ... y_n]$

If we define an operator $A$ as

  $A=YX^†$

where $X^†$ is the pseudo-inverse of $X$, then the Dynamic Mode Decomposition of the pair $(X,Y)$ is given by the eigendecomposition of $A$. That is, the DMD modes and eigenvalues are eigenvectors and eigenvalues of $A$.

The definition above, from Tu et al., is known as the exact DMD. It is currently the most general definition and can be applied to any dataset that satisfies the given requirements. In this post, we are mostly interested in the cases where $A$ satisfies (perhaps approximately) the equation $y_i=Ax_i$ for all $i$. Or, more precisely:

$Y=AX$

Clearly, $X$ is a set of inputs vectors and $Y$ is the corresponding set of output vectors. This particular interpretation of the DMD is extremely powerful, as it provides a convenient method for analyzing (and predicting) dynamical systems for which the governing equations are unknown. More on dynamical systems shortly.

There are a number of theorems that go along with this definition of the DMD 2. One of the more useful theorems states that $Y=AX$ exactly if and only if $X$ and $Y$ are linearly consistent (i.e., whenever $Xv=0$ for some vector $v$, then $Yv=0$ too). Linear consistency is relatively straightforward to test, as we shall see. That being said, linear consistency is not a mandatory prerequisite for using the DMD. Even if the DMD solution for A doesn’t exactly satisfy the equation $Y=AX$, it is still a least-squares solution, minimizing error in an $L_2$ sense.

## Companion Matrix Algorithm (Arnoldi)
1) Define constants $c_j$ such that for vector $r$ satifying $r \bot {X_o X_1 ... X_{n-2}}$

$r = X_{n-1} - \sum_{j=-0}^{n-2} c_j X_j$ 

2) Define the companion matrix C as:

$C :=   \begin{bmatrix}
    0 0 ... 0 c_o \
    1 0 ... 0 c_1 \
    0 1 ... 0 c_2 \
    . . ... . . \
    0 0 ... 1 c_{n-2}    
  \end{bmatrix}$
