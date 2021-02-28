"""A module to obtain the numerical solution of parabolic PDEs.

Various numerical methods are included in this module that can be used
to obtain the solution of the 1D or 2D model equation.

The diffusion equation is a parabolic partial differential equation.

A typical example of a diffusion equation is the unsteady heat
conduction equation, which is expressed in 1D as:

    du/dt = alpha(d2u/dx2)

    where,
                u: measurable quanity
            alpha: a constant physical coefficient

"""
#   ***********************************************************************
#
#   FILE         parabolicsolvers.py
#
#   AUTHOR       Dr. Vishal Sharma
#
#   VERSION      1.0.0-alpha4
#
#   WEBSITE      https://github.com/vxsharma-14/project-NAnPack
#
#   NAnPack Learner's Edition is distributed under the MIT License.
#
#   Copyright (c) 2020 Vishal Sharma
#
#   Permission is hereby granted, free of charge, to any person
#   obtaining a copy of this software and associated documentation
#   files (the "Software"), to deal in the Software without restriction,
#   including without limitation the rights to use, copy, modify, merge,
#   publish, distribute, sublicense, and/or sell copies of the Software,
#   and to permit persons to whom the Software is furnished to do so,
#   subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be
#   included in all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
#   OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
#   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
#   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
#   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#   SOFTWARE.
#
#   You should have received a copy of the MIT License along with
#   NAnPack Learner's Edition.
#
#   ***********************************************************************

from .backend.exceptions import DimensionError
from .tridiagonal import TridiagonalSolver


def FTCS(Uo, diffX, diffY=None):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the explicit Forward Time/Central Space method
    to obtain the solution of the 1D or 2D diffusion equation.

    Call signature:
        FTCS(Uo, diffX, diffY)

    Parameters
    ----------
    Uo : 1D or 2D array (depending on the domain dimensions)
        The dependent variable from time level (n) within the domain.
    diffX : float
        Diffusion number for x-component of the parabolic/diffusion
        equation.
    diffY : float, Default=None for 1-D applications
        Diffusion number for y-component of the parabolic/diffusion
        equation.
    Returns
    -------
    U : 1D or 2D array (depending on the domain dimensions)
        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    shapeU = Uo.shape  # Obtain Dimension
    U = Uo.copy()  # Initialize U

    if len(shapeU) == 1:
        U[1:-1] = (
            Uo[1:-1] + diffX*(Uo[2:] - 2.0*Uo[1:-1] + Uo[0:-2])
            )

    elif len(shapeU) == 2:
        U[1:-1, 1:-1] = (
            Uo[1:-1, 1:-1]
            + diffX*(Uo[2:, 1:-1] - 2.0*Uo[1:-1, 1:-1] + Uo[0:-2, 1:-1])
            + diffY*(Uo[1:-1, 2:] - 2.0*Uo[1:-1, 1:-1] + Uo[1:-1, 0:-2])
            )

    return U


def DuFortFrankel(Uo, Uo2, diffX, diffY=None):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the explicit DuFort-Frankel method
    to obtain the solution of the 1D or 2D diffusion equation.

    Call signature:
        DuFortFrankel(Uo, Uo2, diffX, diffY)

    Parameters
    ----------
    Uo : 1D or 2D array (depending on the domain dimensions)
        The dependent variable from time level (n) within the domain.
    Uo2 : 1D or 2D array
        The dependent variable at time level (n-1).
    diffX : float
        Diffusion number for x-component of the parabolic/diffusion
        equation.
    diffY : float, Default=None for 1-D applications
        Diffusion number for y-component of the parabolic/diffusion
        equation.

    Returns
    -------
    U : 1D or 2D array (depending on the domain dimensions)
        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    shapeU = Uo.shape  # Obtain Dimension
    U = Uo.copy()  # Initialize U

    if len(shapeU) == 1:
        U[1:-1] = (
            ((1.0 - 2.0*diffX)*Uo2[1:-1] + 2.0*diffX*(Uo[2:] + Uo[0:-2]))
            / (1.0 + 2.0*diffX)
            )

    elif len(shapeU) == 2:
        U[1:-1, 1:-1] = (
            ((1.0 - 2.0*diffX - 2.0*diffY)*Uo2[1:-1, 1:-1]
             + 2.0*diffX*(Uo[2:, 1:-1] + Uo[0:-2, 1:-1])
             + 2.0*diffY*(Uo[1:-1, 2:] + Uo[1:-1, 0:-2]))
            / (1.0 + 2.0*diffX + 2.0*diffY)
             )

    return U


def Laasonen(Uo, diffX):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the implicit Laasonen method
    to obtain the solution of the 1D diffusion equation.
    Note: Use only for 1D applications as this implicit formulation is not
    efficient for 2D. Use function ADI() instead for 2D implicit method.

    Call signature:
        Laasonen(Uo, diffX)

    Parameters
    ----------
    Uo : 1D array
        The dependent variable from time level (n) within the domain.
    diffX : float
        Diffusion number for x-component of the parabolic/diffusion
        equation.

    Returns
    -------
    U : 1D array
        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    shapeU = Uo.shape  # Obtain Dimension

    if len(shapeU) == 2:
        print("Laasonen method is inefficient for 2D applications.\
 \nUse Alternating Direction Implicit method instead.")
        raise DimensionError("2D", "diffusion", "Laasonen")
    iMax, = shapeU
    U = Uo.copy()  # Initialize U
    A = [-diffX for i in range(iMax)]
    B = [1.0 + 2.0*diffX for i in range(iMax)]
    C = [-diffX for i in range(iMax)]
    D = U
    UU = Uo.copy()

    U = TridiagonalSolver(iMax, A, B, C, D, UU)

    return U


def CrankNicolson(Uo, diffX):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the implicit CrankNicolson method
    to obtain the solution of the 1D diffusion equation.
    Note: Use only for 1D applications as this implicit formulation is not
    efficient for 2D. Use function ADI() instead for 2D implicit method.

    Call signature:
        CrankNicolson(Uo, diffX)

    Parameters
    ----------
    Uo : 1D array
        The dependent variable from time level (n) within the domain.
    diffX : float
        Diffusion number for x-component of the parabolic/diffusion
        equation.

    Returns
    -------
    U : 1D array
        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    shapeU = Uo.shape  # Obtain Dimension

    if len(shapeU) == 2:
        print("Crank-Nicolson method is inefficient for 2D applications.\
 \nUse Alternating Direction Implicit method instead.")
        raise DimensionError("2D", "diffusion", "Crank-Nicolson")
    iMax, = shapeU
    U = Uo.copy()  # Initialize U
    dd = 0.5*diffX
    A = [-dd for i in range(iMax)]
    B = [(1.0 + 2.0*dd) for i in range(iMax)]
    C = [-dd for i in range(iMax)]
    D = [0 for i in range(iMax)]
    D[1:-1] = dd*Uo[2:] + (1.0 - 2.0*dd)*Uo[1:-1] + dd*Uo[0:-2]
    UU = Uo.copy()

    U = TridiagonalSolver(iMax, A, B, C, D, UU)

    return U


def ADI(Uo, diffX, diffY):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the Alternating Direction Implicit method
    to obtain the solution of the 2D diffusion equation.
    Note: This formulation solves only 2D model equation.

    Call signature:
        ADI(Uo, diffX, diffY)

    Parameters
    ----------
    Uo : 2D array
        The dependent variable from time level (n) within the domain.
    diffX : float
        Diffusion number for x-component of the parabolic/diffusion
        equation.
    diffY : float
        Diffusion number for y-component of the parabolic/diffusion
        equation.

    Returns
    -------
    U : 2D array
        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    shapeU = Uo.shape  # Obtain Dimension

    if len(shapeU) == 1:
        print("Alternating Direction Implicit method is for 2D applicat\
 ions only.")
        raise DimensionError("1D", "diffusion", "ADI")
    iMax, jMax = shapeU
    U = Uo.copy()  # Initialize U
    Uhalf = Uo.copy()  # Uhalf is U at time level (n + 1/2)

    d1 = 0.5*diffX
    d2 = 0.5*diffY

    # *****************************************************************
    # This block of codes solves for U
    # at time level n + 1/2 (i.e. = Uhalf) along constant j line
    # Eq. 5.24 using Tridiagonal system Appendix B in
    # Hoffmann CFD Vol.1
    # *****************************************************************
    A = [-d1 for i in range(iMax)]
    B = [(1.0 + 2.0*d1) for i in range(iMax)]
    C = [-d1 for i in range(iMax)]
    D = [0 for i in range(iMax)]
    UU = [0 for i in range(iMax)]
    for j in range(1, jMax-1):
        UU[0] = Uo[0][j]
        UU[-1] = Uo[iMax-1][j]
        for i in range(1, iMax-1):
            D[i] = d2*Uo[i][j+1] + (1.0 - 2.0*d2)*Uo[i][j] +\
                   d2*Uo[i][j-1]
        UU = TridiagonalSolver(iMax, A, B, C, D, UU)
        for i in range(1, iMax-1):
            # Alternating Direction Implicit method in x-direction
            Uhalf[i][j] = UU[i]

    # *****************************************************************
    # This block of codes solves for U at time level n + 1
    # along constant i line
    # Eq. 5.25 using Tridiagonal system Appendix B in
    # Hoffmann CFD Vol.1
    # *****************************************************************
    A = [-d2 for j in range(jMax)]
    B = [(1.0 + 2.0*d2) for j in range(jMax)]
    C = [-d2 for j in range(jMax)]
    D = [0 for j in range(jMax)]
    UU = [0 for j in range(jMax)]
    for i in range(1, iMax-1):
        UU[0] = Uo[i][0]
        UU[-1] = Uo[i][jMax-1]
        for j in range(1, jMax-1):
            D[j] = d1*Uhalf[i+1][j] + (1.0 - 2.0*d1)*Uhalf[i][j] +\
                   d1*Uhalf[i-1][j]
        UU = TridiagonalSolver(jMax, A, B, C, D, UU)
        for j in range(1, jMax-1):
            # Alternating Direction Implicit method in y-direction
            U[i][j] = UU[j]

    return U
