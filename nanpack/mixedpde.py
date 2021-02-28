"""A module to obtain the numerical solution of viscous Burgers equation.

Various numerical methods are included in this module that can be used
to obtain the solution of the model equation.

The viscous Burgers equation is a scalar representation of the Navier-
Stokes equations which are classified as the mixed hyperbolic-parabolic
equation.

The non-linear viscous Burgers equation is expressed as:

    du/dt + u(du/dx) = nu(d2u/dx2)
or,
    du/dt + dE/dx = nu(d2u/dx2)
or,
    du/dt + A(du/dx) = nu(d2u/dx2)

    where,
                u: measurable quanity
                nu : diffusion coefficient
                E = u^2/2
                A = dE/du

The equation is non-dimensionalized using the following expression:
    x* = x/L
    u* = uL/nu
    t* = nu(t/L2)
"""
#   ***********************************************************************
#
#   FILE         mixedpdesolvers.py
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
from .tridiagonal import TridiagonalSolver
from .backend.exceptions import DimensionError


def FTCS(Uo, Courant, diffX):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the explicit Forward Time Central Spacing (FTCS)
    method
    to obtain the solution of the
    1-D non-linear viscous Burgers equation.

    Call signature:

        FTCS(Uo, Courant, diffX)

    Parameters
    ----------
    Uo : 1D array

        The dependent variable from time level (n) within the domain.

    Courant : float

        Courant number that appears in the convection component of the PDE.

    diffX : float

        Diffusion number for x-component that appears in the diffusion
        component of the PDE.

    Returns
    -------
    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    shapeU = Uo.shape  # Obtain Dimension
    if len(shapeU) == 2:
        raise DimensionError("2D", "viscous Bergers", "FTCS")

    U = Uo.copy()  # Initialize U
    A = Uo.copy()  # Initialize A
    E = Uo*Uo/2
    iM, = shapeU
    A[1:-1] = (E[2:]-E[1:-1]) / (Uo[2:]-Uo[1:-1])
    U[1:-1] = (
        Uo[1:-1]
        - 0.5*0.5*Courant*(A[2:]+A[0:-2])*(Uo[2:]-Uo[0:-2])
        + diffX*(Uo[2:] - 2.0*Uo[1:-1] + Uo[0:-2])
        )

    return U


def FTBCS(Uo, Courant, diffX):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the explicit FTCS method with backward approximation
    of the convective term
    to obtain the solution of the
    1D non-linear viscous Burgers equation.

    Call signature:

        FTBCS(Uo, Courant, diffX)

    Parameters
    ----------
    Uo : 1D array

        The dependent variable from time level (n) within the domain.

    Courant : float

        Courant number that appears in the convection component of the PDE.

    diffX : float

        Diffusion number for x-component that appears in the diffusion
        component of the PDE.

    Returns
    -------
    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    shapeU = Uo.shape  # Obtain Dimension
    if len(shapeU) == 2:
        raise DimensionError("2D", "viscous Bergers", "FTBCS")

    U = Uo.copy()  # Initialize U
    A = Uo.copy()  # Initialize A
    E = Uo*Uo/2
    A[1:-1] = (E[2:]-E[1:-1]) / (Uo[2:]-Uo[1:-1])
    U[1:-1] = (
        Uo[1:-1]
        - 0.5*Courant*(A[2:]+A[0:-2])*(Uo[1:-1]-Uo[0:-2])
        + diffX*(Uo[2:] - 2.0*Uo[1:-1] + Uo[0:-2])
        )

    return U


def DuFortFrankel(Uo, Uo2, Courant, diffX):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the explicit DuFort-Frankel method
    to obtain the solution of the
    1D non-linear viscous Burgers equation.

    Call signature:

        DuFortFrankel(Uo, Uo2, Courant, diffX)

    Parameters
    ----------
    Uo : 1D array

        The dependent variable from time level (n) within the domain.

    Courant : float

        Courant number that appears in the convection component of the PDE.

    diffX : float

        Diffusion number for x-component that appears in the diffusion
        component of the PDE.

    Returns
    -------
    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    shapeU = Uo.shape  # Obtain Dimension
    if len(shapeU) == 2:
        raise DimensionError("2D", "viscous Bergers", "DuFort-Frankel")

    U = Uo.copy()  # Initialize U
    E = Uo*Uo/2
    A = (E[2:]-E[1:-1])/(Uo[2:]-Uo[1:-1])
    c = A*Courant
    U[1:-1] = (
        ((1.0 - 2.0*diffX) / (1.0 + 2.0*diffX))*Uo2[1:-1]
        + ((c[1:-1] + 2.0*diffX) / (1.0 + 2.0*diffX))*Uo[0:-2]
        - ((c[1:-1] - 2.0*diffX) / (1.0 + 2.0*diffX))*Uo[2:]
        )

    return U


def MacCormack(Uo, Courant, diffX=None):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the explicit MacCormack method
    to obtain the solution of the
    1D non-linear viscous Burgers equation.

    Call signature:

        MacCormack(Uo, Courant, diffX)

    Parameters
    ----------
    Uo : 1D array

        The dependent variable from time level (n) within the domain.

    Courant : float

        Courant number that appears in the convection component of the PDE.

    diffX : float

        Diffusion number for x-component that appears in the diffusion
        component of the PDE.

    Returns
    -------
    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    shapeU = Uo.shape  # Obtain Dimension
    if len(shapeU) == 2:
        raise DimensionError("2D", "viscous Bergers", "expl. MacCormack")

    iMax, = shapeU
    U = Uo.copy()  # Initialize U
    Utemp = Uo.copy()
    Etemp = Uo.copy()

    E = Uo*Uo/2

    for i in range(1, iMax-1):
        dUo = (-Courant*(E[i+1] - E[i])
               + diffX*(Uo[i+1] - 2.0*Uo[i] + Uo[i-1]))
        Utemp[i] = Uo[i] + dUo[i]  # Predictor step
        Etemp[i] = Utemp[i]*Utemp[i]
        dUtemp = (-Courant*(Etemp[i]-Etemp[i-1])
                  + diffX*(Utemp[i+1] - 2.0*Utemp[i] + Utemp[i-1]))
        U[i] = 0.5 * (Uo[i]+Utemp[i]+dUtemp)  # Corrector step

    return U


def BTCS(Uo, Courant, diffX):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the  implicit Euler's Backward Time Central Space
    (BTCS) method
    to obtain the solution of the
    1D non-linear viscous Burgers equation.

    Call signature:

        EulersBTCS(Uo, Courant, diffX)

    Parameters
    ----------
    Uo : 1D array

        The dependent variable from time level (n) within the domain.

    Courant : float

        Courant number that appears in the convection component of the PDE.

    diffX : float

        Diffusion number for x-component that appears in the diffusion
        component of the PDE.

    Returns
    -------
    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    shapeU = Uo.shape  # Obtain Dimension
    if len(shapeU) == 2:
        raise DimensionError("2D", "viscous Bergers", "BTCS")

    iMax, = shapeU
    U = Uo.copy()  # Initialize U

    E = Uo*Uo/2
    A = Uo.copy()
    cc = 0.5*Courant
    A[1:-1] = (E[2:] - E[1:-1])/(Uo[2:] - Uo[1:-1])
    AA = [(-diffX - A[i]*cc) for i in range(iMax)]
    BB = [(1.0 + 2.0*diffX) for i in range(iMax)]
    CC = [(-diffX + A[i]*cc) for i in range(iMax)]
    DD = [Uo[i] for i in range(iMax)]
    UU = Uo.copy()
    U = TridiagonalSolver(iMax, AA, BB, CC, DD, UU)

    return U


def BTBCS(Uo, Courant, diffX, Accuracy):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the  implicit backward time central spacing (BTBCS)
    method with backward differencing approximation for the convective term
    to obtain the solution of the
    1D non-linear viscous Burgers equation.

    Call signature:

        BTBCS(Uo, Courant, diffX, Accuracy)

    Parameters
    ----------
    Uo : 1D array

        The dependent variable from time level (n) within the domain.

    Courant : float

        Courant number that appears in the convection component of the PDE.

    diffX : float

        Diffusion number for x-component that appears in the diffusion
        component of the PDE.

    Returns
    -------
    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    from .backend.initialize import init

    shapeU = Uo.shape  # Obtain Dimension
    if len(shapeU) == 2:
        raise DimensionError("2D", "viscous Bergers", "BTBCS")

    if Accuracy.lower() == "first-order":
        th = {
            "1": 1.0,
            "2": 1.0,
            "3": 0.0,
            "4": 0.0
            }
    elif Accuracy.lower == "second-order":
        th = {
            "1": 2.0,
            "2": 3.0/2,
            "3": 0.0,
            "4": -1.0/2
            }
    elif Accuracy.lower == "third-order":
        th = {
            "1": 1.0,
            "2": 1.0/2,
            "3": 1.0/3,
            "4": -1.0/6
            }
    else:
        raise Exception("Invalid input for argument - Accuracy")

    iMax, = shapeU
    E = Uo*Uo/2
    A, AA, BB, CC, DD, UU = init(Uo, 5)
    A[1:-1] = (E[2:]-E[0:-2]) / (Uo[2:]-Uo[0:-2])
    AA = [(-diffX - th(1)*A[i]*Courant) for i in range(iMax)]
    BB = [(1.0 + 2.0*diffX + th(2)*A[i]*Courant) for i in range(iMax)]
    CC = [(-diffX + th(3)*A[i]*Courant) for i in range(iMax)]
    DD = [(Uo[i] + th(4)*A[i]*Courant*Uo[i-2]) for i in range(2, iMax)]
    UU = Uo.copy()
    U = TridiagonalSolver(iMax, AA, BB, CC, DD, UU)

    return U


def ModifiedRungeKutta(Uo, Courant, diffX):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the explicit four-stage Modified Runge-Kutta method
    to obtain the solution of the
    1D non-linear viscous Burgers equation.

    Call signature:

        ModifiedRungeKutta(Uo, Courant, diffX)

    Parameters
    ----------
    Uo : 1D array

        The dependent variable from time level (n) within the domain.

    Courant : float

        Courant number that appears in the convection component of the PDE.

    diffX : float

        Diffusion number for x-component that appears in the diffusion
        component of the PDE.

    Returns
    -------
    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    shapeU = Uo.shape  # Obtain Dimension
    if len(shapeU) == 2:
        raise DimensionError("2D", "viscous Bergers", "Modified RK")

    U = Uo.copy()  # Initialize U

    # 1st stage
    E = Uo*Uo/2
    U[1:-1] = Uo[1:-1] - Courant*(E[2:]-E[0:-2])/8.0
    # -- update BC here (required when Neumann BC is used)
    # 2nd stage
    E = U*U/2
    U[1:-1] = Uo[1:-1] - Courant*(E[2:]-E[0:-2])/6.0
    # -- update BC here (required when Neumann BC is used)
    # 3rd stage
    E = U*U/2
    U[1:-1] = Uo[1:-1] - Courant*(E[2:]-E[0:-2])/4.0
    # -- update BC here (required when Neumann BC is used)
    # 4th stage
    E = U*U/2
    U[1:-1] = Uo[1:-1] - Courant*(E[2:]-E[0:-2])/2.0
    # -- update BC here (required when Neumann BC is used)

    # Add the viscous terms in the viscous Burgers equation
    # after the final stage, (pg 291. CFD Vol 1 Hoffmann].
    U[1:-1] = U[1:-1] + diffX*(U[2:]-2.0*U[1:-1]+U[0:-2])

    return U


def SecondOrderTVD(Uo, Courant, diffX, LimiterFunc, Limiter, Eps=0.1):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the explicit second-order TVD method and their
    verious Limiter functions and Limiters
    to obtain the solution of the
    1D non-linear viscous Burgers equation.

    Call signature:

        SecondOrderTVD(Uo, Courant, diffX, LimiterFunc, Limiter, Eps)

    Parameters
    ----------
    Uo : 1D array

        The dependent variable from time level (n) within the domain.

    Courant : float

        Courant number that appears in the convection component of the PDE.

    diffX : float

        Diffusion number for x-component that appears in the diffusion
        component of the PDE.

    LimiterFunc: str

        Flux limiter function.

    Limiter:

        Limiter type.

    Eps: float, Default=0.1

        A positive constant in the entropy correction term, si in Eq. 6-127
        in CFD Vol. 1 by Hoffmann. Its value must be between 0 and 0.125.

    Returns
    -------
    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    from .tvdfunctions import CalculateTVD
    import backend.fetchoptions as fo
    from .backend.exceptions import TVDLimiterFunctionInputError

    shapeU = Uo.shape  # Obtain Dimension
    if len(shapeU) == 2:
        raise DimensionError("2D", "viscous Bergers", "second-order TVD")

    iMax, = shapeU
    U = Uo.copy()  # Initialize U
    E = Uo*Uo/2

    fetch = fo.FetchOptions()
    limfunc_options = fetch.TVDLimiterFunctionOptions()

    if LimiterFunc not in limfunc_options:
        raise TVDLimiterFunctionInputError(LimiterFunc)

    for i in range(2, iMax-2):
        phiPlus, phiMinus = CalculateTVD(i, Uo, E, Eps, Courant,
                                         Limiter, LimiterFunc)

        # Equation 6-124 and 6-125 in Hoffmann Vol. 1
        hPlus = 0.5 * (E[i+1]+E[i]+phiPlus)
        hMinus = 0.5 * (E[i]+E[i-1]+phiMinus)

        # Calculate diffusion terms in the viscous Bergers equation
        # Equation 7-58
        diffusion = diffX*(Uo[i+1] - 2.0*Uo[i] + Uo[i-1])

        # Equation 6-123
        U[i] = Uo[i] - Courant*(hPlus-hMinus) + diffusion

    return U
