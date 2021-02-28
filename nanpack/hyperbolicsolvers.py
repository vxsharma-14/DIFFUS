"""A module to obtain the numerical solution of hyperbolic PDEs.

Various numerical methods are included in this module that can be used
to obtain the solution of the model equation.

The wave equation and the inviscid Burgers equation are the hyperbolic
partial differential equation.

The first-order wave equation is a linear equation which is expressed as:

    du/dt = -a(du/dx)       for a>0

    where,
                u: measurable quanity
                a: constant speed

The first-order inviscid Burgers equation is a non-linear equation
which is expressed as:

    du/dt = -u(du/dx)
or,
    du/dt = -dE/dx

    where,
                E = u^2/2

The non-linear viscous Burgers equation is a scalar representation of
the Navier-Stokes equation. It is expressed as:

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
"""
#   ***********************************************************************
#
#   FILE         hyperbolicsolvers.py
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


def ExplicitFirstUpwind(cfg, Uo, Courant):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the explicit first upwind differencing method
    to obtain the solution of the
    first-order 1D wave equation or inviscid Burgers equation.

    Call signature:

        FirstOrderUpwind(cfg, Uo, Courant)

    Parameters
    ----------
    cfg :

        Class object of RunConfig class which was created at the
        beginning of the simulation.

    Uo : 1D array

        The dependent variable from time level (n) within the domain.

    Courant : float

        Courant number that appears in the convection component of the PDE.

    Returns
    -------
    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    import numpy as np
    shapeU = Uo.shape  # Obtain Dimension
    if len(shapeU) == 2:
        raise DimensionError("2D", "inviscid Bergers", "First Upwind")

    U = Uo.copy()  # Initialize U
    if cfg.Model.upper() == "FO_WAVE":
        positive_a = 1 + np.sign(cfg.conv)  # for positive a
        negative_a = 1 - np.sign(cfg.conv)  # for negative a
        U[1:-1] = (
            Uo[1:-1]
            - 0.5*Courant*positive_a*(Uo[1:-1]-Uo[0:-2])
            - 0.5*Courant*negative_a*(Uo[2:]-Uo[1:-1])
            )

    elif cfg.Model.upper() == "INV_BURGERS":
        E = Uo*Uo/2
        U[1:-1] = Uo[1:-1] - Courant*(E[1:-1]-E[0:-2])

    return U


def Lax(cfg, Uo, Courant):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the explicit Lax method
    to obtain the solution of the
    first-order 1D wave equation or inviscid Burgers equation.

    Call signature:

        Lax(cfg, Uo, Courant)

    Parameters
    ----------
    cfg :

        Class object of RunConfig class which was created at the
        beginning of the simulation.

    Uo : 1D array

        The dependent variable from time level (n) within the domain.

    Courant : float

        Courant number that appears in the convection component of the PDE.

    Returns
    -------
    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    shapeU = Uo.shape  # Obtain Dimension
    if len(shapeU) == 2:
        raise DimensionError("2D", "inviscid Bergers", "Lax")

    U = Uo.copy()  # Initialize U
    if cfg.Model.upper() == "FO_WAVE":
        U[1:-1] = (
            0.5*(Uo[2:]+Uo[0:-2]) - 0.5*Courant*(Uo[2:]-Uo[0:-2])
            )

    elif cfg.Model.upper() == "INV_BURGERS":
        E = Uo*Uo/2
        U[1:-1] = (
            0.5*(Uo[2:]+Uo[0:-2]) - 0.25*Courant*(E[2:]-E[0:-2])
            )

    return U


def MidpointLeapfrog(cfg, Uo, Courant):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the explicit Midpoint Leapfrog method
    to obtain the solution of the
    first-order 1D wave equation or inviscid Burgers equation.

    Call signature:

        MidpointLeapfrog(cfg, Uo, Courant)

    Parameters
    ----------
    cfg :

        Class object of RunConfig class which was created at the
        beginning of the simulation.

    Uo : 1D array

        The dependent variable from time level (n) within the domain.

    Courant : float

        Courant number that appears in the convection component of the PDE.

    Returns
    -------
    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    shapeU = Uo.shape  # Obtain Dimension
    if len(shapeU) == 2:
        raise DimensionError("2D", "inviscid Bergers", "Midpoint Leapfrog")

    if cfg.Model.upper() == "FO_WAVE":
        raise Exception("This formulation is not available for WAVE\
 equation in this version.")

    elif cfg.Model.upper() == "INV_BURGERS":
        raise Exception("This formulation is not available for BURGERS\
 equation in this version.")

    # return U


def LaxWendroff(cfg, Uo, Courant):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the explicit Lax-Wendroff method
    to obtain the solution of the
    first-order 1D wave equation or inviscid Burgers equation.

    Call signature:

        LaxWendroff(cfg, Uo, Courant)

    Parameters
    ----------
    cfg :

        Class object of RunConfig class which was created at the
        beginning of the simulation.

    Uo : 1D array

        The dependent variable from time level (n) within the domain.

    Courant : float

        Courant number that appears in the convection component of the PDE.

    Returns
    -------
    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    shapeU = Uo.shape  # Obtain Dimension
    if len(shapeU) == 2:
        raise DimensionError("2D", "inviscid Bergers", "Lax-Wendroff")

    U = Uo.copy()  # Initialize U
    if cfg.Model.upper() == "FO_WAVE":
        Courant2 = Courant*Courant
        U[1:-1] = (
            Uo[1:-1] - 0.5*Courant*(Uo[2:]-Uo[0:-2])
            + 0.5*Courant2*(Uo[2:] - 2.0*Uo[1:-1] + Uo[0:-2])
            )

    elif cfg.Model.upper() == "INV_BURGERS":
        E = Uo*Uo/2
        Courant2 = Courant*Courant
        U[1:-1] = (
            Uo[1:-1] - 0.5*Courant*(E[2:]-E[0:-2])
            + 0.25*Courant2*((Uo[2:]+Uo[1:-1])*(E[2:]-E[1:-1])
                             - (Uo[1:-1]+Uo[0:-2])*(E[1:-1]-E[0:-2]))
            )

    return U


def LaxWendroffMultiStep(cfg, Uo, Courant):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the explicit multi-step Lax-Wendroff method
    to obtain the solution of the
    first-order 1D wave equation.

    Call signature:

        LaxWendroffMultiStep(cfg, Uo, Courant)

    Parameters
    ----------
    cfg :

        Class object of RunConfig class which was created at the
        beginning of the simulation.

    Uo : 1D array

        The dependent variable from time level (n) within the domain.

    Courant : float

        Courant number that appears in the convection component of the PDE.

    Returns
    -------
    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    shapeU = Uo.shape  # Obtain Dimension
    if len(shapeU) == 2:
        raise DimensionError("2D", "inviscid Bergers", "multi-step L-Wen.")

    if cfg.Model.upper() == "FO_WAVE":
        U = Uo.copy()  # Initialize U
        Uhalf = Uo.copy()
        iMax, = shapeU
        for i in range(1, iMax-1):
            Uhalf[i] = 0.5*(Uo[i+1]+Uo[i]) - 0.5*Courant*(Uo[i+1]-Uo[i])
            U[i] = Uo[i] - Courant*(Uhalf[i]-Uhalf[i-1])

    elif cfg.Model.upper() == "INV_BURGERS":
        raise Exception("This formulation is not available for BURGERS\
 equation in this version.")

    return U


def MacCormack(cfg, Uo, Courant, diffX=None):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the explicit MacCormack method
    to obtain the solution of the
    first-order 1D wave equation or inviscid/viscous Burgers equation.

    Call signature:

        MacCormack(cfg, Uo, Courant, diffX)

    Parameters
    ----------
    cfg :

        Class object of RunConfig class which was created at the
        beginning of the simulation.

    Uo : 1D array

        The dependent variable from time level (n) within the domain.

    Courant : float

        Courant number that appears in the convection component of the PDE.

    Returns
    -------
    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    shapeU = Uo.shape  # Obtain Dimension
    if len(shapeU) == 2:
        raise DimensionError("2D", "inviscid Bergers", "MacCormack")

    iMax, = shapeU
    U = Uo.copy()  # Initialize U
    Utemp = Uo.copy()
    if cfg.Model.upper() == "FO_WAVE":
        for i in range(1, iMax-1):
            Utemp[i] = Uo[i] - Courant*(Uo[i+1]-Uo[i])  # Predictor step
            U[i] = 0.5*((Uo[i]+Utemp[i])
                        - Courant*(Utemp[i]-Utemp[i-1]))  # Corrector step

    elif cfg.Model.upper() == "INV_BURGERS":
        E = Uo*Uo/2
        Etemp = Utemp*Utemp/2
        for i in range(1, iMax-1):
            Utemp[i] = Uo[i] - Courant*(E[i+1]-E[i])  # Predictor step
            Etemp[i] = Utemp[i]*Utemp[i]/2
            U[i] = 0.5*((Uo[i]+Utemp[i])
                        - Courant*(Etemp[i]-Etemp[i-1]))  # Corrector step

    return U


def FourthOrderRungeKutta(cfg, Uo, Courant):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the explicit four-stage Runge-Kutta method
    to obtain the solution of the
    first-order 1D wave equation or inviscid Burgers equation.

    Call signature:

        FourthOrderRungeKutta(cfg, Uo, Courant)

    Parameters
    ----------
    cfg :

        Class object of RunConfig class which was created at the
        beginning of the simulation.

    Uo : 1D array

        The dependent variable from time level (n) within the domain.

    Courant : float

        Courant number that appears in the convection component of the PDE.

    Returns
    -------
    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    shapeU = Uo.shape  # Obtain Dimension
    if len(shapeU) == 2:
        raise DimensionError("2D", "inviscid Bergers", "4th order RK")

    U = Uo.copy()  # Initialize U
    if cfg.Model.upper() == "FO_WAVE":
        U1 = Uo.copy()
        U2 = Uo.copy()
        U3 = Uo.copy()
        # 1st stage
        U1[1:-1] = Uo[1:-1] - 0.5*Courant*(Uo[2:]-Uo[0:-2])/2.0
        # 2nd stage
        U2[1:-1] = Uo[1:-1] - 0.5*Courant*(U1[2:]-U1[0:-2])/2.0
        # 3rd stage
        U3[1:-1] = Uo[1:-1] - Courant*(U2[2:]-U2[0:-2])/2.0
        # 4th stage
        U[1:-1] = (
            Uo[1:-1]
            - 0.5*Courant
            * (((1.0/6)*(Uo[2:] - Uo[0:-2]))
               + ((1.0/3)*(U1[2:] - U1[0:-2]))
               + ((1.0/3)*(U2[2:] - U2[0:-2]))
               + ((1.0/6)*(U3[2:] - U3[0:-2])))
            )

    elif cfg.Model.upper() == "INV_BURGERS":
        # 1st stage
        U1 = Uo.copy()
        U2 = Uo.copy()
        U3 = Uo.copy()
        E1 = Uo*Uo/2
        U1[1:-1] = Uo[1:-1] - 0.5*Courant*(E1[2:]-E1[0:-2])/2.0
        # 2nd stage
        E2 = U1*U1/2
        U2[1:-1] = Uo[1:-1] - 0.5*Courant*(E2[2:]-E2[0:-2])/2.0
        # 3rd stage
        E3 = U2*U2/2
        U3[1:-1] = Uo[1:-1] - Courant*(E3[2:]-E3[0:-2])/2.0
        # 4th stage
        E4 = U3*U3/2
        U[1:-1] = (
            Uo[1:-1]
            - 0.5*Courant
            * (((1.0/6)*(E1[2:] - E1[0:-2]))
               + ((1.0/3)*(E2[2:] - E2[0:-2]))
               + ((1.0/3)*(E3[2:] - E3[0:-2]))
               + ((1.0/6)*(E4[2:] - E4[0:-2])))
            )

    return U


def ModifiedRungeKutta(cfg, Uo, Courant, diffX):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the explicit four-stage Modified Runge-Kutta method
    to obtain the solution of the
    first-order 1D wave equation or inviscid/viscous Burgers equation.

    Call signature:

        ModifiedRungeKutta(cfg, Uo, Courant, diffX)

    Parameters
    ----------
    cfg :

        Class object of RunConfig class which was created at the
        beginning of the simulation.

    Uo : 1D array

        The dependent variable from time level (n) within the domain.

    Courant : float

        Courant number that appears in the convection component of the PDE.

    Returns
    -------
    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    shapeU = Uo.shape  # Obtain Dimension
    if len(shapeU) == 2:
        raise DimensionError("2D", "inviscid Bergers", "Modified RK")

    U = Uo.copy()  # Initialize U
    if cfg.Model.upper() == "FO_WAVE":
        # 1st stage
        U[1:-1] = Uo[1:-1] - Courant*(U[2:]-U[0:-2])/8.0
        # -- update BC here (required when Neumann BC is used)
        # 2nd stage
        U[1:-1] = Uo[1:-1] - Courant*(U[2:]-U[0:-2])/6.0
        # -- update BC here (required when Neumann BC is used)
        # 3rd stage
        U[1:-1] = Uo[1:-1] - Courant*(U[2:]-U[0:-2])/4.0
        # -- update BC here (required when Neumann BC is used)
        # 4th stage
        U[1:-1] = Uo[1:-1] - Courant*(U[2:]-U[0:-2])/2.0
        # -- update BC here (required when Neumann BC is used)

    elif cfg.Model.upper() == "INV_BURGERS":
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

    return U


def EulersBTCS(cfg, Uo, Courant, diffX=None):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the  implicit Euler's Backward Time Central Space
    (BTCS) method
    to obtain the solution of the
    first-order 1D wave equation or non-linear viscous Burgers equation.

    Call signature:

        EulersBTCS(cfg, Uo, Courant)

    Parameters
    ----------
    cfg :

        Class object of RunConfig class which was created at the
        beginning of the simulation.

    Uo : 1D array

        The dependent variable from time level (n) within the domain.

    Courant : float

        Courant number that appears in the convection component of the PDE.

    Returns
    -------
    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    shapeU = Uo.shape  # Obtain Dimension
    if len(shapeU) == 2:
        raise DimensionError("2D", "inviscid Bergers", "Euler BTCS")

    iMax, = shapeU
    U = Uo.copy()  # Initialize U
    if cfg.Model.upper() == "FO_WAVE":
        cc = 0.5*Courant
        A = [cc for i in range(iMax)]
        B = [-1 for i in range(iMax)]
        C = [-cc for i in range(iMax)]
        D = -Uo
        UU = Uo.copy()
        U = TridiagonalSolver(cfg.iMax, A, B, C, D, UU)

    elif cfg.Model.upper() == "INV_BURGERS":
        raise Exception("This formulation is not available for BURGERS\
 equation in this version.")

    return U


def CrankNicolson(cfg, Uo, Courant):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the  implicit Crank-Nicolson method
    to obtain the solution of the
    first-order 1D wave equation.

    Call signature:

        CrankNicolson(cfg, Uo, Courant)

    Parameters
    ----------
    cfg :

        Class object of RunConfig class which was created at the
        beginning of the simulation.

    Uo : 1D array

        The dependent variable from time level (n) within the domain.

    Courant : float

        Courant number that appears in the convection component of the PDE.

    Returns
    -------
    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    shapeU = Uo.shape  # Obtain Dimension
    if len(shapeU) == 2:
        raise DimensionError("2D", "inviscid Bergers", "Crank-Nicolson")

    if cfg.Model.upper() == "FO_WAVE":
        U = Uo.copy()  # Initialize U
        iMax, = shapeU
        cc = 0.25*Courant
        A = [cc for i in range(iMax)]
        B = [-1.0 for i in range(iMax)]
        C = [-cc for i in range(iMax)]
        D = [0 for i in range(iMax)]
        D[1:-1] = -Uo[1:-1] + 0.25*Courant*(Uo[2:]-Uo[0:-2])
        UU = Uo.copy()

        U = TridiagonalSolver(iMax, A, B, C, D, UU)

    elif cfg.Model.upper() == "INV_BURGERS":
        raise Exception("This formulation is not available for BURGERS\
 equation in this version.")

    return U


def BeamAndWarming(cfg, Uo, Courant):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the Beam and Warming implicit method
    to obtain the solution of the
    first-order 1D wave equation or inviscid Burgers equation.

    Call signature:

        BeamAndWarming(cfg, Uo, Courant)

    Parameters
    ----------
    cfg :

        Class object of RunConfig class which was created at the
        beginning of the simulation.

    Uo : 1D array

        The dependent variable from time level (n) within the domain.

    Courant : float

        Courant number that appears in the convection component of the PDE.

    Returns
    -------
    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    shapeU = Uo.shape  # Obtain Dimension
    if len(shapeU) == 2:
        raise DimensionError("2D", "inviscid Bergers", "Beam and Warming")

    if cfg.Model.upper() == "FO_WAVE":
        raise Exception("This formulation is not available for WAVE\
 equation in this version.")

    elif cfg.Model.upper() == "INV_BURGERS":
        U = Uo.copy()  # Initialize U
        iMax, = shapeU
        E = Uo*Uo/2
        cc = 0.25*Courant
        A = [-cc*Uo[i-1] for i in range(iMax)]
        B = [1.0 for i in range(iMax)]
        C = [cc*Uo[i-1] for i in range(iMax)]
        D = [0.0 for i in range(iMax)]
        # A[1:-1] = -cc*Uo[0:-2]
        # C[1:-1] = cc*Uo[2:]
        D[1:-1] = (
            Uo[1:-1]
            - 0.5*Courant*(E[2:]-E[0:-2])
            + 0.25*Courant*(Uo[2:]*Uo[2:] - Uo[0:-2]*Uo[0:-2])
            )
        UU = Uo.copy()

        U = TridiagonalSolver(iMax, A, B, C, D, UU)

    return U


def FirstOrderTVD(cfg, Uo, Courant):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the explicit first-order TVD method
    to obtain the solution of the
    first-order inviscid Burgers equation.

    Call signature:

        FirstOrderTVD(cfg, Uo, Courant)

    Parameters
    ----------
    cfg :

        Class object of RunConfig class which was created at the
        beginning of the simulation.

    Uo : 1D array

        The dependent variable from time level (n) within the domain.

    Courant : float

        Courant number that appears in the convection component of the PDE.

    Returns
    -------
    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    import nanpack.secondaryfunctions as sec

    shapeU = Uo.shape  # Obtain Dimension
    if len(shapeU) == 2:
        raise DimensionError("2D", "inviscid Bergers", "First-order TVD")

    if cfg.Model.upper() == "FO_WAVE":
        raise Exception("This formulation is not available for WAVE\
 equation in this version.")

    elif cfg.Model.upper() == "INV_BURGERS":
        U = Uo.copy()  # Initialize U
        iMax, = shapeU
        # Utemp = Uo.copy()
        E = Uo*Uo/2
        for i in range(1, iMax-1):
            dUiPlus12 = sec.CalcUi(Uo[i+1], Uo[i])
            dUiMinus12 = sec.CalcUi(Uo[i], Uo[i-1])
            # -- Calculate alpha(i+1/2) using equation 6-98
            if dUiPlus12 == 0:
                alphaiPlus12 = Uo[i]
            else:
                alphaiPlus12 = (E[i+1]-E[i]) / dUiPlus12
            # -- Calculate alpha(i-1/2) using equation 6-100
            if dUiMinus12 == 0:
                alphaiMinus12 = Uo[i]
            else:
                alphaiMinus12 = (E[i]-E[i-1]) / dUiMinus12
            # Equation 6-119 and 6-120 in CFD Vol. 1 by Hoffmann
            phiPlus = abs(alphaiPlus12)*dUiPlus12
            phiMinus = abs(alphaiMinus12)*dUiMinus12
            # Equation 6-117 and 6-118 in CFD Vol. 1 by Hoffmann
            Utemp = Uo[i] - 0.5*Courant*(E[i+1]-E[i-1])
            U[i] = Utemp + 0.5*Courant*(phiPlus-phiMinus)

    return U


def SecondOrderTVD(cfg, Uo, Courant, LimiterFunc, Limiter, Eps=0.1):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the explicit second-order TVD method and their
    verious Limiter functions and Limiters
    to obtain the solution of the
    first-order inviscid Burgers equation.

    Call signature:

        SecondOrderTVD(cfg, Uo, Courant)

    Parameters
    ----------
    cfg :

        Class object of RunConfig class which was created at the
        beginning of the simulation.

    Uo : 1D array

        The dependent variable from time level (n) within the domain.

    Courant : float

        Courant number that appears in the convection component of the PDE.

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
        raise DimensionError("2D", "inviscid Bergers", "Second-order TVD")

    if cfg.Model.upper() == "FO_WAVE":
        raise Exception("This formulation is not available for WAVE\
 equation in this version.")

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

        # Equation 6-123
        U[i] = Uo[i] - Courant*(hPlus-hMinus)

    return U


# def FluxCorrectedTransport(cfg, Uo, Courant, Damp1, Damp2):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the explicit Flux Corrected Transport for the Lax-
    Wendroff method
    to obtain the solution of the
    first-order 1D wave equation.

    Call signature:

        FluxCorrectedTransport(cfg, Uo, Courant, Damp1, Damp2)

    Parameters
    ----------
    cfg :

           Class object of RunConfig class which was created at the
           beginning of the simulation.

    Uo : 1D array

         The dependent variable from time level (n) within the domain.

    Courant : float

        Courant number that appears in the convection component of the PDE.

    Damp1 : float

            Damping term which is added to the predictor step.

    Damp2 : float

            Antii-diffusive term which is added to the corrector step
            to remove excessive damping.

    Returns
    -------
    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """


'''    import fluid.secondaryfunctions as sf
    shapeU = Uo.shape # Obtain Dimension
    if len(shapeU) == 2:
        raise Exception("This formulation is only available for 1D first\
 order wave equation or the inviscid Burgers equation in this version.")

    U = Uo.copy() # Initialize U
    Utemp = Uo.copy()
    if init.Model.upper() == 'FO_WAVE':
        Courant2 = Courant*Courant
        for i in range(1, iMax-1):
            # Predictor Step
            Utemp[i] = Uo[i]\
                       - 0.5*Courant*(Uo[i+1] - Uo[i-1])\
                       + (Damp1 + 0.5*Courant2)*\
                           (Uo[i+1] - 2.0*Uo[i] + Uo[i-1])
        for i in range(2,iMax-2):
            # Corrector step
            U[i] = Utemp[i]\
                   - Damp2*(Utemp[i+1] - 2.0*Utemp[i] + Utemp[i-1])

    elif init.Model.upper() == 'BURGERS':
        raise Exception("This formulation is not available for BURGERS\
 equation in this version.")

    return U'''
