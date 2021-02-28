"""A module for computing exact solutions.

This module contains several known analytical solutions for the simple
physical processes that are useful to validate the developed codes and
numerical methods.
"""
#   ***********************************************************************
#
#   FILE         benchmark.py
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


def HeatConduction(X, Y, T1, T2, T3, T4, Inf):
    """Return the analytical solution of the steady-state heat conduction.

    Use the results to validate the numerical solution of
    - the parabolic 2D unsteady heat conduction PDE, or
    - the elliptic 2D steady state heat conduction PDE

    Case description: A rectangular bar is initially heated to a
    temperature T0 (or initial conditions are guessed for steady state
                    simulation).
    Subsequently, its surfaces are subject to the constant temperatures of
    T1, T2, T3 and T4 as depicted below.

                            Y-axis
                            |
                            |          T3
                             _______________________
                            |                       |
                            |                       |
                            |                       |
                            |                       |
                        T2  |                       | T4
                            |                       |
                            |                       |
                            |                       |
                            |                       |
                            |_______________________|   _____ X-axis
                                       T1

    Call signature:
        HeatConduction(X, Y, T1, T2, T3, T4, Inf)

    Parameters
    ----------
    X, Y: 2D array
        Mesh data containing X and Y locations at all grid points.
    T1, T2, T3, T4: float
        Boundary temperatures (in degree Kelvin units) as specified in the
        above schematic.
    Inf: int
        This number corresponds to infinity in the summation series in
        exact solution.
        Start with a smaller number (approx. 20) and user larger value
        if the solution is not converged when Inf=20.

    Returns
    -------
    Tanaly: 2D array
        Analytical solution of the steady state heat conduction within the
        entire domain.
    """
    import numpy as np
    import math

    print("Starting calculations to obtain analytical solution.")

    L = X[-1, 0]
    W = Y[0, -1]

    iMax, jMax = X.shape

    Ta = np.zeros((iMax, jMax))
    Tb = np.zeros((iMax, jMax))
    Tc = np.zeros((iMax, jMax))
    Td = np.zeros((iMax, jMax))

    Ta1 = np.zeros((iMax, jMax))
    Tb1 = np.zeros((iMax, jMax))
    Tc1 = np.zeros((iMax, jMax))
    Td1 = np.zeros((iMax, jMax))

    Ta2 = np.zeros((iMax, jMax))
    Tb2 = np.zeros((iMax, jMax))
    Tc2 = np.zeros((iMax, jMax))
    Td2 = np.zeros((iMax, jMax))

    Ta3 = np.zeros((iMax, jMax))
    Tb3 = np.zeros((iMax, jMax))
    Tc3 = np.zeros((iMax, jMax))
    Td3 = np.zeros((iMax, jMax))

    TA = np.zeros((iMax, jMax))
    TB = np.zeros((iMax, jMax))
    TC = np.zeros((iMax, jMax))
    TD = np.zeros((iMax, jMax))
    TAnaly = np.zeros((iMax, jMax))

    m = 0
    while m <= Inf:

        m = m + 1
        mPi = m*math.pi

        for i in range(1, iMax-1):
            for j in range(1, jMax-1):

                # *********************************************************
                # Compute Ta
                # *********************************************************
                Ta1[i][j] = (1.0 - math.cos(mPi))/mPi
                Ta2[i][j] = (
                    math.sinh(mPi*(W - Y[i][j])/L) / math.sinh(mPi*W/L)
                    )
                Ta3[i][j] = math.sin(mPi*X[i][j]/L)

                Ta[i][j] = Ta[i][j] + Ta1[i][j]*Ta2[i][j]*Ta3[i][j]

                # *********************************************************
                # Compute Tb
                # *********************************************************
                Tb1[i][j] = (1.0 - math.cos(mPi))/mPi
                Tb2[i][j] = math.sinh(mPi*Y[i][j]/L)/math.sinh(mPi*W/L)
                Tb3[i][j] = math.sin(mPi*X[i][j]/L)

                Tb[i][j] = Tb[i][j] + Tb1[i][j]*Tb2[i][j]*Tb3[i][j]

                # *********************************************************
                # Compute Tc
                # *********************************************************
                Tc1[i][j] = (1.0 - math.cos(mPi))/mPi
                Tc2[i][j] = (
                    math.sinh(mPi*(L - X[i][j])/W) / math.sinh(mPi*L/W)
                    )
                Tc3[i][j] = math.sin(mPi*Y[i][j]/W)

                Tc[i][j] = Tc[i][j] + Tc1[i][j]*Tc2[i][j]*Tc3[i][j]

                # *********************************************************
                # Compute Td
                # *********************************************************
                Td1[i][j] = (1.0 - math.cos(mPi))/mPi
                Td2[i][j] = math.sinh(mPi*X[i][j]/W)/math.sinh(mPi*L/W)
                Td3[i][j] = math.sin(mPi*Y[i][j]/W)

                Td[i][j] = Td[i][j] + Td1[i][j]*Td2[i][j]*Td3[i][j]

    TA[1:-1, 1:-1] = T1*2.0*Ta[1:-1, 1:-1]
    TB[1:-1, 1:-1] = T3*2.0*Tb[1:-1, 1:-1]
    TC[1:-1, 1:-1] = T2*2.0*Tc[1:-1, 1:-1]
    TD[1:-1, 1:-1] = T4*2.0*Td[1:-1, 1:-1]

    TAnaly[1:-1, 1:-1] = TA[1:-1, 1:-1] + TB[1:-1, 1:-1]\
        + TC[1:-1, 1:-1] + TD[1:-1, 1:-1]

    # At boundaries
    TAnaly[:, 0] = T1   # At Y=0
    TAnaly[:, -1] = T3  # At Y=H
    TAnaly[0, :] = T2   # At X=0
    TAnaly[-1, :] = T4  # At X=L

    # At corners
    TAnaly[0, 0] = T1  # At X=0, Y=0
    TAnaly[-1, 0] = T1  # At X=L, Y=0
    TAnaly[0, -1] = T2  # At X=0, Y=H
    TAnaly[-1, -1] = T3  # At x=L, Y=H

    print("Calculating analytical solution: Completed.")

    return TAnaly
#   ***********************************************************************


def ParallelPlateFlow(Uref, X, nu, t, Inf):
    """Return analytical solution of the velocity between parallel plates.

    Call signature:
        ParallelPlateFlow(Uref, X, nu, t, Inf)

    Parameters
    ----------
    Uref: float
        Initial velocity of the lower plate.
    X: 1D array
        Mesh data containing X locations at all grid points between the
        plates.
    nu: float
        Kinematic viscosity of the fluid in the governing N-S equation.
    t: float
        Time at which the analytical solution is required.
    Inf: int
        This number corresponds to infinity in the summation series in
        exact solution. Use its value as 20 which is a tested value for
        converged results.

    Returns
    -------
    Ua: 1D array
        Analytical solution of the fluid flow distribution between the
        plates.
    """
    import math
    import numpy as np

    iM, = X.shape

    Ua = np.zeros((iM))
    Ua1 = np.zeros((iM))
    Ua2 = np.zeros((iM))

    eta1 = X[-1]/(2*math.sqrt(nu*t))
    m = 0
    while m <= Inf:
        for i in range(1, iM-1):
            eta = X[i]/(2*math.sqrt(nu*t))
            Ua1[i] = Ua1[i] + math.erfc(2*m*eta1 + eta)
            Ua2[i] = Ua2[i] + math.erfc(2*eta1*(m+1) - eta)
        m = m + 1  # Increment m

    Ua[1:-1] = Uref * (Ua1[1:-1]-Ua2[1:-1])
    Ua[0] = Uref
    Ua[-1] = 0.0

    return Ua


def ViscousBurgersSolution(X, t):
    """Return the analytical solution of a stationary Burgers equation.

    Call signature:
        ViscousBurgersSolution(X, t)

    Parameters
    ----------
    X: 1D array
        Mesh data containing X locations at all grid points within the
        domain.
    t: float
        Time at which the analytical solution is required.

    Returns
    -------
    Uana: 1D array
        Analytical solution of the viscous Burgers equation.
    """
    import math

    Uana = X.copy()
    shapeX = X.shape
    iMax, = shapeX

    for i in range(0, iMax):
        Uana[i] = -2.0 * math.sinh(X[i]) / (math.cosh(X[i])
                                            - math.exp(-t))
    return Uana

#   ***********************************************************************
