"""A module to obtain the numerical solution of elliptic PDEs.

Various numerical methods are included in this module that can be used
to obtain the solution of the model equation.

The 2D Poisson's equation is an elliptic partial differential equation.

A typical example of a Poisson's equation is the steady state heat
conduction equation, which is expressed as:

    (d2u/dx2) + (d2u/dy2) = f(x,y)

    where,
                u: measurable quanity

    When f(x,y) = 0, the Poisson's equation reduces to Laplace's equations.
    The nummerical formulations in this modules solves the Laplace's eqn.
"""
#   ***********************************************************************
#
#   FILE         ellipticsolvers.py
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


def PointGaussSeidel(Uo, Beta):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the Point-Gauss Seidel method
    to obtain the solution of the Laplace's equation.

    Call signature:
        PointGaussSeidel(Uo, Beta)

    Parameters
    ----------
    Uo : 2D array
        The dependent variable obtained from the previous iteration
        level, n.
    Beta : float
        Coefficient in the Laplace's finite difference approximation.
        Beta = dX/dY

    Returns
    -------
    U : 2D array
        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    shapeU = Uo.shape  # Obtain Dimension
    if len(shapeU) == 1:
        raise DimensionError("1D", "Laplace's", "Point Gauss-Seidel")
    # Proceed to numerical solution
    U = Uo.copy()  # Initialize U
    iMax, jMax = shapeU
    B2 = Beta*Beta
    A = 0.5/(1.0 + B2)

    # Python numpy array slicing operation cannot be used in PGS method
    # because PGS utilizes solution at k+1 level at points (i-1,j) and
    # (i,j-1) which can only be taken into account using FOR loops.
    # Numpy slicing operation below will result in Jacobi iteration method.
    # U[1:-1,1:-1] = A*(U[2:,1:-1] + U[0:-2,1:-1] +\
    #                   B2*(U[1:-1,2:] + U[1:-1,0:-2]))
    # Another point to note that in the Laplace's SOLVERS formulation, the
    # dependent variable U is used on RHS of discretized eqn instead of Uo
    # as in other MODELS, which is due to the formulation requirement to
    # use
    # values of dependent variable from advanced time step (k+1) at points
    # (i-1,j) or (i,j-1).
    for i in range(1, iMax-1):
        for j in range(1, jMax-1):
            U[i][j] = A*(U[i+1][j] + U[i-1][j]
                         +
                         B2*(U[i][j+1] + U[i][j-1]))

    return U
# *************************************************************************


def LineGaussSeidel_i(Uo, Beta):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the Line-Gauss Seidel method along constant i
    direction (parallel to y-axis)
    to obtain the solution of the Laplace's equation.

    Call signature:
        LineGaussSeidel_i(Uo, Beta)

    Parameters
    ----------
    Uo : 2D array
        The dependent variable obtained from the previous iteration
        level, n.
    Beta : float
        Coefficient in the Laplace's finite difference approximation.
        Beta = dX/dY

    Returns
    -------
    U : 2D array
        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    shapeU = Uo.shape  # Obtain Dimension
    if len(shapeU) == 1:
        raise DimensionError("1D", "Laplace's", "Line Gauss-Seidel")
    # Proceed to numerical solution
    U = Uo.copy()  # Initialize U
    iMax, jMax = shapeU
    B2 = Beta*Beta
    A = [B2 for j in range(jMax)]
    B = [-2.0*(1.0 + B2) for j in range(jMax)]
    C = [B2 for j in range(jMax)]
    D = [0 for j in range(jMax)]
    UU = [0 for j in range(jMax)]
    # NOTE that in the Laplace's SOLVERS formulation, the dependent
    # variable U
    # is used on RHS of discretized eqn instead of Uo as in other MODELS,
    # which is due to the formulation requirement to use values of
    # dependent
    # variable from advanced time steps (k+1) at points (i-1,j) or (i,j-1).
    for i in range(1, iMax-1):
        UU[0] = U[i][0]  # Convert U to 1-D array for Tridiagonal solver
        UU[-1] = U[i][jMax-1]
        for j in range(1, jMax-1):
            D[j] = -(U[i+1][j] + U[i-1][j])
        UU = TridiagonalSolver(jMax, A, B, C, D, UU)
        for j in range(1, jMax-1):
            U[i][j] = UU[j]

    return U
#   ***********************************************************************


def LineGaussSeidel_j(Uo, Beta):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the Line-Gauss Seidel method along constant j
    direction (parallel to x-axis)
    to obtain the solution of the Laplace's equation.

    Call signature:
        LineGaussSeidel_j(Uo, Beta)

    Parameters
    ----------
    Uo : 2D array
        The dependent variable obtained from the previous iteration
        level, n.
    Beta : float
        Coefficient in the Laplace's finite difference approximation.
        Beta = dX/dY

    Returns
    -------
    U : 2D array
        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    shapeU = Uo.shape  # Obtain Dimension
    if len(shapeU) == 1:
        raise DimensionError("1D", "Laplace's", "Line Gauss-Seidel")
    # Proceed to numerical solution
    U = Uo.copy()  # Initialize U
    iMax, jMax = shapeU
    B2 = Beta*Beta
    A = [B2 for i in range(iMax)]
    B = [-2.0*(1.0 + B2) for i in range(iMax)]
    C = [B2 for i in range(iMax)]
    D = [0 for i in range(iMax)]
    UU = [0 for i in range(iMax)]
    # NOTE that in the Laplace's SOLVERS formulation, the dependent
    # variable U
    # is used on RHS of discretized eqn instead of Uo as in other MODELS,
    # which is due to the formulation requirement to use values of
    # dependent
    # variable from advanced time steps (k+1) at points (i-1,j) or (i,j-1).
    for j in range(1, jMax-1):
        UU[0] = U[0][j]  # Convert Uo to 1-D array for Tridiagonal solver
        UU[-1] = U[iMax-1][j]
        for i in range(1, iMax-1):
            D[i] = -(U[i][j+1] + U[i][j-1])
        UU = TridiagonalSolver(iMax, A, B, C, D, UU)
        for i in range(1, iMax-1):
            U[i][j] = UU[i]

    return U
#   ***********************************************************************


def PSOR(Uo, Beta, RelaxParam=1.78):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the Point Successive Over-Relaxation (PSOR) method
    to obtain the solution of the Laplace's equation.

    Call signature:
        PSOR(Uo, Beta, RelaxParam)

    Parameters
    ----------
    Uo : 2D array
        The dependent variable obtained from the previous iteration
        level, n.
    Beta : float
        Coefficient in the Laplace's finite difference approximation.
        Beta = dX/dY
    RelaxParam : float, Default = 1.78
        Relaxation Parameter is used for faster convergence of PSOR method.
        Specify RelaxParam values between 0 and 2.0 to obtain convergence.
        If 0 < RelaxParam < 1: it is called UNDER-RELAXATION.
        If RelaxParam = 1: Point Gauss-Seidel method is recovered.
        An optimum value is determined by performing numerical
        experimentations.
        RelaxParam = 1.78 was found to be an optimum value for PSOR method
        for the problem with a rectangular domain having uniform grid step
        with the Dirichlet BC imposed (see Hoffmann Vol. 1, pg 164, 170).
    Returns
    -------
    U : 2D array
        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    shapeU = Uo.shape  # Obtain Dimension
    if len(shapeU) == 1:
        raise DimensionError("1D", "Laplace's", "Point S Over-Relaxation.")
    # Proceed to numerical solution
    U = Uo.copy()  # Initialize U
    iMax, jMax = shapeU
    B2 = Beta*Beta
    A = 0.5/(1.0 + B2)
    # Python numpy array slicing operation cannot be used in PSOR method
    # because PGS utilizes solution at k+1 level at points (i-1,j) and
    # (i,j-1) which can only be taken into account using FOR loops.
    # Numpy slicing operation is given below but it is not PSOR method.
    # U[1:-1,1:-1] = (1.0 - RelaxParam)*Uo[1:-1,1:-1] +\
    #                 RelaxParam*A*(Uo[2:,1:-1] + Uo[0:-2,1:-1] +\
    #                               B2*(Uo[1:-1,2:] + Uo[1:-1,0:-2]))
    # Another point to note that in the Laplace's SOLVERS formulation, the
    # dependent variable U is used on RHS of discretized eqn instead of Uo
    # as in other MODELS, which is due to the formulation requirement to
    # use the
    # values of dependent variable from advanced time step (k+1) at points
    # (i-1,j) or (i,j-1).
    for i in range(1, iMax-1):
        for j in range(1, jMax-1):
            U[i][j] = (
                (1 - RelaxParam)*U[i][j]
                + RelaxParam*A*(U[i+1][j] + U[i-1][j]
                                + B2*(U[i][j+1] + U[i][j-1]))
                )

    return U
#   ***********************************************************************


def LSOR_i(Uo, Beta, RelaxParam=1.265):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the Line Successive-Over Relaxation (LSOR) method
    along constant i direction (parallel to y-axis).
    to obtain the solution of the Laplace's equation.

    Call signature:
        LSOR_i(Uo, Beta, RelaxParam)

    Parameters
    ----------
    Uo : 2D array
        The dependent variable obtained from the previous iteration
        level, n.
    Beta : float
        Coefficient in the Laplace's finite difference approximation.
        Beta = dX/dY
    RelaxParam : float, Default = 1.265
        Relaxation Parameter is used for faster convergence of LSOR method.
        Specify RelaxParam values between 0 and 2.0 to obtain convergence.
        If 0 < RelaxParam < 1: it is called UNDER-RELAXATION.
        If RelaxParam = 1: Line Gauss-Seidel method is recovered.
        An optimum value is determined by performing numerical
        experimentations.
        RelaxParam = 1.265 was found to be an optimum value for LSOR_i
        method for the problem with a rectangular domain having uniform
        grid step with the Dirichlet BC imposed
        (see Hoffmann Vol. 1, pg 165, 171).

    Returns
    -------
    U : 2D array
        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    shapeU = Uo.shape  # Obtain Dimension

    if len(shapeU) == 1:
        raise DimensionError("1D", "Laplace's", "Line S Over-Relaxation")
    # Proceed to numerical solution
    U = Uo.copy()  # Initialize U
    iMax, jMax = shapeU
    B2 = Beta*Beta
    A = [RelaxParam*B2 for j in range(jMax)]
    B = [-2.0*(1.0 + B2) for j in range(jMax)]
    C = [RelaxParam*B2 for j in range(jMax)]
    D = [0 for j in range(jMax)]
    UU = [0 for j in range(jMax)]
    # NOTE that in the Laplace's SOLVERS formulation, the dependent
    # variable U
    # is used on RHS of discretized eqn instead of Uo as in other MODELS,
    # which is due to the formulation requirement to use values of
    # dependent
    # variable from advanced time steps (k+1) at points (i-1,j) or (i,j-1).
    for i in range(1, iMax-1):
        UU[0] = U[i][0]  # Convert Uo to 1-D array for Tridiagonal solver
        UU[-1] = U[i][jMax-1]
        for j in range(1, jMax-1):
            D[j] = (
                (1.0 - RelaxParam)*B[j]*U[i][j]
                -
                RelaxParam*(U[i+1][j] + U[i-1][j])
                )
        UU = TridiagonalSolver(jMax, A, B, C, D, UU)
        for j in range(1, jMax-1):
            U[i][j] = UU[j]

    return U
#   ***********************************************************************


def LSOR_j(Uo, Beta, RelaxParam=1.265):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the Line Successive-Over Relaxation (LSOR) method
    along constant j direction (parallel to x-axis).
    to obtain the solution of the Laplace's equation.

    Call signature:
        LSOR_j(Uo, Beta, RelaxParam)

    Parameters
    ----------
    Uo : 2D array
        The dependent variable obtained from the previous iteration
        level, n.
    Beta : float
        Coefficient in the Laplace's finite difference approximation.
        Beta = dX/dY
    RelaxParam : float, Default = 1.265
        Relaxation Parameter is used for faster convergence of LSOR method.
        Specify RelaxParam values between 0 and 2.0 to obtain convergence.
        If 0 < RelaxParam < 1: it is called UNDER-RELAXATION.
        If RelaxParam = 1: Line Gauss-Seidel method is recovered.
        An optimum value is determined by performing numerical
        experimentations.
        RelaxParam = 1.265 was found to be an optimum value for LSOR_i
        method for the problem with a rectangular domain having uniform
        grid step with the Dirichlet BC imposed
        (see Hoffmann Vol. 1, pg 165, 171).

    Returns
    -------
    U : 2D array
        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    shapeU = Uo.shape  # Obtain Dimension
    if len(shapeU) == 1:
        raise DimensionError("1D", "Laplace's", "Line S Over-Relaxation")
    # Proceed to numerical solution
    U = Uo.copy()  # Initialize U
    iMax, jMax = shapeU
    B2 = Beta*Beta
    A = [RelaxParam for i in range(iMax)]
    B = [-2.0*(1.0 + B2) for i in range(iMax)]
    C = [RelaxParam for i in range(iMax)]
    D = [0 for i in range(iMax)]
    UU = [0 for i in range(iMax)]
    # NOTE that in the Laplace's SOLVERS formulation, the dependent
    # variable U
    # is used on RHS of discretized eqn instead of Uo as in other MODELS,
    # which is due to the formulation requirement to use values of
    # dependent
    # variable from advanced time steps (k+1) at points (i-1,j) or (i,j-1).
    for j in range(1, jMax-1):
        UU[0] = U[0][j]  # Convert Uo to 1-D array for Tridiagonal solver
        UU[-1] = U[iMax-1][j]
        for i in range(1, iMax-1):
            D[i] = (
                (1.0 - RelaxParam)*B[i]*U[i][j]
                -
                RelaxParam*B2*(U[i][j+1] + U[i][j-1])
                )
        UU = TridiagonalSolver(iMax, A, B, C, D, UU)
        for i in range(1, iMax-1):
            U[i][j] = UU[i]

    return U
#   ***********************************************************************


def ADI(Uo, Beta):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the Alternating Direction Implicit (ADI) method
    to obtain the solution of the Laplace's equation.

    Call signature:
        ADI(Uo, Beta)

    Parameters
    ----------
    Uo : 2D array
        The dependent variable obtained from the previous iteration
        level, n.
    Beta : float
        Coefficient in the Laplace's finite difference approximation.
        Beta = dX/dY

    Returns
    -------
    U : 2D array
        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    shapeU = Uo.shape  # Obtain Dimension
    if len(shapeU) == 1:
        raise DimensionError("1D", "Laplace's", "Alternating Dir. Implict")
    # Proceed to numerical solution
    Uhalf = Uo.copy()  # Uhalf is U at time level (n + 1/2)
    U = Uo.copy()  # Initialize U
    iMax, jMax = shapeU
    B2 = Beta*Beta
    # NOTE that in the Laplace's SOLVERS formulation, the dependent
    # variable U
    # is used on RHS of discretized eqn instead of Uo as in other MODELS,
    # which is due to the formulation requirement to use values of
    # dependent
    # variable from advanced time steps (k+1) at points (i-1,j) or (i,j-1).
    # *********************************************************************
    # This block of codes solves for U at
    # time level n + 1/2 (i.e. = Uhalf) along constant j line
    # Eq. 5.22 using Tridiagonal system Appendix B in
    # Hoffmann CFD Vol.1
    # *********************************************************************
    A = [1.0 for i in range(iMax)]
    B = [-2.0*(1.0 + B2) for i in range(iMax)]
    C = [1.0 for i in range(iMax)]
    D = [0 for i in range(iMax)]
    UU = [0 for i in range(iMax)]
    for j in range(1, jMax-1):
        UU[0] = U[0][j]  # Convert Uo to 1-D array for Tridiagonal solver
        UU[-1] = U[iMax-1][j]
        for i in range(1, iMax-1):
            D[i] = -B2*(U[i][j+1] + Uhalf[i][j-1])
        UU = TridiagonalSolver(iMax, A, B, C, D, UU)
        for i in range(1, iMax-1):
            # Alternating Direction Implicit method in x-direction
            Uhalf[i][j] = UU[i]

    # *********************************************************************
    # This block of codes solves for U at time level n + 1
    # along constant i line
    # Eq. 5.23 using Tridiagonal system Appendix B in
    # Hoffmann CFD Vol.1
    # *********************************************************************
    A = [B2 for j in range(jMax)]
    B = [-2.0*(1 + B2) for j in range(jMax)]
    C = [B2 for j in range(jMax)]
    D = [0 for j in range(jMax)]
    UU = [0 for j in range(jMax)]
    for i in range(1, iMax-1):
        UU[0] = U[i][0]  # Convert Uo to 1-D array for Tridiagonal solver
        UU[-1] = U[i][jMax-1]
        for j in range(1, jMax-1):
            D[j] = -(Uhalf[i+1][j] + U[i-1][j])
        UU = TridiagonalSolver(jMax, A, B, C, D, UU)
        for j in range(1, jMax-1):
            # Alternating Direction Implicit method in y-direction
            U[i][j] = UU[j]

    return U
#   ***********************************************************************


def ADISOR(Uo, Beta, RelaxParam=1.27):
    """Return the numerical solution of dependent variable in the model eq.

    This routine uses the Alternating Direction Implicit Successive Over-
    Relaxation method to obtain the solution of the Laplace's equation.

    Call signature:
        ADISOR(Uo, Beta, RelaxParam)

    Parameters
    ----------
    Uo : 2D array
        The dependent variable obtained from the previous iteration
        level, n.
    Beta : float
        Coefficient in the Laplace's finite difference approximation.
        Beta = dX/dY
    RelaxParam : float, Default = 1.27
        Relaxation Parameter is used for faster convergence of ADISOR
        method.
        Specify RelaxParam values between 0 and 2.0 to obtain convergence.
        If 0 < RelaxParam < 1: it is called UNDER-RELAXATION.
        If RelaxParam = 1: Alternating Direction Implicit method is
        recovered.
        An optimum value is determined by performing numerical
        experimentations.
        RelaxParam = 1.27 was found to be an optimum value for ADISOR
        method for the problem with a rectangular domain having uniform
        grid step with the Dirichlet BC imposed
        (see Hoffmann Vol. 1, pg 172, 183).
    Returns
    -------
    U : 2D array
        The dependent variable calculated at time level (n+1) within the
        entire domain.
    """
    shapeU = Uo.shape  # Obtain Dimension
    if len(shapeU) == 1:
        raise DimensionError("1D", "Laplace's", "ADI S. Over-Relaxation")
    # Proceed to numerical solution
    Uhalf = Uo.copy()  # Uhalf is U at time level (n + 1/2)
    U = Uo.copy()  # Initialize U
    iMax, jMax = shapeU
    B2 = Beta*Beta

    # *********************************************************************
    # This block of codes solves for U
    # at time level n + 1/2 (i.e. = Uhalf) along constant j line
    # Eq. 5.24 using Tridiagonal system Appendix B in
    # Hoffmann CFD Vol.1
    # *********************************************************************
    A = [RelaxParam for i in range(iMax)]
    B = [-2.0*(1.0 + B2) for i in range(iMax)]
    C = [RelaxParam for i in range(iMax)]
    D = [0 for i in range(iMax)]
    UU = [0 for i in range(iMax)]
    for j in range(1, jMax-1):
        UU[0] = U[0][j]  # Convert Uo to 1-D array for Tridiagonal solver
        UU[-1] = U[iMax-1][j]
        for i in range(1, iMax-1):
            D[i] = (
                (1.0 - RelaxParam)*B[i]*U[i][j]
                -
                RelaxParam*B2*(U[i][j+1] + Uhalf[i][j-1])
                )
        # Alternating Direction Implicit SOR method in x-direction
        UU = TridiagonalSolver(iMax, A, B, C, D, UU)
        for i in range(1, iMax-1):
            Uhalf[i][j] = UU[i]

    # *********************************************************************
    # This block of codes solves for U at time level n + 1
    # along constant i line
    # Eq. 5.25 using Tridiagonal system Appendix B in
    # Hoffmann CFD Vol.1
    # *********************************************************************
    A = [RelaxParam*B2 for j in range(jMax)]
    B = [-2.0*(1.0 + B2) for j in range(jMax)]
    C = [RelaxParam*B2 for j in range(jMax)]
    D = [0 for j in range(jMax)]
    UU = [0 for j in range(jMax)]
    for i in range(1, iMax-1):
        UU[0] = U[i][0]  # Convert Uo to 1-D array for Tridiagonal solver
        UU[-1] = U[i][jMax-1]
        for j in range(1, jMax-1):
            D[j] = (
                (1.0 - RelaxParam)*B[j]*Uhalf[i][j]
                -
                RelaxParam*(Uhalf[i+1][j] + U[i-1][j])
                )
        # Alternating Direction Implicit SOR method in y-direction
        UU = TridiagonalSolver(jMax, A, B, C, D, UU)
        for j in range(1, jMax-1):
            U[i][j] = UU[j]

    return U
