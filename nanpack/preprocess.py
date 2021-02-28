"""A module for post-processing operations."""
#   ***********************************************************************
#
#   FILE         preprocess.py
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

from .backend._readconfig import _ReadConfig
from .backend.checkconfig import CheckSections, CheckSetupSection


class RunConfig(_ReadConfig):
    """Class to read inputs from configuration file and set-up variables.

    This class is depenedent on the configparser package to read the
    inputs from the file and therefore, for compability,
    always use the most recent version of the config file provided
    with the package.

    Attributes
    ----------
    InFileName: str, Default= "./input/config.ini".

        The string value representing the file to be read for simulation
        inputs. Allowed file extensions: .ini
    """

    def __init__(self, InFileName="./input/config.ini"):
        """Class constructor for the RunConfig class.

        The contructor method calls the members functions to input
        parameters for the simulation set-up from
        the configuration file.

        The input properties for the simulation are stored in the class
        member variables that can be accessed using the
        RunConfig class constructor object.


        Parameters
        ----------
        InFileName: str
            Path to configuration file.
        """
        self.File = InFileName

        print("*******************************************************")
        print("*******************************************************")
        print("Starting configuration.")
        print()
        print("Searching for simulation configuration file in path:")
        print(f'"{InFileName}"')
        super().__init__(InFileName)

        self._CheckSections()  # Check sections in config file
        print("SUCCESS: Configuration file parsing.")
        print("Accessing numerical setup. Completed")
        print("Accessing domain geometry configuration: Completed")
        print("Accessing meshing configuration: Completed.")
        print("Accessing initial condition settings: Completed.")
        print("Accessing boundary condition settings: Completed")
        print("Accessing constant data: Completed.")
        print("Accessing simulation stop settings: Completed.")
        print("Accessing settings for storing outputs: Completed.")

        # Calculate required quantities.
        self.Mesh()
        self.Initial()
        self.BC()
        self.SimStop()
        self.DisplayConfig()

    def _CheckSections(self):
        """Check sections in the configuration file."""
        # Upon initialization -
        #       1. check - if all sections exist
        CheckSections(self.config, self.File)
        CheckSetupSection(self.config, self.State, self.Model,
                          self.Dimension, self.File)

    def Mesh(self):
        """Access meshing inputs from the DOMAIN and MESH sections.

        Call signature :
            self.Mesh()
        """
        from . import mesh

        if self.GridAutoCalc.upper() == 'YES':
            self.iMax, self.jMax = mesh.CalcGridPoints(self.Dimension,
                                                       self.Length,
                                                       self.dX,
                                                       self.Height,
                                                       self.dY)
        elif self.GridAutoCalc.upper() == 'NO':
            self.dX, self.dY = mesh.CalcGridStepsize(self.Dimension,
                                                     self.Length,
                                                     self.iMax,
                                                     self.Height,
                                                     self.jMax)

    def Initial(self):
        """Access intial condition inputs from the IC section.

        Returns the dependent variables with the initial values.

        Call signature :
            RunConfig.Initial()

        Returns
        -------
        U: 1D or 2D array
            Initial conditions for the dependent variable.
        """
        import nanpack.backend.initialize as init
        # *********** INITIAL CONDITIONS *************
        # Firstly, initialize the dependent variable U which is
        # required before assigning any type of intial conditions,
        # whether cold-start or restart.
        self.U = init.InitialCondition(self.Dimension, self.iMax,
                                       self.jMax)
        return self.U

    def BC(self):
        """Access boundary condition inputs from the BC section.

        Update the boundary conditions of the dependent variable and
        return.

        Call signature :
            RunConfig.BC()

        Returns
        -------
        U: 2D array
            Updated boundary conditions for the dependent variable.
        """
        import nanpack.backend.boundary as bound
        # *********** BOUNDARY CONDITIONS *************
        if self.BCfromFile.upper() == 'YES':
            self.BC = bound.ReadBCfromFile(self.BCFileName)
            # Call function to assign 2D BC
            self.U = bound.BC2D(self.U, self.BC, self.dX, self.dY)
            print("Boundary conditions assignment: Completed.")
        elif self.BCfromFile.upper() == 'NO':
            self.U = self.U

        return self.U

    def SimStop(self):
        """Access simulation stop setting inputs from the STOP section.

        Call signature :
            RunConfig.SimStop()
        """
        import nanpack.mesh as mesh
        # *********** SIM STOP SETTINGS ***********
        # Execute this block for:
        # diffusion eq., first-order wave eq. and Burgers eq.
        if not self.Model.upper() == 'LAPLACE':
            self.dT = mesh.CalcTimeStep(self.CFL, self.diff, self.conv,
                                        self.dX, self.dY,
                                        self.Dimension, self.Model)
            self.nMax = mesh.CalcMaxSteps(self.State, self.nMax, self.dT,
                                          self.totTime)
        # Execute this block for Laplace's eq.
        elif self.Model.upper() == 'LAPLACE':
            self.nMax = self.nMax

    def DisplayConfig(self):
        """Display the configuration to the user for verification.

        Call signature :
            RunConfig.DisplayConfig()
        """
        # ************* PRINT CONFIGURATIONS *************
        print()
        print('**********************************************************')
        print(f'CASE DESCRIPTION                {self.Description}')
        print(f'SOLVER STATE                    {self.State}')
        print(f'MODEL EQUATION                  {self.Model}')
        print(f'DOMAIN DIMENSION                {self.Dimension}')
        print(f'    LENGTH                      {self.Length}')
        if self.Dimension.upper() == '2D':
            print(f'    HEIGHT                      {self.Height}')
        print('GRID STEP SIZE')
        print(f'    dX                          {self.dX:5.3f}')
        if self.Dimension.upper() == '2D':
            print(f'    dY                          {self.dY:5.3f}')
        if not self.Model == 'LAPLACE':
            print(f'TIME STEP                       {self.dT:5.3f}')
        print('GRID POINTS')
        print(f'    along X                     {self.iMax}')
        if self.Dimension.upper() == '2D':
            print(f'    along Y                     {self.jMax}')
        if self.Model.upper() == 'DIFFUSION':
            print(f'DIFFUSION CONST.                {self.diff:6.4e}')
            print(f'DIFFUSION NUMBER                {self.CFL}')
        elif self.Model.upper() == 'FO_WAVE':
            print(f'CONVECTION CONST.               {self.conv}')
            print(f'COURANT NUMBER                  {self.CFL}')
        elif self.Model.upper() == 'BURGERS':
            print(f'COURANT NUMBER                  {self.CFL}')
        if self.State.upper() == 'STEADY':
            print(f'CONVERGENCE CRIT                {self.ConvCrit}')
            print('MAXIMUM ITERATIONS')
            print(f'IF CONVERGENCE NOT OBSERVED     {self.nMax}')
        elif self.State.upper() == 'TRANSIENT':
            print(f'TOTAL SIMULATION TIME           {self.totTime}')
            print(f'NUMBER OF TIME STEPS            {self.nMax}')
        if self.BCfromFile.upper() == 'YES':
            print(f'BC COFIGURATION FILE            "{self.BCFileName}"')
        print(f'START CONDITION                 {self.StartOpt}')
        print('**********************************************************')
        print('SUCEESS: Configuration completed.')
        print()


def BC2D(U, BC, dX, dY):
    """Assign boundary conditions. Add more info."""
    from .backend.boundary import BC2D
    u = BC2D(U, BC, dX, dY)

    return u


def CourantNumber(CFL, Dimension):
    """Return the Courant Number.(this function needs modification)."""
    if Dimension.upper() == '1D':
        convX = CFL  # convection coefficient for the x-term.

    elif Dimension.upper() == '2D':
        convX = CFL  # convection coefficient for the x-term.
        convY = CFL  # convection coefficient for the y-term.
    print('Calculating coefficient in the convective equation: Completed.')

    return convX, convY


def DiffusionNumbers(Dimension, diff, dT, dX, dY=None, Scaling=False):
    """Return the diffusion numbers along X and Y directions.

    Diffusion number is expressed as
        d_x = (nu)dT/dX^2
    In non-dimensional equations, diffusion number is
        d_x = (dT*)/(dX*^2)

    Call signature:
        DiffusionNumbers(Dimension, diff, dT, dX, dY)

    Parameters
    ----------
    Dimension : string
        Dimension of the simulation domain.
        Available Options are "1D" or "2D"
    diff: float
        Constant coefficient in the diffusion equation, such as kinematic
        viscosity.
    dT: float
        Time step size.
    dX: float
        Grid step size in x-direction.
    dY: float
        Grid step size in y-direction. Default=None for 1D applications.
    Scaling: bool, Default=False
        User choice for returning dimensional or non-dimensional diffusion
        numbers.

    Returns
    -------
    diffX, diffY : float values
        Diffusion numbers along X and Y axis, respectively.
    """
    if Scaling is False:
        dX2 = dX*dX
        # Diffusion coefficient for the x term.
        diffX = diff*dT/dX2
        diffY = None
        if Dimension.upper() == "2D":
            dY2 = dY*dY
            diffY = diff*dT/dY2  # diffusion coefficient for y term.
        print("Calculating diffusion numbers: Completed.")

    elif Scaling is True:
        dX2 = dX*dX
        diffX = dT/dX2
        diffY = None
        if Dimension.upper() == "2D":
            dY2 = dY*dY
            diffY = dT/dY2

    return diffX, diffY


def NonDimensionalizeTime(dTime, RefLength, Diff):
    """Return the non-dimensionalized time or time step size.

    Using the expression:
        t* = (nu)(t/L^2)

    Call signature:
        NondimensionalizeTime(dTime, RefLength, Diff)

    Parameters
    ----------
    dTime : float
        Dimensional time or the time step.
    RefLength: float
        Reference or characteristic length.
    Diff: float
        Diffusion coefficient.

    Returns
    -------
    Tstar : float
        Non-dimensional quantity of time.
    """
    from .backend.dimensionalize import NonDimensionalize
    nd = NonDimensionalize()
    Tstar = nd.ndTime(dTime, RefLength, Diff)
    return Tstar


def NonDimensionalizeMesh(dXgrid, RefLength):
    """Return the non-dimensionalized grid locations or grid step size.

    Using the expression:
        x* = x/L

    Call signature:
        NondimensionalizeTime(dTime, refLength, Diff)

    Parameters
    ----------
    dXgrid : float
        Dimensional grid points locations or grid step size.
    RefLength: float
        Reference or characteristic length.

    Returns
    -------
    Xstar : float
        Non-dimensional quantity of grid.
    """
    from .backend.dimensionalize import NonDimensionalize
    nd = NonDimensionalize()
    Xstar = nd.ndGrid(dXgrid, RefLength)
    return Xstar
