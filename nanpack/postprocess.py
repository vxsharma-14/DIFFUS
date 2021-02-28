"""A module for post-processing operations."""
#   ***********************************************************************
#
#   FILE         postprocess.py
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


def DimensionalizeSolution(Ustar, RefLength, Diff):
    """Return the dimensional values of the velocity.

    Using the expression:

        u = (u*)(nu/L)

    Call signature:

        DimensionalizeSolution(Ustar, RefLength, Diff)

    Parameters
    ----------
    Ustar : 1D or 2D array

        Non-dimensional velocities.

    RefLength: float

        Reference or characteristic length.

    Diff: float

        Diffusion coefficient.

    Returns
    -------
    Udim : 1D or 2D araay

        Dimensional velocities.
    """
    from .backend.dimensionalize import NonDimensionalize
    nd = NonDimensionalize()
    U = nd.dimVelU(Ustar, RefLength, Diff)
    return U


def ComputeErrorTerm(uAna, uNum, ReturnType="abs"):
    """Return the error term to calculate accuracy of a numerical scheme.

    Calculated using one of the following expression:

        Absolute Error Term = abs(U_analytical - U_numerical)

                        U_analytical - U_numerical
         % Error Term = -------------------------- * 100
                              U_numerical

    Call signature:

        ComputeErrorTerm(uAna, uNum)

    Parameters
    ----------
    uAna: 1D or 2D array

        The values of the dependent variable obtained from the known exact
        solution within the entire domain.

    uNum: 1D or 2D array

        The values of the dependent variable calculated using an
        approximation method within the entire domain.

    ReturnType: str, Default="abs"

        Allowed inputs = "%" or "abs"
        Specify the return type of error term. It can be calculated as a
        percentage difference
        or as an
        absolute difference
        between the numerical and the analytical solution.

    Returns
    -------
    ErrorTerm: 1D or 2D array

        Error term.
    """
    shapeU = uNum.shape  # Obtain Dimension
    ErrorTerm = uNum.copy()  # Initialize ErrorTerm

    if len(shapeU) == 1:
        ErrorTerm[:] = abs(uAna[:] - uNum[:])
        if ReturnType == "%":
            ErrorTerm[:] = (ErrorTerm[:]/uAna[:])*100.0

    elif len(shapeU) == 2:
        ErrorTerm[:, :] = abs(uAna[:, :] - uNum[:, :])
        if ReturnType == "%":
            ErrorTerm[:, :] = (ErrorTerm[:, :]/uAna[:, :])*100.0

    return ErrorTerm


def FourthOrderDamping(U, DampCoeff):
    """Return fourth-order damping term.

    Calculated using Equation 6-77 in CFD Vol. 1 by Hoffmann.

    Call signature:

        FourthOrderDamping(U, DampCoeff)

    Parameters
    ----------
    U : 1D or 2D array

        The dependent variable from time level (n) within the domain.

    DampCoeff : float

        Damping coefficient. The value must be selected in the range of
        0 to 0.125 for a stable solution.

    Returns
    -------
    D : 1D or 2D array

        The fourth-order damping term within the entire domain.
    """
    D = U.copy()  # Initialize D
    D[2:-2] = (
        -DampCoeff*(U[0:-4] - 4.0*U[1:-3] + 6.0*U[2:-2]
                    - 4.0*U[3:-1] + U[4:])
        )
    return D


def SecondOrderDamping(U, DampCoeff):
    """Return second-order damping term.

    Calculated using Equation 6-80 in CFD Vol. 1 by Hoffmann.

    Call signature:

        SecondOrderDamping(U, DampCoeff)

    Parameters
    ----------
    U : 1D or 2D array

        The dependent variable from time level (n) within the domain.

    DampCoeff: float

        Damping coefficient.

    Returns
    -------
    D : 1D or 2D array

        The second-order damping term within the entire domain.
    """
    D = U.copy()  # Initialize D
    D[1:-1] = -DampCoeff*(U[0:-2] - 2.0*U[1:-1] + U[2:])
    return D


def LInfNormError(U, Uold):
    """Return L-infinity Norm error."""
    shapeU = U.shape  # Obtain shape for Dimension
    if len(shapeU) == 2:  # Dimension = 2D
        Abs_err = abs(Uold[1:, 1:] - U[1:, 1:]).sum(axis=1)
        Error = max(Abs_err)

    elif len(shapeU) == 1:  # Dimension = 1D
        Abs_err = abs(Uold[1:] - U[1:]).sum(axis=1)
        Error = max(Abs_err)

    return Error


def AbsoluteError(U, Uold):
    """Return absolute error.

    Absolute error is calculated to track the changes or deviation in
    the numerical solution. This is done by subtracting the solution
    obtained at the current time step from the solution at previous
    time step within the entire domain.

    The computed error is a scalar floating-point value.

    It is a good practice to observe errors during simulation as it helps
    in identifying whether solution is converging or diverging.
    The computed error is also used in stopping the simulation once the
    errors have converged to the set minimum criteria for convergence.

    This error must not be confused with the ERROR TERM which is used for
    calculating accuracy of the numerical method.
    """
    shapeU = U.shape  # Obtain shape for Dimension
    if len(shapeU) == 2:  # Dimension = 2D
        Error = abs(Uold[1:-1, 1:-1] - U[1:-1, 1:-1]).sum()

    elif len(shapeU) == 1:  # Dimension = 1D
        Error = abs(Uold[1:-1] - U[1:-1]).sum()

    return Error


def MonitorConvergence(CfgClsObj, n, Error):
    """Display error data to monitor solution convergence.

    Parameters
    ----------
    CfgClsObj :

        Class object of RunConfig class which was created at the beginning
        of the simulation. Users may also define a new class containing all
        the user-defined paramaters and create an object for the user-
        defined class and provide the object as a parameter of CfgClsObj in
        WriteConvHistToFile function.

    n: int
        Current time-step or iteration level.

    Error: float
        Convergence log at each iteration level calculated in terms of
        absolute error or L2 inf error.

    Returns
    -------
    none
    """
    nDisp = CfgClsObj.nDisplay
    if n == 0 or n == 1:
        print(f'{"ITER":>7} {"ERROR":>15}')
        print(f'{"----":>7} {"-----":>15}')

    if n % nDisp == 0:
        print(f"{n:>7} {Error:>15.8f}")


def WriteSolutionToFile(U, n, nWrite, nMax, OutFileName, dX, dY=None):
    """Write simulation results of the entire domain to output file.

    Call Signature:

        WriteSolutionToFile(U, n, nWrite, nMax, OutFileName, dX, dY=None)

    Parameters
    ----------
    U : 1D or 2D array, float

        Solution of the dependent variable to be stored.

    n : int

        Iteration level/time step level.

    nWrite : int

        After every nWrite values, the solution is stored in the files.
        This helps saving computational processing requirements at every
        time step or iteration level. This value is provided as the user
        input in the configuration file.

    nMax : int

        Maximum allowed time-step or iterations at which the solution is
        stopped. This value is provided as the user input in the
        configuration file.

    OutFileName : str, Default=None.

        File name to store numerical solutions.
    """
    if n % nWrite == 0 or n == nMax:
        OutFile = open(OutFileName, "w")
        shapeU = U.shape  # Obtain shape for Dimension

        if len(shapeU) == 1:  # Dimension = 1D
            iMax, = shapeU
            print(f"IMAX= {iMax}", file=OutFile)
            print(f"Grid Step= {dX:>6.4f}", file=OutFile)
            print(f'{"X":^14} {"U":^14}', file=OutFile)
            for i in range(0, iMax):
                print(f"{dX*i:>12.8e} {U[i]:>12.8e}", file=OutFile)
        elif len(shapeU) == 2:  # Dimension = 2D
            iMax, jMax = shapeU
            print(f"IMAX= {iMax}, JMAX= {jMax}", file=OutFile)
            print(f"Grid Steps= {dX:>6.4f} {dY:>6.4f}",
                  file=OutFile)
            print(f'{"X":^14} {"Y":^14} {"U":^14}', file=OutFile)
            for i in range(0, iMax):
                for j in range(0, jMax):
                    print(f"{dX*i:>12.8e} {dY*j:>12.8e} {U[i][j]:>12.8e}",
                          file=OutFile)

        OutFile.close()


def WriteConvHistToFile(CfgClsObj, n, Error, HistFName=None):
    """Write convergence history log.

    Call Signature:

        WriteConvHistToFile(CfgClsObj, n, Error, HistFName=None)

    Parameters
    ----------
    CfgClsObj :

        Class object of RunConfig class which was created at the beginning
        of the simulation. Users may also define a new class containing all
        the user-defined paramaters and create an object for the user-
        defined class and provide the object as a parameter of CfgClsObj in
        WriteConvHistToFile function.

    n : int

        Iteration level (time level).

    Error : float

        Error in the solution at each iteration level.

    HistFName : str, Default=None

        File name to store log of convergence history.
    """
    if not HistFName:
        HistFName = CfgClsObj.HistFileName
    # Write convergence data to file
    if n == 0 or n == 1:
        HistFile = open(HistFName, "w")
        print(f'{"ITER":>7} {"ERROR":>15}', file=HistFile)
        print(f'{"----":>7} {"-----":>15}', file=HistFile)

    if n % CfgClsObj.nDisplay == 0 or n == CfgClsObj.nMax:
        HistFile = open(HistFName, "a")
        print(f"{n:>7} {Error:>15.8f}", file=HistFile)

        HistFile.close()

    # Write convergence status to output file at the end of file
    if CfgClsObj.State.upper() == "STEADY":
        if n == CfgClsObj.nMax and Error > CfgClsObj.ConvCrit:
            PrintMsg = f"\nSTATUS: NOT CONVERGED\nMAX. ITERATIONS=\
 {CfgClsObj.nMax}"
            print(PrintMsg)
            print()
            HistFile = open(HistFName, "a")
            print(PrintMsg, file=HistFile)
            print("Writing convergence log file: Completed.")
            print("Files saved:")
            print(f'"{HistFName}".')

            HistFile.close()

        elif Error < CfgClsObj.ConvCrit:
            PrintMsg = f"\nSTATUS: CONVERGED SOLUTION OBTAINED AT\
\nITERATIONS={n},\nMAX. ERROR= {CfgClsObj.ConvCrit}"
            print(PrintMsg)
            print()

            HistFile = open(HistFName, "a")
            print(PrintMsg, file=HistFile)
            print("Writing convergence log file: Completed.")
            print("Files saved:")
            print(f'"{HistFName}".')

            HistFile.close()

    elif CfgClsObj.State.upper() == "TRANSIENT" and n == CfgClsObj.nMax:
        PrintMsg = f"\nSTATUS: SOLUTION OBTAINED AT\nTIME LEVEL=\
 {CfgClsObj.totTime} s.\nTIME STEPS= {CfgClsObj.nMax}"
        print(PrintMsg)
        print()

        HistFile = open(HistFName, "a")
        print(PrintMsg, file=HistFile)
        print("Writing convergence log file: Completed.")
        print("Files saved:")
        print(f'"{HistFName}".')

        HistFile.close()


def WriteSolutionIn1DFormat(CfgClsObj, U, Out1DFName=None):
    """Write 2D output data in 1D format at locations of X or Y.

    Call Signature:

        WriteSolutionIn1DFormat(CfgClsObj, U)

    Parameters
    ----------
    CfgClsObj :

        Class object of RunConfig class which was created at the beginning
        of the simulation. Users may also define a new class containing all
        the user-defined paramaters and create an object for the user-
        defined class and provide the object as a parameter of CfgClsObj in
        WriteConvHistToFile function.

    U : 2D array, float

        Solution of the dependent variable to be stored.

    Out1DFName : str, Default=None.

        File name to save output.
    """
    if not Out1DFName:
        Out1DFName = CfgClsObj.Out1DFName
    OutFile = open(Out1DFName, "w")
    shapeU = U.shape
    iMax, jMax = shapeU
    line = []  # initialize a var to store linewise data
    # Format line # 1 in output file which print location points
    locations = ["{0:^12.2f}".format(li) for li in CfgClsObj.nodes]
    # Format line # 2 in output file which prints dashes
    dashes = ["{0:^12}".format('------') for li in CfgClsObj.nodes]

    if CfgClsObj.PrintNodesDir.upper() == "X":
        iX = [int(i/CfgClsObj.dX) for i in CfgClsObj.nodes]
        msg = "Writing data in 1D format at X locations parallel to\
 Y-axis."
        print(msg)
        print(f'{"Y":<5} {"X=":^2}', *locations, sep=' ', file=OutFile)
        print(f'{"---":<8}', *dashes, sep=' ', file=OutFile)
        for j in range(0, jMax):
            del line[:]  # Delete elements in the list for next Y location
            for i in iX:
                # Read output at different X locations at a constant Y and
                # store in one line
                line.append(U[i][j])
                # Convert "line" list elements to float & format the output
                lines = ["{0:12.8f}".format(float(li)) for li in (line)]
            print(f"{CfgClsObj.dY*j:<8.4f}", *lines, sep=' ', file=OutFile)

    elif CfgClsObj.Direction.upper() == "Y":
        jY = [int(j/CfgClsObj.dY) for j in CfgClsObj.nodes]
        msg = "Writing data in 1D format at Y locations parallel to\
 X-axis."
        print(msg)
        print(f'{"X":<5} {"Y=":^2}', *locations, sep=' ', file=OutFile)
        print(f'{"---":<8}', *dashes, sep=' ', file=OutFile)
        for i in range(0, iMax):
            del line[:]  # Delete elements in the list for next X location
            for j in jY:
                # Read output at different Y locations at a constant X and
                # store in one line
                line.append(U[i][j])
                # Convert "line" list elements to float & format the output
                lines = ["{0:12.8f}".format(float(li)) for li in (line)]
            print(f"{CfgClsObj.dX*i:<8.4f}", *lines, sep=' ', file=OutFile)

    print("Writing output in 1D format: Completed.")
    print("Files saved:")
    print(f'"{Out1DFName}".')

    OutFile.close()


def Plot1DResults(dataFiles, **kwargs):
    """Plot results along 1D axis.

    This plotting function utilizes the features provided by the
    matplotlib plotting library. It accepts keyword arguments as inputs.
    If the keyword arguments is not provided, the default values will be
    assigned.

    Call Signature:

        Plot1DResults(dataFiles, **kwargs)

    Parameters
    ----------
    dataFiles: str list

        Provide path of the saved files as a list. Function will read and
        plot the stored data in these files. Atleast one file input is
        required.

    **kwargs = [
        uAxis,
        Legend,
        Markers,
        useFileCol,
        Title,
        xLabel,
        yLabel
        ]

    uAxis: str, default="X"

        Provide the axis of U on the plot. Allowed values are "X" or "Y"
        string.
        If
            uAxis = "X"
            -----------
        the stored numerical solution of the dependent variable will be
        plotted on the X-axis as a function of values on Y-axis.
        Vice-versa, if
            uAxis = "Y"
            -----------
        the numerical solution will be plotted on the Y-axis as a function
        of values on X-axis.

    Legend: str list, Default=None

        String list to display legends in the plot. Optional.
        Plot legends are a useful tool to make the plots informative.
        If Legend is provided, the length of the Legend list must be
        equal to or greater than the length of dataFiles.

    Markers: str list, Default=None

        List of markers. Optional. The allowed markers are given in
        matplotlib documentation. The Plot1DResults make use of the
        matplotlib marker, markevery, markersize features.
        If Markers is provided, the length of the Markers list must be
        equal to or greater than the length of dataFiles.

    useFileCol: int, Default=1

        Column number of the data stored in the file to plot.
        If not provided, the column number 1 will be plotted versus
        column number 0.

    Title: str, Default=None

        Title text of the plot.

    xLabel: str, Default=None

        X-axis label in the plot.

    yLabel: str, Default=None

        Y-axis label in the plot.
    """
    from os import path
    from .backend._plots import _Plot1D
    from .backend.exceptions import InputFileError

    uAxis = kwargs.get("uAxis", "X")
    legend = kwargs.get("Legend", None)
    markers = kwargs.get("Markers", "default")
    useFileCol = kwargs.get("useFileCol", 1)
    title = kwargs.get("Title", None)
    xlbl = kwargs.get("xLabel", None)
    ylbl = kwargs.get("yLabel", None)

    for file in dataFiles:
        if(path.exists(file)) is False:
            raise InputFileError("FileNotfound", file)
    # Markers
    if markers.lower() == "none":
        markers = None
    elif markers.lower() == "default":
        markers = ["<", "o", "^", ">", "v", "+", "p", "s",
                   "D", "h"]
    else:
        if not markers:
            print("No marker list provided. Using defaults")
            markers = ["<", "o", "^", ">", "v", "+", "p", "s",
                       "D", "h"]
        else:
            markers = markers
    _Plot1D(dataFiles, uAxis, legend, markers, useFileCol, title,
            xlbl, ylbl)


def Plot2DResults(dataFiles, **kwargs):
    """Plot 2D results within the simulation domain.

    This plotting function utilizes the features provided by the
    matplotlib plotting library. It uses keyword arguments as inputs.
    If the keyword arguments is not provided, the default values will be
    assigned.

    Call Signature:

        Plot2DResults(dataFiles, iMax, jMax, **kwargs)

    Parameters
    ----------
    dataFiles: str list

        Provide path of the saved files as a list. Function will read and
        plot the stored data in these files. Atleast one file input is
        required.

    iMax: int

        Number of grid points along X-axis.

    jMax: int

        Number of grid points along Y-axis.

    **kwargs = [
        Title,
        xLabel,
        yLabel,
        cbarLabel,
        nPlots,
        nRow,
        nCol
        cMap,
        Shading,
        PlotType
        ]

    uAxis: str, default="X"

        Provide the axis of U on the plot. Allowed values are "X" or "Y"
        string.
        If
            uAxis = "X"
            -----------
        the stored numerical solution of the dependent variable will be
        plotted on the X-axis as a function of values on Y-axis.
        Vice-versa, if
            uAxis = "Y"
            -----------
        the numerical solution will be plotted on the Y-axis as a function
        of values on X-axis.

    Legend: str list, Default=None

        String list to display legends in the plot. Optional.
        Plot legends are a useful tool to make the plots informative.
        If Legend is provided, the length of the Legend list must be
        equal to or greater than the length of dataFiles.

    Markers: str list, Default=None

        List of markers. Optional. The allowed markers are given in
        matplotlib documentation. The Plot1DResults make use of the
        matplotlib marker, markevery, markersize features.
        If Markers is provided, the length of the Markers list must be
        equal to or greater than the length of dataFiles.

    useFileCol: int, Default=1

        Column number of the data stored in the file to plot.
        If not provided, the column number 1 will be plotted versus
        column number 0.

    Title: str, Default=None

        Title text of the plot.

    xLabel: str, Default=None

        X-axis label in the plot.

    yLabel: str, Default=None

        Y-axis label in the plot.
    """
    from os import path
    import warnings
    from .backend._plots import _Plot2D, _MultiPlot2D
    from .backend.exceptions import InputFileError

    countFiles = len(dataFiles)

    pTitle = kwargs.get("Title", None)
    xLabel = kwargs.get("xLabel", None)
    yLabel = kwargs.get("yLabel", None)
    cbarLbl = kwargs.get("cbarLabel", None)
    nPlots = kwargs.get("nPlots", 1)
    cMap = kwargs.get("cMap", "viridis")
    Shading = kwargs.get("Shading", "gouraud")
    nRow = kwargs.get("nRow", 1)
    nCol = kwargs.get("nCol", 1)
    PlotType = kwargs.get("PlotType", "pcolormesh")
    cLevel = kwargs.get("cLevel", None)
    Alpha = kwargs.get("Alpha", 0.5)

    if (countFiles) > 1 and nPlots == 1:
        warnings.warn("More than one files input found for 1 plot.\
 Plotting data from the first file.")
    else:
        for file in dataFiles:
            if(path.exists(file)) is False:
                raise InputFileError("FileNotfound", dataFiles)

    if nPlots == 1:
        _Plot2D(dataFiles[0], pTitle, xLabel, yLabel, cbarLbl, PlotType,
                cMap, Shading, cLevel, Alpha)
    elif nPlots > 1:
        # Set title arrays for multiple plots
        if isinstance(pTitle, str) is True:
            title = [pTitle for i in range(nPlots)]
        elif len(pTitle) == nPlots:
            title = pTitle
        _MultiPlot2D(nPlots, nRow, nCol, dataFiles, title,
                     xLabel, yLabel, cbarLbl, PlotType, cMap,
                     Shading, cLevel, Alpha)
