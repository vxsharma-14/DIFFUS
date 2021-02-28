#   ***********************************************************************
#
#   FILE         _plots.py
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
import matplotlib.pyplot as plt
import numpy as np


def _Plot1D(dataFiles, uAxis, legend, markers, useFileCol,
            title, xlbl, ylbl):
    """Plot results along 1D axis."""
    print("Preparing data to plot results...")
    countFiles = len(dataFiles)

    # Assign font family and create axis for plotting
    plt.rc('font', family='sans-serif', size=10)
    fig, ax = plt.subplots(dpi=150)

    data = {}  # Initialize a dictionary

    # Use dictionary to create multiple line plots from saved files
    for i in range(countFiles):
        data[f"x{i}"] = np.loadtxt(dataFiles[i], unpack=True, skiprows=3,
                                   usecols=0)
        data[f"u{i}"] = np.loadtxt(dataFiles[i], unpack=True, skiprows=3,
                                   usecols=useFileCol)

    # Start plotting
    print("Plotting 1D results")
    for i in range(countFiles):
        if markers is not None:
            mrkr = markers[i]
        else:
            mrkr = None
        if uAxis == "Y":
            ax.plot(data[f"x{i}"], data[f"u{i}"], linewidth=0.5,
                    label=legend[i], marker=mrkr, markersize=3,
                    markevery=3)
        elif uAxis == "X":
            ax.plot(data[f"u{i}"], data[f"x{i}"], linewidth=0.5,
                    label=legend[i], marker=mrkr, markersize=3,
                    markevery=3)
    # Format and customize plot
    plt.grid(which='both', axis='both', color='lightgrey', linestyle=':',
             linewidth=0.5)
    plt.xlabel(xlbl)
    plt.ylabel(ylbl)
    plt.title(title, fontsize=10)
    ax.legend(fontsize=9)
    # plt.tight_layout()
    plt.show()


def _Plot2D(dataFile, title, xlbl, ylbl, cbarlbl, ptype, cMap, shade,
            clevel, alpha):
    """Plot results within 2D domain."""
    print("Preparing data to plot results...")

    with open(dataFile) as f:
        grid = f.readline().split(",")
    iM = int((grid[0].split("= "))[1])
    jM = int((grid[1].split("= "))[1])

    # Assign font family and create axis for plotting
    plt.rc('font', family='sans-serif', size=10)
    fig, ax = plt.subplots(dpi=150)

    # Read data from saved file
    x, y, u = np.loadtxt(dataFile, unpack=True, skiprows=3)

    # Reshape data based on iM and jM
    X = np.reshape(x, (iM, jM))
    Y = np.reshape(y, (iM, jM))
    U = np.reshape(u, (iM, jM))

    if ptype == "pcolormesh":
        # pcolormesh plot
        plt.pcolormesh(X, Y, U, cmap=cMap, shading=shade)
    elif ptype == "contourf":
        # filled contour plots
        plt.contourf(X, Y, U, cmap=cMap, levels=clevel)
    elif ptype == "contour":
        # contour plots
        plt.contour(X, Y, U, cmap=cMap, levels=clevel)
    elif ptype == "imshow":
        # image plots
        plt.imshow(U, origin="lower", cmap=cMap, interpolation="bilinear")
    elif ptype == "pcolor":
        #
        plt.pcolor(X, Y, U, alpha=alpha, cmap=cMap)

    # Format and customize plot
    plt.xlabel(xlbl, size=10)
    plt.ylabel(ylbl)
    plt.title(title, fontsize=10)
    if cbarlbl is not None:
        cbar = plt.colorbar(format='%.0e')
        cbar.set_label(cbarlbl, rotation=270, size=10, labelpad=15)

    plt.tight_layout()
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()


def _MultiPlot2D(nplots, nrow, ncol, dataFiles, title,
                 xlbl, ylbl, cbarlbl, ptype, cMap, shade, clevel, alpha):
    """Plot multiple results in a single image."""
    print("Preparing data to plot results...")

    # Assign font family and create axis for plotting
    plt.rc('font', family='sans-serif', size=10)
    fig, ax = plt.subplots(dpi=150)

    data = {}  # Initialize a dictionary
    DATA = {}
    iM = []
    jM = []
    for i in range(nplots):
        with open(dataFiles[i]) as f:
            grid = f.readline().split(",")
        iM.append(int((grid[0].split("= "))[1]))
        jM.append(int((grid[1].split("= "))[1]))

    # Use dictionary to create multiple plots
    for i in range(nplots):
        data[f"x{i}"], data[f"y{i}"], data[f"u{i}"] = (
            np.loadtxt(dataFiles[i], unpack=True, skiprows=3))

        # Reshape data based on iM and jM
        DATA[f"X{i}"] = np.reshape(data[f"x{i}"], (iM[i], jM[i]))
        DATA[f"Y{i}"] = np.reshape(data[f"y{i}"], (iM[i], jM[i]))
        DATA[f"U{i}"] = np.reshape(data[f"u{i}"], (iM[i], jM[i]))

    for i in range(nplots):
        plt.subplot(nrow, ncol, i+1)

        if ptype == "pcolormesh":
            # pcolormesh plot
            img = plt.pcolormesh(DATA[f"X{i}"], DATA[f"Y{i}"],
                                 DATA[f"U{i}"],
                                 cmap=cMap, shading=shade)
        elif ptype == "contourf":
            # filled contour plots
            img = plt.contourf(DATA[f"X{i}"], DATA[f"Y{i}"], DATA[f"U{i}"],
                               cmap=cMap, levels=clevel)
        elif ptype == "contour":
            # contour plots
            img = plt.contour(DATA[f"X{i}"], DATA[f"Y{i}"], DATA[f"U{i}"],
                              cmap=cMap, levels=clevel)
        elif ptype == "imshow":
            # image plots
            img = plt.imshow(DATA[f"U{i}"], origin="lower", cmap=cMap,
                             interpolation="bilinear")
        elif ptype == "pcolor":
            #
            img = plt.pcolor(DATA[f"X{i}"], DATA[f"Y{i}"], DATA[f"U{i}"],
                             alpha=alpha, cmap=cMap)

        _cbar_axismap(img, cbarlbl)

        # Format and customize plot

        plt.xlabel(str(xlbl), size=10)
        plt.ylabel(str(ylbl), size=10)
        plt.title(str(title[i]), size=10)
        plt.tight_layout()
        plt.gca().set_aspect("equal", adjustable="box")
    plt.show()


def _cbar_axismap(mappable, ctitle):
    """Control the size of colorbar (aspect ratio) in multi-plot axis."""
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    last_axes = plt.gca()
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar = fig.colorbar(mappable, cax=cax)
    if ctitle is not None:
        cbar.set_label(ctitle, rotation=270, size=12, labelpad=15)
    plt.sca(last_axes)

    return cbar
