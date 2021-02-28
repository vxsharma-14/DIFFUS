"""Not a public module."""
#   ***********************************************************************
#
#   FILE         gridmetrics.py
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


def _metrics_1d(X, dXi=1.0):
    """Return the metrics and Jacobian of the transformation.

    XiX, JJ.
    This function is not complete or tested for accuracy. Doc incomplete.
    """
    iM, = X.shape
    # Initialize
    XXi = X.copy()
    XiX = X.copy()
    JJ = X.copy()
    # At interior grid points
    XXi[1:-1] = (X[2:]-X[0:-2]) / 2.0 / dXi
    # At boundary X = 0.0
    XXi[0] = (-3.0*X[0] + 4.0*X[1] - X[2]) / 2.0 / dXi
    # At boundary X = L
    XXi[-1] = (3.0*X[-1] - 4.0*X[-2] + X[-3]) / 2.0 / dXi
    # Evaluate metrics and Jacobian
    for i in range(0, iM):
        JJ[i] = 1.0 / XXi[i]
        XiX[i] = JJ[i]
    return XiX, JJ


def _metrics_2d(X, Y, dXi=1.0, dEta=1.0):
    """Return the metrics and Jacobian of the transformation.

    XiX, XiY, EtaX, EtaY, JJ. Documentation incomplete.
    """
    iM, jM = X.shape
    # Initialize
    XXi = X.copy()
    YXi = X.copy()
    XEta = X.copy()
    YEta = X.copy()
    XiX = X.copy()
    XiY = X.copy()
    EtaX = X.copy()
    EtaY = X.copy()
    JJ = X.copy()
    # At interior grid points
    XXi[1:-1, 1:-1] = (X[2:, 1:-1]-X[0:-2, 1:-1]) / 2.0 / dXi
    YXi[1:-1, 1:-1] = (Y[2:, 1:-1]-Y[0:-2, 1:-1]) / 2.0 / dXi
    XEta[1:-1, 1:-1] = (X[1:-1, 2:]-X[1:-1, 0:-2]) / 2.0 / dEta
    YEta[1:-1, 1:-1] = (Y[1:-1, 2:]-Y[1:-1, 0:-2]) / 2.0 / dEta
    # At boundary X = 0.0
    XXi[0, 1:-1] = (
        (-3.0*X[0, 1:-1] + 4.0*X[1, 1:-1] - X[2, 1:-1]) / 2.0 / dXi
        )
    YXi[0, 1:-1] = (
        (-3.0*Y[0, 1:-1] + 4.0*Y[1, 1:-1] - Y[2, 1:-1]) / 2.0 / dXi
        )
    XEta[0, 1:-1] = (X[0, 2:]-X[0, 0:-2]) / 2.0 / dEta
    YEta[0, 1:-1] = (Y[0, 2:]-Y[0, 0:-2]) / 2.0 / dEta
    # At boundary X = L
    XXi[-1, 1:-1] = (
        (3.0*X[-1, 1:-1] - 4.0*X[-2, 1:-1] + X[-3, 1:-1]) / 2.0 / dXi
        )
    YXi[-1, 1:-1] = (
        (3.0*Y[-1, 1:-1] - 4.0*Y[-2, 1:-1] + Y[-3, 1:-1]) / 2.0 / dXi
        )
    XEta[-1, 1:-1] = (X[-1, 2:]-X[-1, 0:-2]) / 2.0 / dEta
    YEta[-1, 1:-1] = (Y[-1, 2:]-Y[-1, 0:-2]) / 2.0 / dEta
    # At boundary Y = 0.0
    XXi[1:-1, 0] = (X[2:, 0]-X[0:-2, 0]) / 2.0 / dXi
    YXi[1:-1, 0] = (Y[2:, 0]-Y[0:-2, 0]) / 2.0 / dXi
    XEta[1:-1, 0] = (
        (-3.0*X[1:-1, 0] + 4.0*X[1:-1, 1] - X[1:-1, 2]) / 2.0 / dEta
        )
    YEta[1:-1, 0] = (
        (-3.0*Y[1:-1, 0] + 4.0*Y[1:-1, 1] - Y[1:-1, 2]) / 2.0 / dEta
        )
    # At boundary Y = H
    XXi[1:-1, -1] = (X[2:, -1]-X[0:-2, -1]) / 2.0 / dXi
    YXi[1:-1, -1] = (Y[2:, -1]-Y[0:-2, -1]) / 2.0 / dXi

    XEta[1:-1, -1] = (
        (3.0*X[1:-1, -1] - 4.0*X[1:-1, -2] + X[1:-1, -3]) / 2.0 / dEta
        )
    YEta[1:-1, -1] = (
        (3.0*Y[1:-1, -1] - 4.0*Y[1:-1, -2] + Y[1:-1, -3]) / 2.0 / dEta
        )
    # At vertices
    # X=0.0, Y=0.0
    XXi[0, 0] = (-3.0*X[0, 0] + 4.0*X[1, 0] - X[2, 0]) / 2.0 / dXi
    YXi[0, 0] = (-3.0*Y[0, 0] + 4.0*Y[1, 0] - Y[2, 0]) / 2.0 / dXi
    XEta[0, 0] = (-3.0*X[0, 0] + 4.0*X[0, 1] - X[0, 2]) / 2.0 / dEta
    YEta[0, 0] = (-3.0*Y[0, 0] + 4.0*Y[0, 1] - Y[0, 2]) / 2.0 / dEta
    # X=L, Y=0.0
    XXi[-1, 0] = (3.0*X[-1, 0] - 4.0*X[-2, 0] + X[-3, 0]) / 2.0 / dXi
    YXi[-1, 0] = (3.0*Y[-1, 0] - 4.0*Y[-2, 0] + Y[-3, 0]) / 2.0 / dXi
    XEta[-1, 0] = (-3.0*X[-1, 0] + 4.0*X[-1, 1] - X[-1, 2]) / 2.0 / dEta
    YEta[-1, 0] = (-3.0*Y[-1, 0] + 4.0*Y[-1, 1] - Y[-1, 2]) / 2.0 / dEta
    # X=0.0, Y=H
    XXi[0, -1] = (-3.0*X[0, -1] + 4.0*X[1, -1] - X[2, -1]) / 2.0 / dXi
    YXi[0, -1] = (-3.0*Y[0, -1] + 4.0*Y[1, -1] - Y[2, -1]) / 2.0 / dXi
    XEta[0, -1] = (3.0*X[0, -1] - 4.0*X[0, -2] + X[0, -3]) / 2.0 / dEta
    YEta[0, -1] = (3.0*Y[0, -1] - 4.0*Y[0, -2] + Y[0, -3]) / 2.0 / dEta
    # X=L, Y=H
    XXi[-1, -1] = (3.0*X[-1, -1] - 4.0*X[-2, -1] + X[-3, -1]) / 2.0 / dXi
    YXi[-1, -1] = (3.0*Y[-1, -1] - 4.0*Y[-2, -1] + Y[-3, -1]) / 2.0 / dXi
    XEta[-1, -1] = (3.0*X[-1, -1] - 4.0*X[-1, -2] + X[-1, -3]) / 2.0 / dEta
    YEta[-1, -1] = (3.0*Y[-1, -1] - 4.0*Y[-1, -2] + Y[-1, -3]) / 2.0 / dEta
    # Evaluate metrics and Jacobian
    for i in range(0, iM):
        for j in range(0, jM):
            JJ[i][j] = 1.0 / (XXi[i][j]*YEta[i][j] - YXi[i][j]*XEta[i][j])
            XiX[i][j] = JJ[i][j] * YEta[i][j]
            XiY[i][j] = -JJ[i][j] * XEta[i][j]
            EtaX[i][j] = -JJ[i][j] * YXi[i][j]
            EtaY[i][j] = JJ[i][j] * XXi[i][j]
    return XiX, XiY, EtaX, EtaY, JJ
