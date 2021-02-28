"""Not a public module in version 1.0.0a4."""
#   ***********************************************************************
#
#   FILE         dimensionalize.py
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


class NonDimensionalize:
    """A class to non-dimensionalize various numerical parameters.

    This is not complete or tested for accuracy.
    """

    def __init__(self):
        """Class construcotr."""

    def ndGrid(self, X, L):
        """Non-dimensionalize x point locations along X axis."""
        self.xstar[:] = X[:]/L
        return self.xstar

    def ndVelU(self, U, L, nu):
        """Non-dimensionalize velocity U."""
        shapeU = U.shape
        if len(shapeU) == 1:
            im = shapeU
            self.ustar = [U[i]*L/nu for i in range(im)]
        elif len(shapeU) == 2:
            im, jm = shapeU
            self.ustar = [[U[i][j]*L/nu for j in range(jm)]
                          for i in range(im)]
        return self.ustar

    def ndTime(self, t, L, nu):
        """Non-dimensionalize time t."""
        L2 = L*L
        self.tstar = nu*t/L2
        return self.tstar


class Dimensionalize:
    """A class to dimensionalize various numerical parameters.

    This is not complete or tested for accuracy.
    """

    def __init__(self):
        """Class constructor."""

    def dimGrid(self, xstar, L):
        """Dimensionalize x point locations along X axis."""
        self.x[:] = xstar[:]*L
        return self.x

    def dimVelU(self, ustar, L, nu):
        """Dimensionalize velocity U."""
        import numpy as np
        shapeU = ustar.shape
        if len(shapeU) == 1:
            iM, = shapeU
            self.u = np.zeros((iM), dtype="float")
            for i in range(0, iM):
                self.u[i] = ustar[i]*nu/L
        elif len(shapeU) == 2:
            iM, jM = shapeU
            self.u = np.zeros((iM), dtype="float")
            for i in range(0, iM):
                for j in range(0, jM):
                    self.u[i][j] = ustar[i][j]*nu/L
        return self.u

    def dimTime(self, tstar, L, nu):
        """Non-dimensionalize time t."""
        L2 = L*L
        self.t = tstar*L2/nu
        return self.t
