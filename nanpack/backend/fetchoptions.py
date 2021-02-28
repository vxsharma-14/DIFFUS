"""Not a public module."""

#   ***********************************************************************
#
#   FILE         fetchoptions.py
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


class _FetchOptions:
    """A class to fetch allowed input options for the user."""

    def __init__(self):
        """Class constructor."""

    def UnitOptions(self):
        """Return a list of allowed inputs for UNIT key in config file."""
        self.options = ["SI",
                        "BRITISH"
                        ]
        return self.options

    def StateOptions(self):
        """Return a list of allowed inputs for STATE key in config file."""
        self.options = ["STEADY",
                        "TRANSIENT"
                        ]
        return self.options

    def ModelOptions(self):
        """Return a list of allowed inputs for MODEL key in config file."""
        self.options = ["DIFFUSION",
                        "LAPLACE",
                        "FO_WAVE",
                        "INV_BURGERS",
                        "VISC_BURGERS"]
        return self.options

    def SolverOptions(self):
        """Return a list of allowed inputs for SOLV key in config file."""

    def DimensionOptions(self):
        """Return a list of allowed inputs for DIM key in config file."""
        self.options = ["1D",
                        "2D"]
        return self.options

    def TVDLimiterFunctionOptions(self):
        """Return a list of allowed inputs for Limiter functions.

        The LIMITERFUNC argument is required in the call to function
        SecondOrderTVD().
        """
        self.options = ["Harten-Yee-Upwind",
                        "Modified-Harten-Yee-Upwind",
                        "Roe-Sweby-Upwind",
                        "Davis-Yee-Symmetric"]
        return self.options
