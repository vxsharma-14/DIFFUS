"""Not a public module."""
#   ***********************************************************************
#
#   FILE         exceptions.py
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


class NumericalMethodError(Exception):
    """Raise exception for wrong numerical method and model combination."""

    def __init__(self, method, model):
        self.message = f"Numerical Method '{method}' invalid for the\
 {model} equation."
        super().__init__(self.message)


class DimensionError(Exception):
    """Raise exception for wrong dimension and model/method combination."""

    def __init__(self, dim, model, scheme):
        self.message = f"Dimension '{dim}' of dependent variable is\
 invalid. In this version, {scheme} method can solve {dim} {model}\
 equation."
        super().__init__(self.message)


class InputFileError(Exception):
    """Raise exception when invalid/no file path is provided.

    Parameters
    ----------
    which:

        Allowed inputs for this paramter are: "NoInput" or "FileNotFound".

    file: Default=None

        File name entered by the user.
    """

    def __init__(self, which, file=None):
        if which == "NoInput":
            msg = "No input files provided."
        elif which == "FileNotFound":
            msg = f'File "{file}" not found in the directory.'
        self.message = msg
        super().__init__(self.message)


class InvalidValueError(Exception):
    """Raise exception when invalid key value is encountered."""

    def __init__(self, key, value):
        self.message = f"Invalid value '{value}' to key '{key}'\
 in configuration file."
        super().__init__(self.message)


class SectionNotFoundError(Exception):
    """Raise exception when a section is not found in the config files."""

    def __init__(self, section, file):
        self.message = f"Section [{section}] not found\nin file: {file}."
        super().__init__(self.message)


class TVDLimiterInputError(Exception):
    """Raise exception when an invalid Limiter is entered."""

    def __init__(self, lmtr, limfunc):
        self.msg = f"Invalid TVD limiter {lmtr} provided for {limfunc}."
        super().__init__(self.msg)


class TVDLimiterFunctionInputError(Exception):
    """Raise exception when an invalid Limiter function is entered."""

    def __init__(self, limfunc):
        self.msg = f"Invalid TVD limiter function {limfunc} provided in\
 the call to SecondOrderTVD() function."
        super().__init__(self.msg)


class MeshingInputError(Exception):
    """Raise exception for invalid inputs in the  grid functions."""

    def __init__(self, arg, text):
        self.msg = f"{arg} : {text}."
        super().__init__(self.msg)
