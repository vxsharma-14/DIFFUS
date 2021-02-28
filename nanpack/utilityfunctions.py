"""Not a public function in version 1.0.0a4."""
#   ***********************************************************************
#
#   FILE         utilityfunctions.py
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


def EntropyCorrectionFunction(Alpha, Eps):
    """Return the value of the entropy correction term.

    Solve equation 6-127 or Equation 6-121) defined in CFD Vol. 1
    by Hoffmann.

    Call signature:

        EntropyCorrectionFunction(Alpha, Eps)

    Parameters
    ----------
    Alpha : float

        The entroopy correction is a function of alpha to prevent the
        difficulty that may appear when alpha becomes zero in TVD
        formulation.

    Eps : float

        A positive constant similar to the damping coefficient. Its value
        must be selected within the range 0f 0.0 to 0.125.

    Returns
    -------
    si : float

        Entropy correction term.
    """
    AlphaAbs = abs(Alpha)
    if AlphaAbs >= Eps:
        si = AlphaAbs
    else:
        si = (Alpha*Alpha + Eps*Eps)/(2*Eps)

    return si


def _CalcAlpha(E1, E2, dU, U1, U2):
    """Return the results of Equation 6-128 in CFD Vol.1 by Hoffmann."""
    return ((E1 - E2)/dU if dU != 0 else 0.5*(U1 + U2))


def _CalcUi(U2, U1):
    """Return the backward differencing result of the dependent var."""
    return (U2 - U1)


def _CalcE(U):
    """Return the result of the non-linear variable E = 0.5U**2."""
    return (0.5*U*U)
