#**************************************************************************
#
#   FILE         viscous-burgers-analytical.py
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
#**************************************************************************
def test_analyt():
    '''Test viscous Burgers analytical solution routine.'''
    import nanpack.grid as grid
    import nanpack.benchmark as bm
    import matplotlib.pyplot as plt
    import numpy as np
    x, _ = grid.RectangularGrid(0.2, 91)
    Uanaly = bm.ViscousBurgersSolution(1, x, 0.1)
    ax = plt.subplot()
    plt.plot(x,Uanaly)
    plt.xticks(np.arange(-9.0, 10.0, 3.0))
    plt.yticks(np.arange(round(min(Uanaly)), round(max(Uanaly)+1), 1.0))
    plt.show()

if __name__ == '__main__':
    test_analyt()
