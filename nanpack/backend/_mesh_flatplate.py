from ._util_clustering import _calculate_beta1, _calculate_eta1


def _flateplate(A, B, x, y, dX, dY, iM, jM, Xi, Eta, beta):
    beta1 = _calculate_beta1(beta)
    for j in range(0, jM):
        # Grid points along wall X = 0.0
        eta1 = _calculate_eta1(Eta[0, j], Eta[0, -1])
        x[0][j] = 0.0
        y[0][j] = B * (((beta+1.0) - (beta-1.0)*(beta1**eta1))
                       / ((beta1**eta1)+1.0))
        # Grid points along wall X = A
        eta1 = _calculate_eta1(Eta[-1, j], Eta[-1, -1])
        x[-1][j] = A
        y[-1][j] = B * (((beta+1) - (beta-1)*(beta1**eta1))
                        / ((beta1**eta1)+1.0))

        for i in range(1, iM-1):
            x[i][j] = x[i-1][j] + dX
            y[i][j] = y[0][j]

    return x, y
