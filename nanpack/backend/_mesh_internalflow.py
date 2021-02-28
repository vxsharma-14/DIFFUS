from ._util_clustering import _calculate_beta1, _calculate_eta1


def _cavity(A, B, x, y, iM, jM, Xi, Eta, alphax, alphay, betax, betay):
    betax1 = _calculate_beta1(betax)
    betay1 = _calculate_beta1(betay)

    for i in range(0, iM):
        # Grid points along wall Y = 0.0
        xi1 = _calculate_eta1(Xi[i, 0], Xi[-1, 0], False)
        x[i][0] = _get_clust_wall(xi1, alphax, betax, betax1, A)
        y[i][0] = 0.0
        # Grid points along wall Y = B
        xi1 = _calculate_eta1(Xi[i, -1], Xi[-1, -1], False)
        x[i][-1] = _get_clust_wall(xi1, alphax, betax, betax1, A)
        y[i][-1] = B

    for j in range(0, jM):
        # Grid points along wall X = 0.0
        x[0][j] = 0.0
        eta1 = _calculate_eta1(Eta[0, j], Eta[0, -1], False)
        y[0][j] = _get_clust_wall(eta1, alphay, betay, betay1, B)
        # Grid points along wall X = A
        x[-1][j] = A
        eta1 = _calculate_eta1(Eta[-1, j], Eta[-1, -1], False)
        y[-1][j] = _get_clust_wall(eta1, alphay, betay, betay1, B)

    for i in range(1, iM-1):
        for j in range(1, jM-1):
            x[i][j] = x[i][0]
            y[i][j] = y[0][j]

    return x, y


def _duct(A, B, x, y, dX, dY, Xi, Eta, iM, jM, alphax, alphay, betax,
          betay):
    if betax is not None:
        beta1 = _calculate_beta1(betax)
        for i in range(0, iM):
            # Grid points along wall Y = 0.0
            xi1 = _calculate_eta1(Xi[i, 0], Xi[-1, 0], False)
            x[i][0] = _get_clust_wall(xi1, alphax, betax, beta1, A)
            y[i][0] = 0.0
            # Grid points along wall Y = B
            xi1 = _calculate_eta1(Xi[i, -1], Xi[-1, -1], False)
            x[i][-1] = _get_clust_wall(xi1, alphax, betax, beta1, A)
            y[i][-1] = B

            for j in range(1, jM-1):
                x[i][j] = x[i][0]
                y[i][j] = y[i][j-1] + dY

    elif betay is not None:
        beta1 = _calculate_beta1(betay)
        for j in range(0, jM):
            # Grid points along wall X = 0.0
            x[0][j] = 0.0
            eta1 = _calculate_eta1(Eta[0, j], Eta[0, -1], False)
            y[0][j] = _get_clust_wall(eta1, alphay, betay, beta1, B)
            # Grid points along wall X = A
            x[-1][j] = A
            eta1 = _calculate_eta1(Eta[-1, j], Eta[-1, -1], False)
            y[-1][j] = _get_clust_wall(eta1, alphay, betay, beta1, B)

            for i in range(1, iM-1):
                x[i][j] = x[i-1][j] + dX
                y[i][j] = y[0][j]

    return x, y


def _get_clust_wall(metric1, alpha, beta, beta1, dom_size):
    A = dom_size
    met_e = (metric1-alpha) / (1.0-alpha)
    axis = A * (((2.0*alpha+beta)*(beta1**met_e)
                + 2.0*alpha - beta)
                / ((2.0*alpha+1.0)*(1.0+(beta1**met_e))))

    return axis
