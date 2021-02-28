def _clustering_parameter(eta_i_j, eta_i_jm, beta, del_i):
    """Return the clustering function at each grid point i, j."""
    beta1 = _calculate_beta1(beta)
    eta1 = _calculate_eta1(eta_i_j, eta_i_jm)
    S = del_i * (1.0 - beta * (beta1**eta1-1.0)/(beta1**eta1+1.0))
    return S


def _calculate_beta1(beta):
    return (beta+1.0) / (beta-1.0)


def _calculate_eta1(eta_i_j, eta_i_jm, sub_from_1=True):
    if sub_from_1 is False:
        return (eta_i_j/eta_i_jm)
    else:
        return (1.0 - eta_i_j/eta_i_jm)
