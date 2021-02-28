from .exceptions import MeshingInputError
from ._util_clustering import _calculate_beta1, _calculate_eta1


def _meshing_func(clust_loc, Xi, dX, iM, alpha, beta):
    """Return a meshing function requested by the user."""
    if beta is None:
        raise MeshingInputError("Beta", "clustering value not provided")
    clust_region = {
        "left": _clustering_left_wall,
        "right": _clustering_right_wall,
        "both": _clustering_both_walls,
        "middle": _clustering_mid_axis,
        }
    SelectedFunction = clust_region.get(clust_loc)

    return SelectedFunction(Xi, dX, iM, alpha, beta)


def _clustering_left_wall(Xi, dX, iM, alpha, beta):
    # clust_opt option must always be ON in this type of mesh.
    x = Xi.copy()
    A = dX * (iM-1)
    beta1 = _calculate_beta1(beta)
    for i in range(0, iM):
        # Grid clustering near X = 0.0
        xi1 = _calculate_eta1(Xi[i], Xi[-1])
        x[i] = A * (((beta+1.0) - (beta-1.0)*(beta1**xi1))
                    / ((beta1**xi1)+1.0))
    return x


def _clustering_right_wall(Xi, dX, iM, alpha, beta):
    # clust_opt option must always be ON in this type of mesh.
    x = Xi.copy()
    A = dX * (iM-1)
    beta1 = _calculate_beta1(beta)
    alpha = 0.0
    for i in range(0, iM):
        # Grid clustering near X = A
        xi1 = _calculate_eta1(Xi[i], Xi[-1], False)
        x[i] = _get_clust_wall(xi1, alpha, beta, beta1, A)
    return x


def _clustering_both_walls(Xi, dX, iM, alpha, beta):
    # clust_opt option must always be ON in this type of mesh.
    x = Xi.copy()
    A = dX * (iM-1)
    beta1 = _calculate_beta1(beta)
    for i in range(0, iM):
        # Grid clustering equally distributed between X = 0 and X = A
        xi1 = _calculate_eta1(Xi[i], Xi[-1], False)
        x[i] = _get_clust_wall(xi1, alpha, beta, beta1, A)
    return x


def _clustering_mid_axis(clust_opt, Xi, dX, iM):
    return


def _get_clust_wall(metric1, alpha, beta, beta1, dom_size):
    A = dom_size
    met_e = (metric1-alpha) / (1.0-alpha)
    axis = A * (((2.0*alpha+beta)*(beta1**met_e)
                + 2.0*alpha - beta)
                / ((2.0*alpha+1.0)*(1.0+(beta1**met_e))))

    return axis
