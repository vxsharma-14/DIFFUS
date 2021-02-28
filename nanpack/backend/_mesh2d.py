"""."""
from ._util_meshkwargs import _unpack_clust_kwargs
from ._util_meshkwargs import _unpack_geometric_kwargs
from ._mesh_internalflow import _duct, _cavity
from ._mesh_flatplate import _flateplate
from ._mesh_otype import _ogrid_airfoil, _ogrid_cylinder
from ._mesh_bluntbody import _blunt_body_cone, _blunt_body_ellipt
from .exceptions import MeshingInputError


def _meshing_func(geo_type, clust_opt, Xi, Eta, dX, dY, iM, jM,
                  **mesh_kw):
    """Return a meshing function requested by the user."""
    geo_temp = {
        "blunt-body-cone": _bluntcone_mesh,
        "blunt-body-ellipse": _bluntellip_mesh,
        # "C-grid": _cgrid_grid,
        "cavity": _cavity_mesh,
        "duct": _duct_mesh,
        "flat-plate": _flatplate_mesh,
        "o-grid-airfoil": _o_airfoil_mesh,
        "o-grid-cylinder": _o_cyl_mesh,
        # "Wind-tunnel-airfoil": _wt_airfoil_grid,
        # "Wind-tunnel-cylinder": _wt_cyl_grid,
        }
    SelectedGeometry = geo_temp.get(geo_type)

    return SelectedGeometry(clust_opt, Xi, Eta, dX, dY, iM, jM,
                            **mesh_kw)


def _bluntcone_mesh(clust_opt, Xi, Eta, dX, dY, iM, jM, **mesh_kw):
    # clust_opt option may be ON or OFF in this type of mesh.
    _, beta = _unpack_clust_kwargs(**mesh_kw)
    beta = beta["beta"]
    if beta is None and clust_opt is True:
        raise MeshingInputError("Beta", "clustering value not provided")
    geom = _unpack_geometric_kwargs(**mesh_kw)
    length = geom["length"]
    angle = geom["cangle"]
    radius = geom["cradius"]
    major = geom["major_out"]
    minor = geom["minor_out"]
    i1loc = geom["i1_location"]
    x, y = _blunt_body_cone(length, angle, radius, major, minor,
                            Xi, Eta, iM, jM, clust_opt, beta, i1loc)
    return x, y


def _bluntellip_mesh(clust_opt, Xi, Eta, dX, dY, iM, jM, **mesh_kw):
    # clust_opt option may be ON or OFF in this type of mesh.
    _, beta = _unpack_clust_kwargs(**mesh_kw)
    beta = beta["beta"]
    if beta is None and clust_opt is True:
        raise MeshingInputError("Beta", "clustering value not provided")
    geom = _unpack_geometric_kwargs(**mesh_kw)
    major1 = geom["major_out"]
    major2 = geom["major_in"]
    minor1 = geom["minor_out"]
    minor2 = geom["minor_in"]
    i1loc = geom["i1_location"]
    x, y = _blunt_body_ellipt(major1, major2, minor1, minor2, Xi, Eta,
                              iM, jM, clust_opt, beta, i1loc)
    return x, y


def _cgrid_mesh(clust_opt, Xi, Eta, dX, dY, iM, jM, **mesh_kw):
    # clust_opt option be ON or OFF in this type of mesh.
    return


def _duct_mesh(clust_opt, Xi, Eta, dX, dY, iM, jM, **mesh_kw):
    """Return a mesh in a duct with clustering along X or Y walls."""
    alpha, beta = _unpack_clust_kwargs(**mesh_kw)
    alpha_x = alpha["alphaX"]
    alpha_y = alpha["alphaY"]
    beta_x = beta["betaX"]
    beta_y = beta["betaY"]
    # clust_opt option must always be ON in this type of mesh.
    if beta_x is None and beta_y is None:
        raise MeshingInputError("BetaX or BetaY", "value not entered")
    elif beta_x is not None and beta_y is not None:
        raise MeshingInputError("BetaX, BetaY", "only 1 must be entered.")
    geom = _unpack_geometric_kwargs(**mesh_kw)
    A = geom["length"]
    B = geom["height"]
    if A is None:
        A = (iM-1) * dX
    if B is None:
        B = (jM-1) * dY
    x = Xi.copy()
    y = Eta.copy()
    x, y = _duct(A, B, x, y, dX, dY, Xi, Eta, iM, jM, alpha_x, alpha_y,
                 beta_x, beta_y)

    return x, y


def _cavity_mesh(clust_opt, Xi, Eta, dX, dY, iM, jM, **mesh_kw):
    """Return a mesh in a cavity with clustering at all walls."""
    alpha, beta = _unpack_clust_kwargs(**mesh_kw)
    alpha_x = alpha["alphaX"]
    alpha_y = alpha["alphaY"]
    beta_x = beta["betaX"]
    beta_y = beta["betaY"]
    # clust_opt option must always be ON in this type of mesh.
    if beta_x is None or beta_y is None:
        raise MeshingInputError("BetaX, BetaY", "both must be entered")
    geom = _unpack_geometric_kwargs(**mesh_kw)
    A = geom["length"]
    B = geom["height"]
    if A is None:
        A = (iM-1) * dX
    if B is None:
        B = (jM-1) * dY
    x = Xi.copy()
    y = Eta.copy()
    x, y = _cavity(A, B, x, y, iM, jM, Xi, Eta, alpha_x, alpha_y,
                   beta_x, beta_y)

    return x, y


def _flatplate_mesh(clust_opt, Xi, Eta, dX, dY, iM, jM, **mesh_kw):
    # clust_opt option must always be ON in this type of mesh.
    """Return a mesh in a cavity with clustering at all walls."""
    _, beta = _unpack_clust_kwargs(**mesh_kw)
    beta = beta["beta"]
    # clust_opt option must always be ON in this type of mesh.
    if beta is None:
        raise MeshingInputError("Beta", "clustering value not provided")
    geom = _unpack_geometric_kwargs(**mesh_kw)
    A = geom["length"]
    B = geom["height"]
    if A is None:
        A = (iM-1) * dX
    if B is None:
        B = (jM-1) * dY
    x = Xi.copy()
    y = Eta.copy()
    x, y = _flateplate(A, B, x, y, dX, dY, iM, jM, Xi, Eta, beta)

    return x, y


def _o_cyl_mesh(clust_opt, Xi, Eta, dX, dY, iM, jM, **mesh_kw):
    # clust_opt option may be ON or OFF in this type of mesh.
    _, beta = _unpack_clust_kwargs(**mesh_kw)
    beta = beta["beta"]
    if beta is None and clust_opt is True:
        raise MeshingInputError("Beta", "clustering value not provided")
    geom = _unpack_geometric_kwargs(**mesh_kw)
    orad = geom["radius_out"]
    irad = geom["radius_in"]
    x = Xi.copy()
    y = Xi.copy()
    x, y = _ogrid_cylinder(x, y, orad, irad, Xi, Eta, iM, jM,
                           clust_opt, beta)

    return x, y


def _o_airfoil_mesh(clust_opt, Xi, Eta, dX, dY, iM, jM, **mesh_kw):
    # clust_opt option may be ON or OFF in this type of mesh.
    _, beta = _unpack_clust_kwargs(**mesh_kw)
    beta = beta["beta"]
    if beta is None and clust_opt is True:
        raise MeshingInputError("Beta", "clustering value not provided")
    geom = _unpack_geometric_kwargs(**mesh_kw)
    ch = geom["chord"]
    thick = geom["thickness"]
    rad = geom["radius_out"]
    x = Xi.copy()
    y = Xi.copy()
    x, y = _ogrid_airfoil(x, y, rad, ch, thick, Xi, Eta, iM, jM,
                          clust_opt, beta)

    return x, y


def _wt_airfoil_mesh():
    return


def _wt_cyl_mesh():
    return
