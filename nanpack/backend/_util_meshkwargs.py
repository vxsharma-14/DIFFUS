"""Not a public module."""


def _unpack_clust_kwargs(**grid_kw):
    alpha = grid_kw.get("Alpha", None)
    beta = grid_kw.get("Beta", None)
    alphaX = grid_kw.get("AlphaX", 0.5)
    alphaY = grid_kw.get("AlphaY", 0.5)
    betaX = grid_kw.get("BetaX", None)
    betaY = grid_kw.get("BetaY", None)
    alpha = _pack_alpha_into_dict(alphaX, alphaY, alpha)
    beta = _pack_beta_into_dict(betaX, betaY, beta)
    return alpha, beta


def _unpack_geometric_kwargs(**grid_kw):
    lgt = grid_kw.get("Length", None)
    hgt = grid_kw.get("Height", None)
    chd = grid_kw.get("Chord", 1.0)
    th = grid_kw.get("Thickness", 0.2)
    orad = grid_kw.get("outRad", 3.0)
    irad = grid_kw.get("inRad", 1.0)
    cang = grid_kw.get("cAngle", 5.0)
    crad = grid_kw.get("cRadius", 1.0)
    omajor = grid_kw.get("outMajor", 4.0)
    ominor = grid_kw.get("outMinor", 3.5)
    imajor = grid_kw.get("inMajor", 3.0)
    iminor = grid_kw.get("inMinor", 2.0)
    i1loc = grid_kw.get("i1Location", 2.0)
    geo = _pack_geom_into_dict(lgt, hgt, chd, th, orad, irad, cang,
                               crad, omajor, ominor, imajor, iminor, i1loc)
    return geo


def _pack_alpha_into_dict(alphaX, alphaY, alpha):
    """Return a dictionary for alpha."""
    alp = {
        "alphaX": alphaX,
        "alphaY": alphaY,
        "alpha": alpha
        }

    return alp


def _pack_beta_into_dict(betaX, betaY, beta):
    """Return a dictionary for beta."""
    bet = {
        "betaX": betaX,
        "betaY": betaY,
        "beta": beta
        }

    return bet


def _pack_geom_into_dict(lgt, hgt, ch, th, orad, irad, cang, crad, omaj,
                         ominr, imaj, iminr, i1loc):
    """Return a dictionary for geometric data."""
    geo = {
        "length": lgt,
        "height": hgt,
        "chord": ch,
        "thickness": th,
        "radius_out": orad,
        "radius_in": irad,
        "cangle": cang,
        "cradius": crad,
        "major_out": omaj,
        "minor_out": ominr,
        "major_in": imaj,
        "minor_in": iminr,
        "i1_location": i1loc
        }
    return geo
