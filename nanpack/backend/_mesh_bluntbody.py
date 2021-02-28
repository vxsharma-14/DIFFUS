import math
from ._util_clustering import _clustering_parameter


def _blunt_body_cone(length, cone_angle, cone_radius, major, minor,
                     Xi, Eta, iM, jM, clust_opt, beta,
                     i1_loc="below-stagnation"):
    """Return mesh for the axis-symmetric blunt cone at zero degree AoA."""
    x = Xi.copy()
    y = Xi.copy()

    # Calculate cone dimensions
    cone_ang_degree = cone_angle*math.pi/180.0
    # cone_sh = slant height of the cone (straight line segment)
    cone_sh = (length-cone_radius) / math.cos(cone_ang_degree)
    # cone_cir = circumference of the cone front curvature
    # (quarter circle segment)
    cone_cir = math.pi*cone_radius/2
    # height = vertical distance between cone trailing edge and x-axis
    # height = cone_radius + cone_sh*math.sin(cone_ang_degree)
    # del_i= step size along the circumference of the blunt body
    del_i = (cone_cir+cone_sh) / (iM-1)
    # i1 = grid point location of the leading edge circle
    i1 = int(cone_cir/del_i)
    # angle of each grid lines on blunt cone quarter circle
    # cthet = [i*(math.pi/2.0/(i1-1)) for i in range(0, i1)]
    cthet = _calculate_angle_theta(i1_loc, i1, geo="cone")
    # distance from coordinate origin to center of blunt cone
    # quarter circle on x-axis
    o_dist = major - length + cone_radius
    # Define grid points on the conical blunt body surface
    for i in range(0, i1):
        k = i1 - i - 1
        x[i][0] = o_dist - cone_radius*math.cos(cthet[i])
        y[i][0] = cone_radius*math.sin(cthet[i])
    dh = cone_sh / (iM-i1)
    for i in range(i1, iM):
        k = i - i1 + 1
        x[i][0] = o_dist + k*dh
        y[i][0] = cone_radius + k*dh*math.sin(cone_ang_degree)
    thet = _calculate_angle_theta(i1_loc, iM, geo="ellip")
    for i in range(0, iM):
        k = iM - i - 1
        # Equation after 9-56 in CFD Vol.1 by Dr. Hoffmann
        ro1 = math.sin(thet[k])*math.sin(thet[k]) / (major*major)
        ro2 = math.cos(thet[k])*math.cos(thet[k]) / (minor*minor)
        r_o = 1.0 / math.sqrt(ro1+ro2)
        # Define grid points on the outer ellipse surface
        x[i][-1] = major - r_o*math.sin(thet[k])
        y[i][-1] = r_o*math.cos(thet[k])

    x, y = _calculate_interior_mesh(x, y, Eta, iM, jM, beta, clust_opt)

    return x, y


def _blunt_body_ellipt(major1, major2, minor1, minor2,
                       Xi, Eta, iM, jM, clust_opt, beta,
                       i1_loc="below-stagnation"):
    """Return mesh for the axis-symmetric blunt cone at zero degree AoA.

    This formulation is for an elliptical body not a conical body.
    make a new function for the cone and change description.
    """
    x = Xi.copy()
    y = Xi.copy()

    thet = _calculate_angle_theta(i1_loc, iM, geo="ellip")

    for i in range(0, iM):
        # Equation after 9-56 in CFD Vol.1 by Dr. Hoffmann
        k = iM - i - 1
        ri1 = math.sin(thet[k])*math.sin(thet[k]) / (major2*major2)
        ri2 = math.cos(thet[k])*math.cos(thet[k]) / (minor2*minor2)
        r_i = 1.0 / math.sqrt(ri1+ri2)
        # Define grid points on the elliptical blunt body surface
        x[i][0] = major1 - r_i*math.sin(thet[k])
        y[i][0] = r_i*math.cos(thet[k])

        ro1 = math.sin(thet[k])*math.sin(thet[k]) / (major1*major1)
        ro2 = math.cos(thet[k])*math.cos(thet[k]) / (minor1*minor1)
        r_o = 1.0 / math.sqrt(ro1+ro2)
        # Define grid points on the outer ellipse surface
        x[i][-1] = major1 - r_o*math.sin(thet[k])
        y[i][-1] = r_o*math.cos(thet[k])

    x, y = _calculate_interior_mesh(x, y, Eta, iM, jM, beta, clust_opt)

    return x, y


def _calculate_interior_mesh(x, y, Eta, iM, jM, beta, clust_opt):
    S = Eta.copy()
    # Calculate the distance, delI between j = jM and j = 1
    # Calculate the angle alpI of each line i from j = 1 to j = j jM
    delI = []
    alpha = []
    for i in range(0, iM):
        delxx = x[i, 0] - x[i, -1]
        delyy = y[i, -1] - y[i, 0]
        delI.append(math.sqrt(delxx*delxx + delyy*delyy))
        alpha.append(math.asin(delyy/delI[i]))

    # Calculate interior grid point locations using
    # clustering option = True or False.
    for i in range(0, iM):
        for j in range(0, jM):
            if clust_opt is True:
                S[i][j] = _clustering_parameter(Eta[i][j], Eta[i, -1],
                                                beta, delI[i])
            else:
                S[i][j] = delI[i]*(j)/(jM-1)
    for j in range(1, jM-1):
        for i in range(0, iM):
            x[i][j] = x[i][0] - S[i][j]*math.cos(alpha[i])
            y[i][j] = y[i][0] + S[i][j]*math.sin(alpha[i])

    return x, y


def _calculate_angle_theta(i1loc, im, geo):
    """Return the value of angle theta for grid lines i = 1..im."""
    if geo == "ellip":
        d_alp = _calculate_angle_alpha_ellip(i1loc, im)
    elif geo == "cone":
        d_alp = _calculate_angle_alpha_cone(i1loc, im)
    thet = [0.0 for i in range(0, im)]
    # Calculate angle theta along all grid lines i
    for i in range(0, im):
        # t = im - i - 1
        thet[i] = i*d_alp
    return thet


def _calculate_angle_alpha_ellip(i1_location, im):
    """Return angle step size- dalpha.

    This formulation places i = 1 below or along the line of body symmetry.

    if i1_placement = BELOW STAGNATION LINE, i=1 is placed below the
                                             stagnation line such that
                                             i = 1 and i = 2 are
                                             symmetrical about the
                                             stagnation line (along x-axis)
    if i1_placement = ALONG STAGNATION LINE, i=1 is placed along the
                                             stagnation line.

    Place i = 1 below stagnation line to avoid difficulties in
    solving fluid equations.
    """
    if i1_location.lower() == "below-stagnation":
        return (-math.pi / (3-2*im))
    elif i1_location.lower() == "along-stagnation":
        return (math.pi / (2*im-2))


def _calculate_angle_alpha_cone(i1_location, im):
    """Return angle step size- dalpha.

    This formulation places i = 1 below or along the line of body symmetry.

    if i1_location = BELOW STAGNATION LINE, i=1 is placed below the
                                             stagnation line such that
                                             i = 1 and i = 2 are
                                             symmetrical about the
                                             stagnation line (along x-axis)
    if i1_location = ALONG STAGNATION LINE, i=1 is placed along the
                                             stagnation line.

    Place i = 1 below stagnation line to avoid difficulties in
    solving fluid equations.
    """
    if i1_location.lower() == "below-stagnation":
        return (-math.pi / 2.0 / (1.5-im))
    elif i1_location.lower() == "along-stagnation":
        return (math.pi / 2.0 / (im-1))
