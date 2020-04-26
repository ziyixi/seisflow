"""
make_perturbed_cmtsolution.py: make perturbed cmt solution using the frechet information.
"""
import numpy as np
import pyproj


def add_src_frechet(src_frechet, cmtsolution, max_dxs_ratio, fix_location=False, fix_focal=False):
    data = src_frechet
    # convert from dyne*cm to N*m
    dchi_dmt = np.array([
        [data[0], data[3], data[4]],
        [data[3], data[1], data[5]],
        [data[4], data[5], data[2]]
    ]) * 1e7
    # ! note here we have a bug, the output of specfem should be de,dn,-dz, to change to dr,dtheta,dphi, we have to
    # ! change dp,-dt,-dr to dr,dt,dp (-ddep,-dlat,dlon)
    dchi_dxs = np.array([
        -data[8],
        -data[7],
        data[6]
    ])
    # do the normalization as in tao's sem_utils
    # the final result will be the same regardless the factor
    cmt_tensor = cmtsolution.focal_mechanisms[0].moment_tensor.tensor
    mt = np.array([[cmt_tensor.m_rr, cmt_tensor.m_rt, cmt_tensor.m_rp],
                   [cmt_tensor.m_rt, cmt_tensor.m_tt, cmt_tensor.m_tp],
                   [cmt_tensor.m_rp, cmt_tensor.m_tp, cmt_tensor.m_pp]])
    m0 = (0.5*np.sum(mt**2))**0.5
    R_earth = 6371000.0
    dchi_dxs_ratio = R_earth * dchi_dxs
    dchi_dmt_ratio = m0 * dchi_dmt

    # ====== scale CMT gradient
    if (np.sum(dchi_dxs_ratio ** 2) == 0):
        scale_factor_xs = 0.0
    else:
        scale_factor_xs = max_dxs_ratio / (np.sum(dchi_dxs_ratio ** 2)) ** 0.5
    if (np.sum(dchi_dmt_ratio ** 2) == 0):
        scale_factor_mt = 0.0
    else:
        scale_factor_mt = max_dxs_ratio/(np.sum(dchi_dmt_ratio**2))**0.5
    dxs_ratio = scale_factor_xs * dchi_dxs_ratio
    dmt_ratio = scale_factor_mt * dchi_dmt_ratio
    print(dxs_ratio, dmt_ratio, scale_factor_xs, scale_factor_mt)
    dxs = R_earth * dxs_ratio
    dmt = m0 * dmt_ratio

    # * add to the raw CMTSOLUTION
    # firstly we have to rely on x,y,z to convert the coordinate (not a sphere)
    ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    lat = cmtsolution.preferred_origin().latitude
    lon = cmtsolution.preferred_origin().longitude
    alt = -cmtsolution.preferred_origin().depth
    x, y, z = pyproj.transform(lla, ecef, lon, lat, alt)
    r = (x**2 + y**2 + z**2)**0.5
    # get rotation matrix
    theta = np.arccos(z/r)
    phi = np.arctan2(y, x)
    sthe = np.sin(theta)
    cthe = np.cos(theta)
    sphi = np.sin(phi)
    cphi = np.cos(phi)
    # coordinate transformation matrix (r,theta,phi) to (x,y,z)
    a = np.array(
        [[sthe*cphi, cthe*cphi, -1.0*sphi],
         [sthe*sphi, cthe*sphi,      cphi],
            [cthe, -1.0*sthe,      0.0]])
    # dr,dtheta,dphi to dx,dy,dz
    dxs_xyz = np.dot(a, dxs)
    x += dxs_xyz[0]
    y += dxs_xyz[1]
    z += dxs_xyz[2]
    lon, lat, alt = pyproj.transform(ecef, lla, x, y, z)
    # # add dmt
    # mt += dmt
    # we have to get mt at the new position
    mt_xyz = np.dot(np.dot(a, mt), np.transpose(a))
    dmt_xyz = np.dot(np.dot(a, dmt), np.transpose(a))
    mt_xyz += dmt_xyz
    # get new a at new x,y,z
    r = (x**2 + y**2 + z**2)**0.5
    theta = np.arccos(z/r)
    phi = np.arctan2(y, x)
    sthe = np.sin(theta)
    cthe = np.cos(theta)
    sphi = np.sin(phi)
    cphi = np.cos(phi)
    a = np.array(
        [[sthe*cphi, cthe*cphi, -1.0*sphi],
         [sthe*sphi, cthe*sphi,      cphi],
            [cthe, -1.0*sthe,      0.0]])
    # convert back to mt
    mt = np.dot(np.dot(np.transpose(a), mt_xyz), a)

    # write to the new CMTSOLUTION
    cmtsolution_new = cmtsolution.copy()
    if(not fix_location):
        cmtsolution_new.preferred_origin().latitude = lat
        cmtsolution_new.preferred_origin().longitude = lon
        cmtsolution_new.preferred_origin().depth = -alt
    if(not fix_focal):
        cmtsolution_new.focal_mechanisms[0].moment_tensor.tensor.m_rr = mt[0, 0]
        cmtsolution_new.focal_mechanisms[0].moment_tensor.tensor.m_tt = mt[1, 1]
        cmtsolution_new.focal_mechanisms[0].moment_tensor.tensor.m_pp = mt[2, 2]
        cmtsolution_new.focal_mechanisms[0].moment_tensor.tensor.m_rt = mt[1, 0]
        cmtsolution_new.focal_mechanisms[0].moment_tensor.tensor.m_rp = mt[2, 0]
        cmtsolution_new.focal_mechanisms[0].moment_tensor.tensor.m_tp = mt[2, 1]

    return cmtsolution_new
