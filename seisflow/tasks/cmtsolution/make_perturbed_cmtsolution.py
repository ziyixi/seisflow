"""
make_perturbed_cmtsolution_kim2011.py: make perturbed cmtsolution based on doi: 10.1111/j.1365-246X.2011.05027.x
"""


import numpy as np
import pyproj


def add_src_frechet(src_frechet, cmtsolution, max_dxs_ratio, fix_location=False, fix_focal=False):
    data = src_frechet
    # convert from dyne*cm to N*m
    # factor 2 is for converting 9 parametes to 6.
    dchi_dmt = np.array([
        [data[0], data[3]/2, data[4]/2],
        [data[3]/2, data[1], data[5]/2],
        [data[4]/2, data[5]/2, data[2]]
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
    lat = cmtsolution.preferred_origin().latitude
    lon = cmtsolution.preferred_origin().longitude
    alt = -cmtsolution.preferred_origin().depth
    # * we need to transfer everything to x,y,z
    ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
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
    # * dr,dtheta,dphi to dx,dy,dz
    dchi_dxs_xyz = np.dot(a, dchi_dxs)
    mt_xyz = np.dot(np.dot(a, mt), np.transpose(a))
    dchi_dmt_xyz = np.dot(np.dot(a, dchi_dmt), np.transpose(a))
    xs = np.array([x, y, z])
    # get the normalization factor
    delta_mt = np.sqrt(2)/np.sqrt(np.sum(dchi_dmt_xyz**2))
    delta_xs = 1/np.sqrt(np.sum(dchi_dxs_xyz**2))
    # * do normalization with step length
    dchi_dxs_xyz = dchi_dxs_xyz*delta_xs*max_dxs_ratio
    dchi_dmt_xyz = dchi_dmt_xyz*delta_mt*max_dxs_ratio
    xs = xs / delta_xs
    mt_xyz = mt_xyz/delta_mt
    # update source
    xs += dchi_dxs_xyz
    mt_xyz += dchi_dmt_xyz
    # denormalize
    xs = xs*dchi_dxs_xyz
    mt_xyz = mt_xyz * delta_xs
    # convert back to rtp coordinate
    x, y, z = xs
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
    lon, lat, alt = pyproj.transform(ecef, lla, x, y, z)

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
