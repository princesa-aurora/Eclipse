import numpy as np
from astropy.time import Time
import astropy.units as u
from scipy.spatial import KDTree




def j2000_to_utc_datetime(j2000_seconds) :
    # the input is seconds that have passed since the J2000.0 epoch
    # convert it into a standard UTC datetime string

    j2000_epoch = Time('J2000.0', scale='tt')
    time_tt = j2000_epoch + j2000_seconds * u.s
    time_utc = time_tt.utc

    return time_utc.strftime("%d.%m.%Y %H:%M:%S") + " UTC"



def get_nearest_idcs(source_lon, source_lat, target_lon, target_lat, k) :
    # get the indices in source_lon/lat of the k closest points to the points in target_lon/lat
    # i.e. target_lon/lat has shape(n,) , then return has shape(n,k) (with trivial dimensions being flattened)

    assert np.asarray(source_lon).shape == np.asarray(source_lat).shape
    assert np.asarray(target_lon).shape == np.asarray(target_lat).shape

    # convert spherical to cartesian
    X_source = np.stack((np.cos(source_lat) * np.cos(source_lon),
                                    np.cos(source_lat) * np.sin(source_lon),
                                    np.sin(source_lat)), axis=-1)

    X_target = np.stack((np.cos(target_lat) * np.cos(target_lon),
                    np.cos(target_lat) * np.sin(target_lon),
                    np.sin(target_lat)), axis=-1)

    # use KDTree to find the desired indices
    tree = KDTree(X_source)
    _, nearest_idcs = tree.query(X_target, k=k)

    return nearest_idcs



def get_interp_weights(lon_interp, lat_interp, lon_target, lat_target) :
    # using the 3 interpolation points (lon/lat)
    # get the interpolations weights that are used to linearly interpolate the value at the target point (lon/lat)
    # for this to work the interpolation points of course have to be close to the target point

    assert np.asarray(lon_interp).shape == (3,)
    assert np.asarray(lat_interp).shape == (3,)
    assert np.ndim(lon_target) == 0
    assert np.ndim(lat_target) == 0

    # rotation matrix that rotates X_target to the z_axis
    RotMat = np.array([[np.sin(lat_target)*np.cos(lon_target), np.sin(lat_target)*np.sin(lon_target), -np.cos(lat_target)],
                [-np.sin(lon_target), np.cos(lat_target), 0],
                [np.cos(lat_target)*np.cos(lon_target), np.cos(lat_target)*np.sin(lon_target), np.sin(lat_target)]])


    X_interp = np.column_stack((np.cos(lat_interp) * np.cos(lon_interp),
                                    np.cos(lat_interp) * np.sin(lon_interp),
                                    np.sin(lat_interp)))

    # rotate towards z-axis
    X_interp = X_interp@RotMat.T
    # extract x and y coordinates (z coordinate is very close to 1 after the rotation and carries no important information)
    x_interp = X_interp[:,0]
    y_interp = X_interp[:,1]

    # build interpolation weights from x_interp
    interp_weights = np.array([x_interp[1]*y_interp[2] - y_interp[1]*x_interp[2],
                            x_interp[2]*y_interp[0] - y_interp[2]*x_interp[0],
                            x_interp[0]*y_interp[1] - y_interp[0]*x_interp[1]])
    interp_weights /= interp_weights.sum()

    return interp_weights



def interpolate_linearly(x, x1, f1, x2, f2) :
    # do a standard linear interpolation fro the value at x
    # given the interpolation points (x1,f1) and (x2,f2)

    m = (f2-f1)/(x2-x1)
    t = (f1*x2 - f2*x1)/(x2-x1)

    return m*x+t



