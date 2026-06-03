import numpy as np
from astropy.time import Time
import astropy.units as u
from scipy.spatial import KDTree
from global_land_mask import globe
import reverse_geocoder as rg
import pycountry
from PIL import Image



BOLD = '\033[1m'
RESET = '\033[0m'



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


def is_on_land(lon, lat) :
    # return True is the point is on land, False is over water
    assert np.asarray(lon).shape == np.asarray(lat).shape

    lon_deg = lon*180/np.pi
    lat_deg = lat*180/np.pi

    if (np.ndim(lat) == 0):
        if lon_deg > 180: lon_deg -= 360
    else:
        lon_deg[lon_deg > 180] -= 360

    return globe.is_land(lat_deg, lon_deg)



def get_location_name(lon, lat):
    # get a name for a location given by lon/lat
    lon_deg = lon *180/np.pi
    lat_deg = lat *180/np.pi
    if lon_deg > 180: lon_deg -= 360

    # check if the coordinate is over water
    if not is_on_land(lon, lat):
        # figure out what ocean it is in (only roughly, smaller seas like the mediterranean are omitted)
        if lat_deg > 60.0:
            return "Arctic Ocean"
        elif lat_deg < -60.0:
            return "Southern Ocean"
        elif lon_deg < -67.3:
            return "Pacific Ocean"
        elif lon_deg < 20.0:
            return "Atlantic Ocean"
        elif lon_deg < 146.9:
            return "Indian Ocean"
        else:
            return "Pacific Ocean"
    
    # check for Antarctica (since it is uninhabited it wouldn't be recognized by reverse geocoder)
    if lat_deg < -60:
        return "Antarctica"

    # if on land and not in Antarctica, safely run the offline reverse geocoder
    coordinates = (lat_deg, lon_deg)
    results = rg.search(coordinates, verbose=False)

    if results:
        info = results[0]
        name = info.get('name')
        admin1 = info.get('admin1')
        cc = info.get('cc')

        try:
            # convert from country name abbreviation to full name
            country_lookup = pycountry.countries.get(alpha_2=cc)
            country_name = getattr(country_lookup, 'common_name', country_lookup.name)
        except (KeyError, AttributeError):
            country_name = cc 

        if admin1 and admin1 != name:
            return f"{name} ({admin1}), {country_name}"
        else:
            return f"{name}, {country_name}"
        
    return "Unknown Land Location"



def mean_of_lons_and_lats(lons, lats):
    # compute the mean of the positions in cartesian forms
    # and then convert back to spherical

    x = np.stack((np.cos(lats)*np.cos(lons), np.cos(lats)*np.sin(lons), np.sin(lats)), axis=-1)
    x_mean = np.mean(x, axis=0)

    lon_mean = np.arctan2(x_mean[1], x_mean[0])
    lat_mean = np.arctan(x_mean[2]/np.sqrt(x_mean[0]**2+x_mean[1]**2))

    return lon_mean, lat_mean



def is_valid_png(file_path):
    # is the file a valid png?
    try:
        with Image.open(file_path) as img:
            img.verify() # PIL does all the work of checking the file
        return True
    except (FileNotFoundError, IOError, SyntaxError):
        return False



