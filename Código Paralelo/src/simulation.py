from typing import Tuple
from src.data import Ocean
import numpy as np
import math

EARTH_RADIUS = 6378137
DEG2M = np.deg2rad(1) * EARTH_RADIUS


def unit_vector(lat: np.array, lon: np.array) -> np.array:
    """
    Compute the unit vector of the given lat, lon coordinates.

    Parameters
    ----------
    lat : np.array
        Array containing the latitude of the points.
    lon : np.array
        Array containing the longitude of the points.

    Returns
    -------
    np.array
        Array containing the unit vector of the given lat, lon coordinates.
    """
    lat, lon = np.deg2rad(lat), np.deg2rad(lon)
    x = np.cos(lat) * np.cos(lon)
    y = np.cos(lat) * np.sin(lon)
    z = np.sin(lat)

    return np.array([x, y, z])


def distance_from_points(
    lat_start: np.array,
    lon_start: np.array,
    lat_end: np.array,
    lon_end: np.array,
    radius: float = EARTH_RADIUS,
) -> np.array:
    """
    Compute the distance in meters between two points in a sphere.
    The distance is in meters because `EARTH_RADIUS` is in meters.

    Parameters
    ----------
    lat_start : np.array
        Array containing the latitude of the starting points.
    lon_start : np.array
        Array containing the longitude of the starting points.
    lat_end : np.array
        Array containing the latitude of the ending points.
    lon_end : np.array
        Array containing the longitude of the ending points.
    radius : float
        Radius of the sphere. Default value is the Earth's radius in meters.

    Returns
    -------
    np.array
        Array containing the distances in meters between the two points.
    """
    unit_vector_start = unit_vector(lat_start, lon_start)
    unit_vector_end = unit_vector(lat_end, lon_end)

    # Compute the angle between the two unit vectors
    dot_product = np.sum(unit_vector_start * unit_vector_end, axis=0)
    angle = np.arccos(dot_product)  # radians
    dist = radius * angle

    # Compute the distance between the two points
    return dist


def angle_from_points(
    lat_start: np.array,
    lon_start: np.array,
    lat_end: np.array,
    lon_end: np.array,
) -> np.array:
    """
    Compute the angle in degrees (w.r.t. latitude=azimuth)
    between two points in a sphere.

    Parameters
    ----------
    lat_start : np.array
        Array containing the latitude of the starting points in degrees.
    lon_start : np.array
        Array containing the longitude of the starting points in degrees.
    lat_end : np.array
        Array containing the latitude of the ending points in degrees.
    lon_end : np.array
        Array containing the longitude of the ending points in degrees.

    Returns
    -------
    np.array
        Array containing the angles in degrees between the two points.
    """
    x1, y1, z1 = unit_vector(lat_start, lon_start)
    x2, y2, z2 = unit_vector(lat_end, lon_end)

    return np.rad2deg(
        np.arctan2(-x2 * y1 + x1 * y2, -(x1 * x2 + y1 * y2) * z1 + (x1**2 + y1**2) * z2)
    )


def module_angle_from_components(v: np.array, u: np.array) -> Tuple[np.array]:
    """
    Compute the module and the angle in degrees (with respect to latitude=azimuth)
    of the given vector.

    Parameters
    ----------
    v : np.array
        Array containing the v (latitude) component of the vector.
    u : np.array
        Array containing the u (longitude) component of the vector.

    Returns
    -------
    Tuple[np.array]
        Tuple containing the module and the angle of the vector.
    """

    module = np.sqrt(v**2 + u**2)
    angle = np.rad2deg(np.arctan2(u, v))

    return module, angle


def components_from_module_angle(module: np.array, angle: np.array) -> Tuple[np.array]:
    """
    Compute the v and u components of the given vector.

    Parameters
    ----------
    module : np.array
        Array containing the module of the vector.
    angle : np.array
        Array containing the angle in degrees (with respect to latitude).

    Returns
    -------
    Tuple[np.array]
        Tuple containing the v and u components of the vector.
    """
    angle = np.deg2rad(angle)
    v = module * np.cos(angle)  # When parallel to latitude: 0º -> cos(0º) = 1
    u = module * np.sin(angle)  # When perpendicular to latitude: 90º -> sin(90º) = 1

    return v, u


def compute_time(
    lat_start: np.ndarray,
    lon_start: np.ndarray,
    lat_end: np.ndarray,
    lon_end: np.ndarray,
    vel_ship: float,
    ocean_data: Ocean,
) -> float:
    """
    Function to compute time needed to navigate between two points.

    Parameters
    ----------
    lat_start : np.ndarray
        Latitude of the starting point
    lon_start : np.ndarray
        Longitude of the starting point
    lat_end : np.ndarray
        Latitude of the ending point
    lon_end : np.ndarray
        Longitude of the ending point
    distance : np.ndarray
        Distance between the two points
    vel_ship : float
        Ship velocity in m/s.
    ocean_data : Ocean
        Ocean object.
    clip_p : float
        Percentage to set the minimum of the velocity over ground, by default 0.2

    Returns
    -------
    float
        Seconds per segment
    """
    current_v, current_u = ocean_data.get_currents(lat_start, lon_start)
    if (math.isnan(current_v) or math.isnan(current_u)):
        return float('nan')
    
    current_v_f, current_u_f = ocean_data.get_currents(lat_end, lon_end)
    if (math.isnan(current_v_f) or math.isnan(current_u_f)):
        return float('nan')

    cur_mod, cur_ang = module_angle_from_components(current_v, current_u)

    angle = angle_from_points(lat_start, lon_start, lat_end, lon_end)
    distance = distance_from_points(lat_start, lon_start, lat_end, lon_end)

    # Angle of the current projected on the line of movement over ground
    cur_ang_proj = cur_ang - angle

    # u is the components perpendicular to the ship direction over ground
    # v is the components parallel to the ship direction over ground
    cur_u_proj, cur_v_proj = components_from_module_angle(cur_mod, cur_ang_proj)

    # Clip the perpendicular velocity to avoid negative values
    # We assume the vessel will increase its power to compensate
    vel_v2 = np.clip(vel_ship**2 - cur_v_proj**2, a_min=0, a_max=None)
    vel_ground = cur_u_proj + np.sqrt(vel_v2)

    return distance / vel_ground
