import glob
import os
import rasterio
import numpy as np
import geopandas as gpd
import shapely
from shapely import Point
import srtm4
from rpcm import RPCModel, rpc_from_rpc_file
import geopy
import geopy.distance


def localize_row_col_geometry(
    geometry: shapely.geometry.Point,
    rpc_model: RPCModel,
):
    """
    Project points in the row and column space to the geographic space, with altitude correction.
    All the altitudes here are in meters.
    Args:
        geometry: The geometry to add altitudes to, should be a shapely point.
        rpc_model: the RPCModel to base the projection on.
    returns:
        A tuple containing (new longitude, new latitude).
    """
    points = []
    col_idx, row_idx = geometry.coords[0]
    points.append((col_idx, row_idx))
    col_indices, row_indices = np.asarray(points).T
    coords = [
        rpc_model.localization(col_idx, row_idx, 0)
        for col_idx, row_idx in zip(col_indices, row_indices)
    ]
    corrected_lons, corrected_lats = [], []
    for (lon, lat), col_idx, row_idx in zip(coords, col_indices, row_indices):
        
        new_altitude = srtm4.srtm4(lon, lat)
        new_lon, new_lat = 0, 0
        max_iterations = 10
        ground_dist_tol = 50
        ground_distance = ground_dist_tol + 1
        iteration_no = 0
        while (ground_distance > ground_dist_tol) and (
            iteration_no < max_iterations
        ):
            iteration_no += 1
            old_lon, old_lat = new_lon, new_lat
            new_lon, new_lat = rpc_model.localization(
                col_idx, row_idx, new_altitude
            )
            new_altitude = srtm4.srtm4(new_lon, new_lat)
            ground_distance = geopy.distance.distance(
                (old_lat, old_lon), (new_lat, new_lon)
            ).m
        corrected_lons.append(new_lon)
        corrected_lats.append(new_lat)
    points.append(
        ((lon, lat) for (lon, lat) in zip(corrected_lons, corrected_lats))
    )
            
    corrected_lons.append(new_lon)
    corrected_lats.append(new_lat)
    points.append(
        ((lon, lat) for (lon, lat) in zip(corrected_lons, corrected_lats))
    )
    return (new_lon, new_lat)


def read_and_prepare_matches(path_to_all_matches: str, match_columns: list=[0, 1]):
    """
    Loads in the merged matches file and returns a geodataframe containing Point-like geometry.
    
    path_to_all_matches: The path to the text file containing all of the matches.
    match_columns: The columns to be used as x, y for the matches geodataframe output.
    """
    matches = np.loadtxt(path_to_all_matches)
    matches = matches.astype(int)
    matches = matches[:, match_columns]
    matches_x = matches[:, 0]
    matches_y = matches[:, 1]
    matches_gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy(matches_x, matches_y))
    return matches_gdf


def correct_points(matches_gdf: gpd.GeoDataFrame, rpc_model_path: str, every_x_match: int=10):
    """
    matches_gdf: A geodataframe containing point geometries of the matches created by s2p in image coordinates.
    rpc_model_path: A path to the rpc model to be used to calculate real coordinates from image coordinates.
    every_x_match: For decreasing the size of the matches geodataframe, takes every point, i.e 5 is every 5.
    
    return: two numpy arrays, (lons and lats).
    """
    corrected_points = []
    if type(rpc_model_path) == str:
        rpc_model = rpc_from_rpc_file(rpc_model_path)
    else:
        rpc_model = rpc_model_path
    matches_gdf = matches_gdf["geometry"][::every_x_match]
    
    # Iterate over matches and correct them.
    for i, point in enumerate(matches_gdf):
        localized_point = localize_row_col_geometry(point, rpc_model)
        corrected_points.append(list(localized_point))
        
    # Prepare the lons and lats for returning
    corrected_points = np.array(corrected_points)
    lons = corrected_points[:, 0]
    lats = corrected_points[:, 1]
    
    return lons, lats


def points_to_geojson(out_path: str, lons: np.ndarray, lats: np.ndarray):
    """
    Saves the lons and lats from correct_points as a GeoJSON.
    
    out_path: The path to save the output GeoJSON to,
    lons: The longitudes returned by correct_points.
    lats: The latitudes returned by correct_points.
    """
    gdf_corrected = gpd.GeoDataFrame(geometry=gpd.points_from_xy(lons, lats, crs="EPSG:4326"))
    gdf_corrected.to_file(out_path, driver="GeoJSON")
    
    
def matches_to_geojson(path_to_matches: str, path_to_rpc: str, every_x_match: int, match_columns: list, out_path: str):
    matches_gdf = read_and_prepare_matches(path_to_matches, match_columns)
    lons, lats = correct_points(matches_gdf, path_to_rpc, every_x_match)
    points_to_geojson(out_path, lons, lats)
    

    