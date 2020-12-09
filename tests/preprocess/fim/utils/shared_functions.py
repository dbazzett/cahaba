from pathlib import Path
import pandas as pd
import rasterio.shutil
import requests
import numpy as np
import rasterio
import geopandas as gpd
from rasterio.warp import calculate_default_transform, reproject, Resampling
import rasterio.crs
from rasterio.merge import merge
from rasterio import features
from shapely.geometry import shape

##############################################################################
#Function Get Metadata
##############################################################################
def get_metadata(ahps_code, walk_upstream, walk_downstream):
    metadata_url = f'http://nwcal-wrds-ti01.nwc.nws.noaa.gov/api/staging/v2/location/metadata/nws_lid/{ahps_code}/?upstream_trace_distance={walk_upstream}&downstream_trace_distance{walk_downstream}'
    response = requests.get(metadata_url)
    if response.ok:
        metadata_json = response.json()
        locations = metadata_json['locations']
        assert len(locations) == 1
        metadata = locations[0]
    return metadata
#######################################################################
#Function to get datum information
#######################################################################
def get_datum(metadata):
    datum = metadata['altitude']
    if not datum == '':
        datum = float(datum)
    vcs = metadata['alt_datum_code']
    lat = float(metadata['latitude'])
    lon = float(metadata['longitude'])
    return vcs, datum, lat, lon
#######################################################################
#Function to get conversion adjustment NGVD to NAVD in FEET
#######################################################################
def ngvd_to_navd_ft(lat, lon):
    #Get conversion from NGVD to NAVD in feet
    #Get vertical adjustment from noaa API
    vertical_datum = f'https://vdatum.noaa.gov/vdatumweb/api/tidal?lon={lon}&lat={lat}&s_v_frame=NGVD29&s_h_frame=NAD27&t_v_frame=NAVD88'
    response = requests.get(vertical_datum)
    if response.ok:      
        metadata_json = response.json()
        vertcon = float(metadata_json['tar_height'])
        datum_adj_ft = round(vertcon*3.281,2)
        return datum_adj_ft
##############################################################################
#Function to get nwm segments
###############################################################################
def get_nwm_segs(metadata):
    upstream_segs = metadata['upstream_nwm_features']
    segment = [metadata['nwm_feature_id']]
    downstream_segs = metadata['downstream_nwm_features']
    all_segments = upstream_segs + segment + downstream_segs
    return all_segments
###############################################################################
#Function to get thresholds
#There appear to be an original and a calculated value; USGS and NWS sources figure out which is best.
def get_threshold(ahps_code):
#Download threshold information (action,minor,moderate and stages and flows)
    threshold_url = f'http://nwcal-wrds-ti01.nwc.nws.noaa.gov/api/staging/v2/location/threshold/all/all/{ahps_code}/'
    response = requests.get(threshold_url)
    if response.ok:
        threshold_json = response.json()
        threshold_data = threshold_json['thresholds'][0]
        thresholds = ['action', 'minor', 'moderate', 'major']
        stages = {}
        flows = {}
        for threshold in thresholds:
            stage = threshold_data['original_values'][f'{threshold}_stage']
            flow =  threshold_data['calculated_values'][f'{threshold}_flow']            
            if not (stage or flow) == 'None':
                stage = float(stage)
                flow = float(flow)
            stages[threshold] = stage
            flows[threshold] = flow
        
        stages['units'] = threshold_data['metadata']['stage_unit']
        flows['units'] = threshold_data['metadata']['flow_unit']
        return stages, flows
########################################################################
#Function to find all ESRI (AIG) or tif grids that can be opened by rasterio in a subdir
########################################################################
def find_grids(directory):
    #Finds all tif files and esri grids that can be opened by rasterio
    grids = [depthgrid for depthgrid in directory.glob('*') if (depthgrid.is_dir() or '.tif' in depthgrid.suffix) and rasterio.shutil.exists(str(depthgrid))]
    grid_name = [i.stem for i in grids]
    return grids, grid_name
########################################################################
#Function to find all shapefiles in a subdir
########################################################################
def find_shp(directory):
    #Finds all tif files and esri grids that can be opened by rasterio
    shapefiles = [shapefile for shapefile in directory.glob('*') if '.shp' in shapefile.suffix]
    shapefile_name = [i.stem for i in shapefiles]
    return shapefiles, shapefile_name
########################################################################
#Specific to AHPS: get polygon CRS
########################################################################
def get_shp_crs(parent_path):
    for subdir in parent_path.glob('*'):
        if 'polygon' in subdir.name.lower() and subdir.is_dir():
            shapefile_paths, shapefile_names = find_shp(subdir)
    if shapefile_paths:
        shapefile = shapefile_paths[0]
        data = gpd.read_file(shapefile)
        shapefile_crs = data.crs
    return shapefile_crs
#######################################################################
#Function to download rating curve from API
#######################################################################
def get_rating_curve(ahps_code):
    #Download rating curve associated with each ahps 
    rating_curve_url = f'http://nwcal-wrds-ti01.nwc.nws.noaa.gov/api/staging/v2/location/rating_curve/{ahps_code}/'
    response = requests.get(rating_curve_url)
    if response.ok:
        site_json = response.json()
        rating_curves = site_json['rating_curves']
        assert len(rating_curves)<=1
        if rating_curves:          
           rating_curve = pd.DataFrame(rating_curves[0]['rating_curve'], dtype=float)
           stage_units = rating_curves[0]['metadata']['stage_unit']
           flow_units = rating_curves[0]['metadata']['flow_unit']
        else:
            rating_curve = pd.DataFrame(columns = ['flow', 'stage'], dtype=float)
        return rating_curve
#######################################################################
#Function to return a correct maps (TEST THIS GOOD). Deal with case where threshold is blank
########################################################################
def select_grids(dataframe, stages:dict):
    thresholds = ['action', 'minor', 'moderate', 'major']
    maps = {}
    for i,threshold in enumerate(thresholds):
        low_bound = stages[threshold]
        if i<=2:
            upper_bound = stages[thresholds[i+1]]
            value = dataframe.query(f'({low_bound}<=stage) & (stage<{upper_bound})')['stage'].min()
            if np.isfinite(value):
                map_path = dataframe.query(f'stage == {value}').path.item()
            else: 
                map_path = 'No Map'
        else:
            value = dataframe.query(f'{low_bound}<=stage')['stage'].min()
            if np.isfinite(value):
                map_path = dataframe.query(f'stage == {value}').path.item()
        maps[threshold] = map_path    
    
    #Get the maximum map to be used as the extent
    max_value = dataframe['stage'].max()
    map_path = dataframe.query(f'stage == {max_value}').path.item()
    maps['extent'] = map_path 

    return maps 
##############################################################################
# Function to write flow file
##############################################################################
def flow_data(segments, flows):
    #Convert cfs to cms
    cfs_to_cms = 0.3048**3
    flows_cms = round(flows * cfs_to_cms,2)
    flow_data = pd.DataFrame({'feature_id':segments, 'discharge':flows_cms})
    flow_data = flow_data.astype({'feature_id' : int , 'discharge' : float})
    return flow_data   
###############################################################################
# Function to preprocess ahps grid
########################################################
def process_grid(benchmark, benchmark_profile, reference):
    #Determine the new transform and dimensions of reprojected/resampled raster.
    new_transform, new_width, new_height = calculate_default_transform(benchmark_profile['crs'].to_wkt(), reference.crs, benchmark.width, benchmark.height, *benchmark.bounds, resolution = reference.res)

    #Define an empty array that is same dimensions as output by the "calculate_default_transform" command. 
    benchmark_projected = np.empty((new_height,new_width), dtype=np.int32)

    #Read in dataset_object
    benchmark_arr = benchmark.read(1)
    
    #Define NODATA value, since using 
    nodata_value = -99

    #Reproject and resample the benchmark dataset. Bilinear resampling due to continuous depth data.
    reproject(benchmark_arr, 
              destination = benchmark_projected,
              src_transform = benchmark.transform, 
              src_crs = benchmark_profile['crs'],
              src_nodata = benchmark.nodata,
              dst_transform = new_transform, 
              dst_crs = reference.crs,
              dst_nodata = nodata_value,
              dst_resolution = reference.res,
              resampling = Resampling.bilinear)

    #Convert entire depth grid to boolean (1 = Flood, 0 = No Flood)
    boolean_nodata_value = 0
    boolean_benchmark = np.where(benchmark_projected != nodata_value, 1, boolean_nodata_value)
    boolean_benchmark = boolean_benchmark.astype('uint8')
    #Update profile (data type, NODATA, transform, width/height).
    profile = reference.profile
    profile.update(transform = new_transform)
    profile.update(dtype = rasterio.uint8)
    profile.update(nodata = boolean_nodata_value) #Update NODATA to some integer so we can keep int8 datatype. There are no NODATA in the raster dataset.
    profile.update (width = new_width)
    profile.update(height = new_height)
    
    return boolean_benchmark, profile
########################################################################
#Mosaic AHPS
########################################################################    
def mosaic_grids(list_of_grids:list):
    #Open each of the grids in rasterio
    grid_objects = [rasterio.open(grid) for grid in list_of_grids]
    nodata_value = 0
    mosaic, out_trans = merge(grid_objects, nodata = nodata_value, method = 'max')
    #Select profile
    profile = grid_objects[0].profile
    #Modify profile
    new_height = mosaic.shape[1]
    new_width = mosaic.shape[2]
    profile.update(transform = out_trans, width = new_width, height = new_height, nodata = nodata_value)
    return mosaic, profile
########################################################################
#Convert raster to polygon
########################################################################
def raster_to_feature(grid_path):
    #Convert a raster to a shapefile
    dataset = rasterio.open(grid_path)
    data = dataset.read(1)
    msk = dataset.read_masks(1)   
    spatial = []
    values = []
    for geom, val in rasterio.features.shapes(data, mask = msk, transform = dataset.transform):
        spatial.append(shape(geom))
        values.append(val)        
    spatial_geodatabase = gpd.GeoDataFrame({'values': values,'geometry':spatial }, crs = dataset.crs)
    dissolve_geodatabase = spatial_geodatabase.dissolve(by = 'values')
    return dissolve_geodatabase

########################################################################
#Mosaic flow files
########################################################################
def merge_flows(list_of_flow_files: list):
    merged_contents = pd.DataFrame()
    for file in list_of_flow_files:
        contents = pd.read_csv(file)
        merged_contents = merged_contents.append(contents, ignore_index = True)        
    return merged_contents
########################################################################
#Aggregate method
########################################################################
def aggregate_method(method:str, ahps_list:list):
    if ahps_list[0] == 'rfc_forecast_points':
        url = 'http://nwcal-wrds-ti01.nwc.nws.noaa.gov/api/staging/v2/location/metadata/nws_lid/all/?must_include=rfc_forecast_point'   
    else:
        sites = ('%2C').join(ahps_list)
        url = f'http://nwcal-wrds-ti01.nwc.nws.noaa.gov/api/staging/v2/location/metadata/nws_lid/({sites})/'    
    response = requests.get(url)
    if response.ok:
        site_json = response.json()
        metadata = site_json['locations']
        df = pd.json_normalize(metadata)       
    if method == 'huc':
        dictionary = df.groupby('huc')['nws_lid'].apply(list).to_dict()   
    if method == 'rfc':
        dictionary = df.groupby('rfc')['nws_lid'].apply(list).to_dict()        
    if method == 'state':
        dictionary = df.groupby('state')['nws_lid'].apply(list).to_dict()
    
    if '' in dictionary.keys():
        dictionary['not_assigned'] = dictionary.pop('')
    return dictionary