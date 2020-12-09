########################################################################
#This section will mine through the destination directory and get all the grids using the aggregation
destination_dir = Path(r'/path/to/aggregated/output')
ahps_dir = Path(r'/path/to/prepocessed fim')
ahps_codes = [i.name for i in ahps_dir.glob('*') if i.is_dir() and len(i.name) == 5]
thresholds = ['action','minor','moderate','major']
dictionary = aggregate_method('huc', ahps_codes)
for aggregate in dictionary:   
    ahps_sites = dictionary[aggregate]      

    #Loop through each threshold and make a composite grid and flow file. Place the output grid and flow file in a directory based on the aggregation method.
    for threshold in thresholds:       
        all_directories = [(ahps_dir / site / threshold) for site in ahps_sites]
        grid_paths = [i for directory in all_directories for i in find_grids(directory)[0]]
        flowfile_paths = [i for directory in all_directories for i in directory.glob('*.csv')]
        if grid_paths:
            try:
                merged_flowfile = merge_flows(flowfile_paths)
                mosaic, mosaic_profile = mosaic_grids(grid_paths)
                
                #Output grid and flow file to destination
                outputdir = destination_dir / aggregate / threshold
                outputdir.mkdir(parents = True, exist_ok = True)
                output_raster = outputdir / ('ahps_huc_' + aggregate + "_depth_" + threshold + '.tif')
                
                with rasterio.Env():
                    with rasterio.open(output_raster, 'w', **mosaic_profile) as dst:
                        dst.write(mosaic)            
                
                output_flowfile = outputdir / ("ahps_huc_" + aggregate + "_flows_" + threshold + '.csv')
                merged_flowfile.to_csv(output_flowfile, index = False)                            
        
            except:
                print(f'error with merge for {aggregate} and {threshold}')
    
    #Loop through 
    aggregate_dir = destination_dir / aggregate
    if aggregate_dir.is_dir():
        #create extent grid
        extent_grids = [(ahps_dir / site / 'extent.tif') for site in ahps_sites]
        #Create a merged extent layer with attribute for each ahps code.
        all_extents = gpd.GeoDataFrame()
        for extent in extent_grids:
            extent_geodataframe = raster_to_feature(extent)
            site = extent.parent.stem
            extent_geodataframe['code'] = site
            all_extents = all_extents.append(extent_geodataframe)
        
        extent_filepath = aggregate_dir / ("ahps_huc_" + aggregate + "_extent.shp")
        all_extents.to_file(extent_filepath)

#if needed merge all extent shapefiles together
extent_shapefiles = list(destination_dir.rglob('*_extent.shp'))
all_extents = gpd.GeoDataFrame()
for extent in extent_shapefiles:
    layer = gpd.read_file(extent)
    all_extents = all_extents.append(layer)
all_extents.to_file(destination_dir/'ahps_extents.shp')