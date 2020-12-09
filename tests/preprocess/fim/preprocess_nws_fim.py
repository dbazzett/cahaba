########################################################
#Preprocess AHPS NWS
########################################################
source_dir = Path(r'/path/to/ahps/datasets')
destination = Path(r'/path/to/outputs')
reference_raster= Path(r'/path/to/reference raster)'

ahps_codes = [i.name for i in source_dir.glob('*') if i.is_dir() and len(i.name) == 5]
all_df = pd.DataFrame()
#Find depth grid subfolder
for code in ahps_codes:
    print(code)
    if code == 'mnda2':
        print('skipping mnda2')
        continue
    #Get Rating Curve for site and calculate elevation from stage
    rating_curve = get_rating_curve(code)
    #If rating curve is empty exit
    if rating_curve.empty:
        print(f'{code} has no rc')
        continue
        #raise RuntimeError(f'Rating Curve is Empty for {code}')
    #Get metadata of site and search for NWM segments 10 miles upstream/10miles downstream (this number may change)
    metadata = get_metadata(code,2,2)   
    #Get the datum to calculate elevation
    vcs, datum, lat, lon = get_datum(metadata)
    if datum == '':
        print(f'{code} has no datum')
        continue
    #Datum adjustment if necessary
    if vcs == 'NGVD29':
        datum = ngvd_to_navd_ft(lat,lon) + datum

    #Add elevation to rating_curve include datum adjustment
    #Assumes that datum supplied by WRDS is in FEET!!!!!!!!!!! Verify
    rating_curve['elevation'] = rating_curve['stage'] + datum
    
    #Search through ahps directory find depth grid folder
    parent_path = source_dir / code    
    for subdir in parent_path.glob('*'):
        if 'depth_grid' in subdir.name.lower() and subdir.is_dir():
            #find all available grids in depth grid directory (.tif or ESRI grids)
            grid_paths, grid_names = find_grids(subdir)       
            #Check if grids are present
            if grid_paths:
                #Dataframe containing grid paths and grids
                df = pd.DataFrame({'code': code, 'path':grid_paths, 'name': grid_names})
                #Determine elevation from the grid name.
                df['elevation'] = df['name'].str.replace('elev_', '', case = False).str.replace('_','.').astype(float)
                #Add datum to dataframe
                df['datum'] = datum                
                #Add a stage column using the datum. Stage is rounded to the nearest 0.1 ft.
                df['stage'] = round(df['elevation'] - datum,1)                                
                #Next steps interpolate a flow for each supplied elevation
                #Ensure elevation is in ascending order
                df.sort_values(by = 'elevation', ascending = True, inplace = True)
                #Interpolate a flow from the elevation
                df['flow'] = np.interp(df['elevation'], rating_curve['elevation'], rating_curve['flow'], left = np.nan, right = np.nan)

                all_df = all_df.append(df)
            else: 
                print(f'{code} has no grids')
    #Get thresholds for action, minor, moderate, major
    stages, flows = get_threshold(code)
    #Select the grid closest to threshold
    grids = select_grids(df, stages)
    #Obtain NWM segments to apply flows
    segments = get_nwm_segs(metadata)
    #Write raster and flow file
    try:
        for i in ['action', 'minor', 'moderate', 'major']:
            flow = flows[i]
            flow_grid = grids[i]
            if flow !='None' and flow_grid != 'No Map':
                flow_info = flow_data(segments,flow) 
                benchmark = rasterio.open(flow_grid)
                benchmark_profile = benchmark.profile
                benchmark_profile.update(driver = 'Gtiff')
                #if grid doesn't have CRS, then assign CRS using a polygon from the ahps inundation library
                if not benchmark.crs:
                    #Populate crs data using the first polygon layer associated with ahps code
                    shapefile_crs = get_shp_crs(parent_path)
                    benchmark_profile.update(crs = shapefile_crs)
                
                #Here we need to transform the benchmark dataset to be like the reference dataset                        
                reference = rasterio.open(reference_raster)
                boolean_benchmark, boolean_profile = process_grid(benchmark, benchmark_profile, reference)
    
                #Output grid and flow file to destination
                outputdir = destination / code / i
                outputdir.mkdir(parents = True, exist_ok = True)
                output_raster = outputdir / (flow_grid.stem + '.tif')
                with rasterio.Env():
                    with rasterio.open(output_raster, 'w', **boolean_profile) as dst:
                        dst.write(boolean_benchmark,1)
                #Write out the flow file to csv
                output_flow_file = outputdir / (flow_grid.stem + '.csv')
                flow_info.to_csv(output_flow_file, index = False)
    except:
        print(f'issue with {code}')                
    #Process extents, only create extent if ahps code has flow file/grids.
    ahps_directory = destination / code
    if ahps_directory.is_dir():               
        extent_grid = grids['extent']
        benchmark = rasterio.open(extent_grid)
        benchmark_profile = benchmark.profile
        benchmark_profile.update(driver = 'Gtiff')

        #if grid doesn't have CRS, then assign CRS using a polygon from the ahps inundation library
        if not benchmark.crs:
            #Populate crs data using the first polygon layer associated with ahps code
            shapefile_crs = get_shp_crs(parent_path)
            benchmark_profile.update(crs = shapefile_crs)

        #Here we need to transform the benchmark dataset to be like the reference dataset        
        boolean_benchmark, boolean_profile = process_grid(benchmark, benchmark_profile, reference)
        
        output_raster = ahps_directory / ('extent.tif')
        with rasterio.Env():
            with rasterio.open(output_raster, 'w', **boolean_profile) as dst:
                dst.write(boolean_benchmark,1)
