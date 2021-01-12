#!/bin/bash -e

## INITIALIZE TOTAL TIME TIMER ##
T_total_start

echo -e $startDiv"Parameter Values"
echo -e "extent=$extent"
echo -e "negativeBurnValue=$negativeBurnValue"
echo -e "maxSplitDistance_meters=$maxSplitDistance_meters"
echo -e "mannings_n=$manning_n"
echo -e "stage_min_meters=$stage_min_meters"
echo -e "stage_interval_meters=$stage_interval_meters"
echo -e "stage_max_meters=$stage_max_meters"
echo -e "slope_min=$slope_min"
echo -e "ms_buffer_dist=$ms_buffer_dist"
echo -e "ncores_gw=$ncores_gw"
echo -e "ncores_fd=$ncores_fd"
echo -e "defaultMaxJobs=$defaultMaxJobs"
echo -e "memfree=$memfree"$stopDiv

## SET OUTPUT DIRECTORY FOR UNIT ##
# Define and create output directory for each Hydrologic Unit Code (HUC)
hucNumber="$1"
outputHucDataDir=$outputRunDataDir/$hucNumber
mkdir $outputHucDataDir


## SET VARIABLES AND FILE INPUTS ##
# Define variables for HUC ids and some input files
hucUnitLength=${#hucNumber}
huc4Identifier=${hucNumber:0:4}
huc2Identifier=${hucNumber:0:2}
input_NHD_WBHD_layer=WBDHU$hucUnitLength
input_DEM=$inputDataDir/nhdplus_rasters/HRNHDPlusRasters"$huc4Identifier"/elev_cm.tif
input_NLD=$inputDataDir/nld_vectors/huc2_levee_lines/nld_preprocessed_"$huc2Identifier".gpkg

# Define the landsea water body mask using either Great Lakes or Ocean polygon input #
if [[ $huc2Identifier == "04" ]] ; then
  input_LANDSEA=$inputDataDir/landsea/gl_water_polygons.gpkg
  echo -e "Using $input_LANDSEA for water body mask (Great Lakes)"
else
  input_LANDSEA=$inputDataDir/landsea/water_polygons_us.gpkg
fi


## GET WBD ##
# Get the NHD Watershed Boundary Dataset (WBD) for input HUC(s)
echo -e $startDiv"Get WBD $hucNumber"$stopDiv
date -u
Tstart
[ ! -f $outputHucDataDir/wbd.gpkg ] && \
ogr2ogr -f GPKG $outputHucDataDir/wbd.gpkg $input_WBD_gdb $input_NHD_WBHD_layer -where "HUC$hucUnitLength='$hucNumber'"
Tcount


## BUFFER WBD ##
# Buffer WBD by 5000m to ensure all data cover the input HUC(s)
echo -e $startDiv"Buffer WBD $hucNumber"$stopDiv
date -u
Tstart
[ ! -f $outputHucDataDir/wbd_buffered.gpkg ] && \
ogr2ogr -f GPKG -dialect sqlite -sql "select ST_buffer(geom, 5000) from 'WBDHU$hucUnitLength'" $outputHucDataDir/wbd_buffered.gpkg $outputHucDataDir/wbd.gpkg
Tcount


## GET STREAMS ##
# Get Vector Layers (NHD flowlines, NWM lakes, NWM flowlines, NLD lines, LandSea polygon) and subset/clip/mask to the wbd_buffered boundary for the input HUC.
echo -e $startDiv"Get Vector Layers and Subset $hucNumber"$stopDiv
date -u
Tstart
[ ! -f $outputHucDataDir/NHDPlusBurnLineEvent_subset.gpkg ] && \
$libDir/snap_and_clip_to_nhd.py -d $hucNumber -w $input_NWM_Flows -f $input_NWM_Headwaters -s $input_NHD_Flowlines -l $input_NWM_Lakes -r $input_NLD -u $outputHucDataDir/wbd.gpkg -g $outputHucDataDir/wbd_buffered.gpkg -y $inputDataDir/ahp_sites/nws_lid.gpkg -v $input_LANDSEA -c $outputHucDataDir/NHDPlusBurnLineEvent_subset.gpkg -z $outputHucDataDir/nld_subset_levees.gpkg -a $outputHucDataDir/nwm_lakes_proj_subset.gpkg -t $outputHucDataDir/nwm_headwaters_proj_subset.gpkg -m $input_NWM_Catchments -n $outputHucDataDir/nwm_catchments_proj_subset.gpkg -e $outputHucDataDir/nhd_headwater_points_subset.gpkg -x $outputHucDataDir/LandSea_subset.gpkg -b $outputHucDataDir/nwm_subset_streams.gpkg -p $extent
Tcount

if [ "$extent" = "MS" ]; then
  if [[ ! -f $outputHucDataDir/NHDPlusBurnLineEvent_subset.gpkg ]] ; then
    echo "No AHPs point(s) within HUC $hucNumber boundaries. Aborting run_by_unit.sh"
    rm -rf $outputHucDataDir
    exit 0
  fi
fi


## Clip WBD8 ##
# Clip WBD for all HUC 8s within the input HUC list
echo -e $startDiv"Clip WBD8"$stopDiv
date -u
Tstart
ogr2ogr -f GPKG -clipsrc $outputHucDataDir/wbd_buffered.gpkg $outputHucDataDir/wbd8_clp.gpkg $inputDataDir/wbd/WBD_National.gpkg WBDHU8
Tcount


## CLIP DEM ##
# Clip HUC DEM to the buffered WBD(s)
echo -e $startDiv"Clip DEM $hucNumber"$stopDiv
date -u
Tstart
[ ! -f $outputHucDataDir/dem.tif ] && \
gdalwarp -cutline $outputHucDataDir/wbd_buffered.gpkg -crop_to_cutline -ot Int32 -r bilinear -of "GTiff" -overwrite -co "BLOCKXSIZE=512" -co "BLOCKYSIZE=512" -co "TILED=YES" -co "COMPRESS=LZW" -co "BIGTIFF=YES" $input_DEM $outputHucDataDir/dem.tif
Tcount


## GET RASTER METADATA
# Get DEM Metadata (fsize, ncols, nrows, ndv, xmin, ymin, xmax, ymax, cellsize_resx, cellsize_resy)
echo -e $startDiv"Get DEM Metadata $hucNumber"$stopDiv
date -u
Tstart
read fsize ncols nrows ndv xmin ymin xmax ymax cellsize_resx cellsize_resy<<<$($libDir/getRasterInfoNative.py $outputHucDataDir/dem.tif)


## RASTERIZE NLD MULTILINES ##
# Rasterize all NLD polylines using zelev vertices (interpolating raster values btw vertices). Elev units are in feet (NAVD88).
echo -e $startDiv"Rasterize all NLD multilines using zelev vertices"$stopDiv
date -u
Tstart
[ ! -f $outputHucDataDir/nld_rasterized_elev.tif ] && [ -f $outputHucDataDir/nld_subset_levees.gpkg ] && \
gdal_rasterize -l nld_subset_levees -3d -at -init $ndv -te $xmin $ymin $xmax $ymax -ts $ncols $nrows -ot Float32 -of GTiff -co "COMPRESS=LZW" -co "BIGTIFF=YES" -co "TILED=YES" $outputHucDataDir/nld_subset_levees.gpkg $outputHucDataDir/nld_rasterized_elev.tif
Tcount


## CONVERT TO METERS ##
# Convert DEM units from centimeters to meters → pixel value/100. Also converts type from Int32 to Float32
echo -e $startDiv"Convert DEM to Meters $hucNumber"$stopDiv
date -u
Tstart
[ ! -f $outputHucDataDir/dem_meters.tif ] && \
gdal_calc.py --quiet --type=Float32 --co "BLOCKXSIZE=512" --co "BLOCKYSIZE=512" --co "TILED=YES" --co "COMPRESS=LZW" --co "BIGTIFF=YES" -A $outputHucDataDir/dem.tif --outfile="$outputHucDataDir/dem_meters.tif" --calc="A/100" --NoDataValue=$ndv
Tcount


## RASTERIZE REACH BOOLEAN (1 & 0) ##
# Rasterize reach with boolean output - 1’s for streamline pixels and 0’s for non-streamline pixels
echo -e $startDiv"Rasterize Reach Boolean $hucNumber"$stopDiv
date -u
Tstart
[ ! -f $outputHucDataDir/flows_grid_boolean.tif ] && \
gdal_rasterize -ot Int32 -burn 1 -init 0 -co "COMPRESS=LZW" -co "BIGTIFF=YES" -co "TILED=YES" -te $xmin $ymin $xmax $ymax -ts $ncols $nrows $outputHucDataDir/NHDPlusBurnLineEvent_subset.gpkg $outputHucDataDir/flows_grid_boolean.tif
Tcount


## RASTERIZE NHD HEADWATERS (1 & 0) ##
# Rasterize NHD Headwaters
echo -e $startDiv"Rasterize NHD Headwaters $hucNumber"$stopDiv
date -u
Tstart
[ ! -f $outputHucDataDir/headwaters.tif ] && \
gdal_rasterize -ot Int32 -burn 1 -init 0 -co "COMPRESS=LZW" -co "BIGTIFF=YES" -co "TILED=YES" -te $xmin $ymin $xmax $ymax -ts $ncols $nrows $outputHucDataDir/nhd_headwater_points_subset.gpkg $outputHucDataDir/headwaters.tif
Tcount


## RASTERIZE NWM CATCHMENTS ##
# Rasterize NWM Catchments polygons
if [ "$extent" = "FR" ]; then
  echo -e $startDiv"Raster NWM Catchments $hucNumber"$stopDiv
  date -u
  Tstart
  [ ! -f $outputHucDataDir/nwm_catchments_proj_subset.tif ] && \
  gdal_rasterize -ot Int32 -a ID -a_nodata 0 -init 0 -co "COMPRESS=LZW" -co "BIGTIFF=YES" -co "TILED=YES" -te $xmin $ymin $xmax $ymax -ts $ncols $nrows $outputHucDataDir/nwm_catchments_proj_subset.gpkg $outputHucDataDir/nwm_catchments_proj_subset.tif
  Tcount
fi


## BURN LEVEES INTO DEM ##
# Only runs if levees exist in HUC. Burn/merge NLD levee elevation into dem. Replaces DEM values with NLD levee elevation values wherever NLD zelev is greater than original DEM. Also converts NLD elev values from feet to meters.
echo -e $startDiv"Burn nld levees into dem & convert nld elev to meters (*Overwrite dem_meters.tif output) $hucNumber"$stopDiv
date -u
Tstart
[ -f $outputHucDataDir/nld_rasterized_elev.tif ] && \
gdal_calc.py --quiet --type=Float32 --overwrite --NoDataValue $ndv --co "BLOCKXSIZE=512" --co "BLOCKYSIZE=512" --co "TILED=YES" --co "COMPRESS=LZW" --co "BIGTIFF=YES" -A $outputHucDataDir/dem_meters.tif -B $outputHucDataDir/nld_rasterized_elev.tif --outfile="$outputHucDataDir/dem_meters.tif" --calc="maximum(A,(B*0.3048))" --NoDataValue=$ndv
Tcount


## DEM Reconditioning ##
# Using AGREE methodology, hydroenforce the DEM so that it is consistent
# with the supplied stream network. This allows for more realistic catchment
# delineation which is ultimately reflected in the output FIM mapping.
echo -e $startDiv"Creating AGREE DEM using $buffer meter buffer"$stopDiv
date -u
Tstart
[ ! -f $outputHucDataDir/dem_burned.tif ] && \
$libDir/agreedem.py -r $outputHucDataDir/flows_grid_boolean.tif -d $outputHucDataDir/dem_meters.tif -w $outputHucDataDir -g $outputHucDataDir/temp_work -o $outputHucDataDir/dem_burned.tif -b $buffer -sm 10 -sh 1000
Tcount


## PIT REMOVE BURNED DEM ##
# Remove sinks using RichDEM Depression Filling technique.
echo -e $startDiv"Pit remove Burned DEM $hucNumber"$stopDiv
date -u
Tstart
[ ! -f $outputHucDataDir/dem_burned_filled.tif ] && \
rd_depression_filling $outputHucDataDir/dem_burned.tif $outputHucDataDir/dem_burned_filled.tif
Tcount


## D8 FLOW DIR ##
# D8 Flow Directions on Burned DEM. Creates a grid of flow direction from each grid cell to one of its adjacent or diagonal neighbors, calculated using the direction of steepest descent.
echo -e $startDiv"D8 Flow Directions on Burned DEM $hucNumber"$stopDiv
date -u
Tstart
[ ! -f $outputHucDataDir/flowdir_d8_burned_filled.tif ] && \
mpiexec -n $ncores_fd $taudemDir2/d8flowdir -fel $outputHucDataDir/dem_burned_filled.tif -p $outputHucDataDir/flowdir_d8_burned_filled.tif
Tcount


## DINF FLOW DIR ##
# echo -e $startDiv"DINF on Filled Thalweg Conditioned DEM"$stopDiv
# date -u
# Tstart
# [ ! -f $outputHucDataDir/flowdir_dinf_thalwegCond.tif] && \
# mpiexec -n $ncores_fd $taudemDir2/dinfflowdir -fel $outputHucDataDir/dem_thalwegCond_filled.tif -ang $outputHucDataDir/flowdir_dinf_thalwegCond.tif -slp $outputHucDataDir/slopes_dinf.tif
# Tcount


## D8 FLOW ACCUMULATIONS ##
# Calculates a grid of contributing areas using the single direction D8 flow model. 
echo -e $startDiv"D8 Flow Accumulations $hucNumber"$stopDiv
date -u
Tstart
$taudemDir/aread8 -p $outputHucDataDir/flowdir_d8_burned_filled.tif -ad8  $outputHucDataDir/flowaccum_d8_burned_filled.tif -wg  $outputHucDataDir/headwaters.tif -nc
Tcount


# THRESHOLD ACCUMULATIONS ##
# Outputs an indicator (1,0) grid identifying cells with input values >= the threshold value. The standard use is to use an accumulated source area grid to as the input grid to generate a stream raster grid as the output.
echo -e $startDiv"Threshold Accumulations $hucNumber"$stopDiv
date -u
Tstart
$taudemDir/threshold -ssa $outputHucDataDir/flowaccum_d8_burned_filled.tif -src  $outputHucDataDir/demDerived_streamPixels.tif -thresh 1
Tcount


## PREPROCESSING FOR LATERAL THALWEG ADJUSTMENT ###
# Assign a unique ID for each stream pixel and writes to file. It then uses this raster to run GRASS r.grow.distance tool to create the allocation and proximity rasters required to complete the lateral thalweg conditioning.
echo -e $startDiv"Preprocessing for lateral thalweg adjustment $hucNumber"$stopDiv
date -u
Tstart
$libDir/unique_pixel_and_allocation.py -s $outputHucDataDir/demDerived_streamPixels.tif -o $outputHucDataDir/demDerived_streamPixels_ids.tif -g $outputHucDataDir/temp_grass
Tcount


## ADJUST THALWEG MINIMUM USING LATERAL ZONAL MINIMUM ##
# Performing lateral thalweg adjustment. Using 50m search radius. Algorithm searches for the zonal minimum elevation in each pixel catchment. It updates the catchment_min_dict with this zonal minimum elevation value.
echo -e $startDiv"Performing lateral thalweg adjustment $hucNumber"$stopDiv
date -u
Tstart
$libDir/adjust_thalweg_lateral.py -e $outputHucDataDir/dem_meters.tif -s $outputHucDataDir/demDerived_streamPixels.tif -a $outputHucDataDir/demDerived_streamPixels_ids_allo.tif -d $outputHucDataDir/demDerived_streamPixels_ids_dist.tif -t 50 -o $outputHucDataDir/dem_lateral_thalweg_adj.tif
Tcount


## MASK BURNED DEM FOR STREAMS ONLY ###
# Mask Burned DEM for Thalweg Only
echo -e $startDiv"Mask Burned DEM for Thalweg Only $hucNumber"$stopDiv
date -u
Tstart
[ ! -f $outputHucDataDir/flowdir_d8_burned_filled_flows.tif ] && \
gdal_calc.py --quiet --type=Int32 --overwrite --co "COMPRESS=LZW" --co "BIGTIFF=YES" --co "TILED=YES" -A $outputHucDataDir/flowdir_d8_burned_filled.tif -B $outputHucDataDir/demDerived_streamPixels.tif --calc="A/B" --outfile="$outputHucDataDir/flowdir_d8_burned_filled_flows.tif" --NoDataValue=0
Tcount


## FLOW CONDITION STREAMS ##
# Flow Condition Thalweg. Soft burns in elevations by tracking down D8 flow directions to ensure there is no uphill elevation.
echo -e $startDiv"Flow Condition Thalweg $hucNumber"$stopDiv
date -u
Tstart
[ ! -f $outputHucDataDir/dem_thalwegCond.tif ] && \
$taudemDir/flowdircond -p $outputHucDataDir/flowdir_d8_burned_filled_flows.tif -z $outputHucDataDir/dem_lateral_thalweg_adj.tif -zfdc $outputHucDataDir/dem_thalwegCond.tif
Tcount


## D8 SLOPES ##
# D8 Slopes from DEM. Contain the slope, as evaluated in the direction of steepest descent, and is reported as drop/distance, i.e. tan of the angle
echo -e $startDiv"D8 Slopes from DEM $hucNumber"$stopDiv
date -u
Tstart
mpiexec -n $ncores_fd $taudemDir2/d8flowdir -fel $outputHucDataDir/dem_lateral_thalweg_adj.tif -sd8 $outputHucDataDir/slopes_d8_dem_meters.tif
Tcount


## STREAMNET FOR REACHES ##
# This tool orders the stream network according to the Strahler ordering system. When two reaches of equal order join the downstream reach order is increased by 1.
echo -e $startDiv"Stream Net for Reaches $hucNumber"$stopDiv
date -u
Tstart
[ ! -f $outputHucDataDir/demDerived_reaches.shp ] && \
$taudemDir/streamnet -p $outputHucDataDir/flowdir_d8_burned_filled.tif -fel $outputHucDataDir/dem_thalwegCond.tif -ad8 $outputHucDataDir/flowaccum_d8_burned_filled.tif -src $outputHucDataDir/demDerived_streamPixels.tif -ord $outputHucDataDir/streamOrder.tif -tree $outputHucDataDir/treeFile.txt -coord $outputHucDataDir/coordFile.txt -w $outputHucDataDir/sn_catchments_reaches.tif -net $outputHucDataDir/demDerived_reaches.shp
Tcount


## SPLIT DERIVED REACHES ##
# 1) split stream segments based on lake boundaries and input threshold distance
# 2) calculate channel slope, manning's n, and LengthKm for each segment
# 3) create unique ids using HUC8 boundaries (and unique 'fossid' column)
# 4) create network traversal attribute columns (To_Node, From_Node, NextDownID)
# 5) create points layer with segment verticies encoded with HydroID's (used for catchment delineation in next step)
echo -e $startDiv"Split Derived Reaches $hucNumber"$stopDiv
date -u
Tstart
[ ! -f $outputHucDataDir/demDerived_reaches_split.gpkg ] && \
$libDir/split_flows.py $outputHucDataDir/demDerived_reaches.shp $outputHucDataDir/dem_thalwegCond.tif $outputHucDataDir/demDerived_reaches_split.gpkg $outputHucDataDir/demDerived_reaches_split_points.gpkg $maxSplitDistance_meters $slope_min $outputHucDataDir/wbd8_clp.gpkg $outputHucDataDir/nwm_lakes_proj_subset.gpkg $lakes_buffer_dist_meters
Tcount

if [[ ! -f $outputHucDataDir/demDerived_reaches_split.gpkg ]] ; then
  echo "No AHPs point(s) within HUC $hucNumber boundaries. Aborting run_by_unit.sh"
  rm -rf $outputHucDataDir
  exit 0
fi


## MASK RASTERS BY MS BUFFER ##
# Mainstem Only: Mask Rasters with Stream Buffer
if [ "$extent" = "MS" ]; then
  echo -e $startDiv"Mask Rasters with Stream Buffer $hucNumber"$stopDiv
  date -u
  Tstart
  $libDir/fr_to_ms_raster_mask.py $outputHucDataDir/demDerived_reaches_split.gpkg $outputHucDataDir/flowdir_d8_burned_filled.tif $outputHucDataDir/dem_thalwegCond.tif $outputHucDataDir/slopes_d8_dem_meters.tif $outputHucDataDir/flowdir_d8_MS.tif $outputHucDataDir/dem_thalwegCond_MS.tif $outputHucDataDir/slopes_d8_dem_metersMS.tif $outputHucDataDir/demDerived_streamPixels.tif $outputHucDataDir/demDerived_streamPixelsMS.tif $ms_buffer_dist
  Tcount

  if [[ ! -f $outputHucDataDir/dem_thalwegCond_MS.tif ]] ; then
    echo "No AHPs point(s) within HUC $hucNumber boundaries. Aborting run_by_unit.sh"
    rm -rf $outputHucDataDir
    exit 0
  fi

  dem_thalwegCond=$outputHucDataDir/dem_thalwegCond_MS.tif
  slopes_d8_dem_meters=$outputHucDataDir/slopes_d8_dem_metersMS.tif
  flowdir_d8_burned_filled=$outputHucDataDir/flowdir_d8_MS.tif
  demDerived_streamPixels=$outputHucDataDir/demDerived_streamPixelsMS.tif
else
  dem_thalwegCond=$outputHucDataDir/dem_thalwegCond.tif
  slopes_d8_dem_meters=$outputHucDataDir/slopes_d8_dem_meters.tif
  flowdir_d8_burned_filled=$outputHucDataDir/flowdir_d8_burned_filled.tif
  demDerived_streamPixels=$outputHucDataDir/demDerived_streamPixels.tif
fi


## GAGE WATERSHED FOR REACHES ##
# Calculates Gage watersheds grid. Each grid cell is labeled with the identifier (from column id) of the gage to which it drains directly without passing through any other gages.
echo -e $startDiv"Gage Watershed for Reaches $hucNumber"$stopDiv
date -u
Tstart
[ ! -f $outputHucDataDir/gw_catchments_reaches.tif ] && \
mpiexec -n $ncores_gw $taudemDir/gagewatershed -p $flowdir_d8_burned_filled -gw $outputHucDataDir/gw_catchments_reaches.tif -o $outputHucDataDir/demDerived_reaches_split_points.gpkg -id $outputHucDataDir/idFile.txt
Tcount


## VECTORIZE FEATURE ID CENTROIDS ##
# Vectorize Pixel Centroids
echo -e $startDiv"Vectorize Pixel Centroids $hucNumber"$stopDiv
date -u
Tstart
[ ! -f $outputHucDataDir/flows_points_pixels.gpkg ] && \
$libDir/reachID_grid_to_vector_points.py $demDerived_streamPixels $outputHucDataDir/flows_points_pixels.gpkg featureID
Tcount


## GAGE WATERSHED FOR PIXELS ##
# Calculates Gage watersheds grid. Each grid cell is labeled with the identifier (from column id) of the gage to which it drains directly without passing through any other gages.
echo -e $startDiv"Gage Watershed for Pixels $hucNumber"$stopDiv
date -u
Tstart
[ ! -f $outputHucDataDir/gw_catchments_pixels.tif ] && \
mpiexec -n $ncores_gw $taudemDir/gagewatershed -p $flowdir_d8_burned_filled -gw $outputHucDataDir/gw_catchments_pixels.tif -o $outputHucDataDir/flows_points_pixels.gpkg -id $outputHucDataDir/idFile.txt
Tcount


# D8 REM ##
# Calculates REM/HAND/Detrended DEM
echo -e $startDiv"D8 REM $hucNumber"$stopDiv
date -u
Tstart
[ ! -f $outputHucDataDir/rem.tif ] && \
$libDir/rem.py -d $dem_thalwegCond -w $outputHucDataDir/gw_catchments_pixels.tif -o $outputHucDataDir/rem.tif -t $demDerived_streamPixels
Tcount


## DINF DISTANCE DOWN ##
# echo -e $startDiv"DINF Distance Down on Filled Thalweg Conditioned DEM $hucNumber"$stopDiv
# date -u
# Tstart
# [ ! -f $outputHucDataDir/flowdir_dinf_thalwegCond.tif] && \
# mpiexec -n $ncores_fd $taudemDir/dinfdistdown -ang $outputHucDataDir/flowdir_dinf_thalwegCond.tif -fel $outputHucDataDir/dem_thalwegCond_filled.tif -src $demDerived_streamPixels -dd $outputHucDataDir/rem.tif -m ave h
# Tcount


## BRING DISTANCE DOWN TO ZERO ##
# Zero out negative values in distance down grid
echo -e $startDiv"Zero out negative values in distance down grid $hucNumber"$stopDiv
date -u
Tstart
[ ! -f $outputHucDataDir/rem_zeroed.tif ] && \
gdal_calc.py --quiet --type=Float32 --overwrite --co "COMPRESS=LZW" --co "BIGTIFF=YES" --co "TILED=YES" -A $outputHucDataDir/rem.tif --calc="(A*(A>=0))" --NoDataValue=$ndv --outfile=$outputHucDataDir/"rem_zeroed.tif"
Tcount


## POLYGONIZE REACH WATERSHEDS ##
# Creates vector polygons for all connected regions of pixels in the raster sharing a common pixel value.
echo -e $startDiv"Polygonize Reach Watersheds $hucNumber"$stopDiv
date -u
Tstart
[ ! -f $outputHucDataDir/gw_catchments_reaches.gpkg ] && \
gdal_polygonize.py -8 -f GPKG $outputHucDataDir/gw_catchments_reaches.tif $outputHucDataDir/gw_catchments_reaches.gpkg catchments HydroID
Tcount


## PROCESS CATCHMENTS AND MODEL STREAMS STEP 1 ##
# Filters the reaches and catchments to only include features associated with current huc (using “fosid” of wbdhu8).
echo -e $startDiv"Process catchments and model streams step 1 $hucNumber"$stopDiv
date -u
Tstart
[ ! -f $outputHucDataDir/gw_catchments_reaches_filtered_addedAttributes.gpkg ] && \
$libDir/filter_catchments_and_add_attributes.py $outputHucDataDir/gw_catchments_reaches.gpkg $outputHucDataDir/demDerived_reaches_split.gpkg $outputHucDataDir/gw_catchments_reaches_filtered_addedAttributes.gpkg $outputHucDataDir/demDerived_reaches_split_filtered.gpkg $outputHucDataDir/wbd8_clp.gpkg $hucNumber

if [[ ! -f $outputHucDataDir/gw_catchments_reaches_filtered_addedAttributes.gpkg ]] ; then
  echo "No relevant streams within HUC $hucNumber boundaries. Aborting run_by_unit.sh"
  rm -rf $outputHucDataDir
  exit 0
fi
Tcount


## GET RASTER METADATA ## *****
# Get Clipped Raster Metadata
echo -e $startDiv"Get Clipped Raster Metadata $hucNumber"$stopDiv
date -u
Tstart
read fsize ncols nrows ndv_clipped xmin ymin xmax ymax cellsize_resx cellsize_resy<<<$($libDir/getRasterInfoNative.py $outputHucDataDir/gw_catchments_reaches.tif)
Tcount


## RASTERIZE NEW CATCHMENTS AGAIN ##
# Rasterize filtered catchments
echo -e $startDiv"Rasterize filtered catchments $hucNumber"$stopDiv
date -u
Tstart
[ ! -f $outputHucDataDir/gw_catchments_reaches_filtered_addedAttributes.tif ] && \
gdal_rasterize -ot Int32 -a HydroID -a_nodata 0 -init 0 -co "COMPRESS=LZW" -co "BIGTIFF=YES" -co "TILED=YES" -te $xmin $ymin $xmax $ymax -ts $ncols $nrows $outputHucDataDir/gw_catchments_reaches_filtered_addedAttributes.gpkg $outputHucDataDir/gw_catchments_reaches_filtered_addedAttributes.tif
Tcount


## RASTERIZE LANDSEA (OCEAN AREA) POLYGON (IF APPLICABLE) ##
# Rasterize filtered/dissolved ocean/glake polygon. Assigning ndv to all water body pixels (Ocean & Great Lakes) and a value of 1 everywhere else.
echo -e $startDiv"Rasterize filtered/dissolved ocean/Glake polygon $hucNumber"$stopDiv
date -u
Tstart
[ -f $outputHucDataDir/LandSea_subset.gpkg ] && [ ! -f $outputHucDataDir/LandSea_subset.tif ] && \
gdal_rasterize -ot Int32 -burn $ndv -a_nodata $ndv -init 1 -co "COMPRESS=LZW" -co "BIGTIFF=YES" -co "TILED=YES" -te $xmin $ymin $xmax $ymax -ts $ncols $nrows $outputHucDataDir/LandSea_subset.gpkg $outputHucDataDir/LandSea_subset.tif
Tcount


## MASK SLOPE RASTER ##
# Masking Slope Raster to HUC
echo -e $startDiv"Masking Slope Raster to HUC $hucNumber"$stopDiv
date -u
Tstart
[ ! -f $outputHucDataDir/slopes_d8_dem_meters_masked.tif ] && \
gdal_calc.py --quiet --type=Float32 --overwrite --co "COMPRESS=LZW" --co "BIGTIFF=YES" --co "TILED=YES" -A $slopes_d8_dem_meters -B $outputHucDataDir/gw_catchments_reaches_filtered_addedAttributes.tif --calc="(A*(B>0))+((B<=0)*-1)" --NoDataValue=-1 --outfile=$outputHucDataDir/"slopes_d8_dem_meters_masked.tif"
Tcount


## MASK REM RASTER ##
# Masking REM Raster to HUC
echo -e $startDiv"Masking REM Raster to HUC $hucNumber"$stopDiv
date -u
Tstart
[ ! -f $outputHucDataDir/rem_zeroed_masked.tif ] && \
gdal_calc.py --quiet --type=Float32 --overwrite --co "COMPRESS=LZW" --co "BIGTIFF=YES" --co "TILED=YES" -A $outputHucDataDir/rem_zeroed.tif -B $outputHucDataDir/gw_catchments_reaches_filtered_addedAttributes.tif --calc="(A*(B>0))" --NoDataValue=$ndv --outfile=$outputHucDataDir/"rem_zeroed_masked.tif"
Tcount


## MASK REM RASTER TO REMOVE OCEAN AREAS ##
# Additional masking to REM raster to remove ocean/glake areas in HUC
echo -e $startDiv"Additional masking to REM raster to remove ocean/Glake areas in HUC $hucNumber"$stopDiv
date -u
Tstart
[ -f $outputHucDataDir/LandSea_subset.tif ] && \
gdal_calc.py --quiet --type=Float32 --overwrite --co "COMPRESS=LZW" --co "BIGTIFF=YES" --co "TILED=YES" -A $outputHucDataDir/rem_zeroed_masked.tif -B $outputHucDataDir/LandSea_subset.tif --calc="(A*B)" --NoDataValue=$ndv --outfile=$outputHucDataDir/"rem_zeroed_masked.tif"
Tcount


## MAKE CATCHMENT AND STAGE FILES ##
# Generate Catchment List and Stage List Files. Stage values in meters.
echo -e $startDiv"Generate Catchment List and Stage List Files $hucNumber"$stopDiv
date -u
Tstart
$libDir/make_stages_and_catchlist.py $outputHucDataDir/demDerived_reaches_split_filtered.gpkg $outputHucDataDir/gw_catchments_reaches_filtered_addedAttributes.gpkg $outputHucDataDir/stage.txt $outputHucDataDir/catchment_list.txt $stage_min_meters $stage_interval_meters $stage_max_meters
Tcount


## HYDRAULIC PROPERTIES ##
# Hydraulic Properties
echo -e $startDiv"Hydraulic Properties $hucNumber"$stopDiv
date -u
Tstart
[ ! -f $outputHucDataDir/src_base.csv ] && \
$taudemDir/catchhydrogeo -hand $outputHucDataDir/rem_zeroed_masked.tif -catch $outputHucDataDir/gw_catchments_reaches_filtered_addedAttributes.tif -catchlist $outputHucDataDir/catchment_list.txt -slp $outputHucDataDir/slopes_d8_dem_meters_masked.tif -h $outputHucDataDir/stage.txt -table $outputHucDataDir/src_base.csv
Tcount


## FINALIZE CATCHMENTS AND MODEL STREAMS ##
# Finalize catchments and model streams. Generate a library csv to “crosswalk” (match) hydroids to NWM ids. FR: crosswalk using majority catchment method. MS: crosswalk using stream segment midpoint method.
echo -e $startDiv"Finalize catchments and model streams $hucNumber"$stopDiv
date -u
Tstart
[ ! -f $outputHucDataDir/gw_catchments_reaches_filtered_addedAttributes_crosswalked.gpkg ] && \
$libDir/add_crosswalk.py -d $outputHucDataDir/gw_catchments_reaches_filtered_addedAttributes.gpkg -a $outputHucDataDir/demDerived_reaches_split_filtered.gpkg -s $outputHucDataDir/src_base.csv -l $outputHucDataDir/gw_catchments_reaches_filtered_addedAttributes_crosswalked.gpkg -f $outputHucDataDir/demDerived_reaches_split_filtered_addedAttributes_crosswalked.gpkg -r $outputHucDataDir/src_full_crosswalked.csv -j $outputHucDataDir/src.json -x $outputHucDataDir/crosswalk_table.csv -t $outputHucDataDir/hydroTable.csv -w $outputHucDataDir/wbd8_clp.gpkg -b $outputHucDataDir/nwm_subset_streams.gpkg -y $outputHucDataDir/nwm_catchments_proj_subset.tif -m $manning_n -z $input_NWM_Catchments -p $extent
Tcount


## CLEANUP OUTPUTS ##
# Cleaning up outputs for different run types based on input from user
echo -e $startDiv"Cleaning up outputs $hucNumber"$stopDiv
args=()
[[ ! -z "$whitelist" ]] && args+=( "-w$whitelist" )
(( production == 1 )) && args+=( '-p' )
(( viz == 1 )) && args+=( '-v' )
date -u
Tstart
$libDir/output_cleanup.py $hucNumber $outputHucDataDir "${args[@]}"
Tcount
