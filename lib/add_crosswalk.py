#!/usr/bin/env python3

import geopandas as gpd
import pandas as pd
from numpy import unique
from rasterstats import zonal_stats
import json
import argparse
import sys

def add_crosswalk(input_catchments_fileName,input_flows_fileName,input_srcbase_fileName,output_catchments_fileName,output_flows_fileName,output_src_fileName,output_src_json_fileName,output_crosswalk_fileName,output_hydro_table_fileName,input_huc_fileName,input_nwmflows_fileName,input_nwmcatras_fileName,mannings_n,input_nwmcat_fileName,extent,calibration_mode=False):

    input_catchments = gpd.read_file(input_catchments_fileName)
    input_flows = gpd.read_file(input_flows_fileName)
    input_huc = gpd.read_file(input_huc_fileName)
    input_nwmflows = gpd.read_file(input_nwmflows_fileName)

    if extent == 'FR':
        ## crosswalk using majority catchment method

        # calculate majority catchemnts
        majority_calc = zonal_stats(input_catchments, input_nwmcatras_fileName, stats=['majority'], geojson_out=True)
        input_majorities = gpd.GeoDataFrame.from_features(majority_calc)
        input_majorities = input_majorities.rename(columns={'majority' : 'feature_id'})

        input_majorities = input_majorities[:][input_majorities['feature_id'].notna()]
        if input_majorities.feature_id.dtype != 'int': input_majorities.feature_id = input_majorities.feature_id.astype(int)
        if input_majorities.HydroID.dtype != 'int': input_majorities.HydroID = input_majorities.HydroID.astype(int)

        input_nwmflows = input_nwmflows.rename(columns={'ID':'feature_id'})
        if input_nwmflows.feature_id.dtype != 'int': input_nwmflows.feature_id = input_nwmflows.feature_id.astype(int)
        relevant_input_nwmflows = input_nwmflows[input_nwmflows['feature_id'].isin(input_majorities['feature_id'])]
        relevant_input_nwmflows = relevant_input_nwmflows.filter(items=['feature_id','order_'])

        if input_catchments.HydroID.dtype != 'int': input_catchments.HydroID = input_catchments.HydroID.astype(int)
        output_catchments = input_catchments.merge(input_majorities[['HydroID','feature_id']],on='HydroID')
        output_catchments = output_catchments.merge(relevant_input_nwmflows[['order_','feature_id']],on='feature_id')

        if input_flows.HydroID.dtype != 'int': input_flows.HydroID = input_flows.HydroID.astype(int)
        output_flows = input_flows.merge(input_majorities[['HydroID','feature_id']],on='HydroID')
        if output_flows.HydroID.dtype != 'int': output_flows.HydroID = output_flows.HydroID.astype(int)
        output_flows = output_flows.merge(relevant_input_nwmflows[['order_','feature_id']],on='feature_id')
        output_flows = output_flows.merge(output_catchments.filter(items=['HydroID','areasqkm']),on='HydroID')

    elif extent == 'MS':
        ## crosswalk using stream segment midpoint method
        input_nwmcat = gpd.read_file(input_nwmcat_fileName, mask=input_huc)
        input_nwmcat = input_nwmcat.rename(columns={'ID':'feature_id'})
        if input_nwmcat.feature_id.dtype != 'int': input_nwmcat.feature_id = input_nwmcat.feature_id.astype(int)
        input_nwmcat=input_nwmcat.set_index('feature_id')

        input_nwmflows = input_nwmflows.rename(columns={'ID':'feature_id'})
        if input_nwmflows.feature_id.dtype != 'int': input_nwmflows.feature_id = input_nwmflows.feature_id.astype(int)

        # Get stream midpoint
        stream_midpoint = []
        hydroID = []
        for i,lineString in enumerate(input_flows.geometry):
            hydroID = hydroID + [input_flows.loc[i,'HydroID']]
            stream_midpoint = stream_midpoint + [lineString.interpolate(0.05,normalized=True)]

        input_flows_midpoint = gpd.GeoDataFrame({'HydroID':hydroID, 'geometry':stream_midpoint}, crs=input_flows.crs, geometry='geometry')
        input_flows_midpoint = input_flows_midpoint.set_index('HydroID')

        # Create crosswalk
        crosswalk = gpd.sjoin(input_flows_midpoint, input_nwmcat, how='left', op='within').reset_index()
        crosswalk = crosswalk.rename(columns={"index_right": "feature_id"})
        crosswalk = crosswalk.filter(items=['HydroID', 'feature_id'])
        crosswalk = crosswalk.merge(input_nwmflows[['feature_id','order_']],on='feature_id')

        if input_catchments.HydroID.dtype != 'int': input_catchments.HydroID = input_catchments.HydroID.astype(int)
        output_catchments = input_catchments.merge(crosswalk,on='HydroID')

        if input_flows.HydroID.dtype != 'int': input_flows.HydroID = input_flows.HydroID.astype(int)
        output_flows = input_flows.merge(crosswalk,on='HydroID')
        output_flows = output_flows.merge(output_catchments.filter(items=['HydroID','areasqkm']),on='HydroID')

    # read in manning's n values
    if calibration_mode == False:
        with open(mannings_n, "r") as read_file:
            mannings_dict = json.load(read_file)
    else:
        mannings_dict = {}
        for cnt,value in enumerate(mannings_n.split(",")[2:]):
            streamorder = cnt+1
            mannings_dict[str(streamorder)] = value

    output_flows['ManningN'] = output_flows['order_'].astype(str).map(mannings_dict)

    # Adjust short model reach rating curves
print("Adjusting model reach rating curve")
maxArea   = 0.25
maxLength = 0.5
where = "{} < {} and {} < {}".format("areasqkm", maxArea, "LengthKm", maxLength)
sml_segs = pd.DataFrame()

    # if segment is too short, find upstream segment and apply it's channel geometry
for tab_ind in output_flows.index:
    if output_flows["areasqkm"][tab_ind] < 0.25 and output_flows["LengthKm"][tab_ind] < 0.5:
      short_id = output_flows['HydroID'][tab_ind]
      if (len(output_flows.loc[output_flows['NextDownID'] == short_id]['HydroID']) > 1):
          # max_da = max(output_flows.loc[output_flows['NextDownID'] == short_id][FN_PUOrder])
          # up_id = output_flows.loc[(output_flows[FN_PUOrder] == max_da) & (output_flows['NextDownID'] == short_id)]['HydroID'].item()
          print ("multiple up ids")
      else:
          up_id = output_flows.loc[output_flows['NextDownID'] == short_id]['HydroID'].item()
      sml_segs = sml_segs.append({'short_id':short_id, 'up_id':up_id}, ignore_index=True)
    print("Number of short reaches [{}] = {}".format(where, len(sml_segs)))

    # calculate src_full
    input_src_base = pd.read_csv(input_srcbase_fileName, dtype= object)
    if input_src_base.CatchId.dtype != 'int': input_src_base.CatchId = input_src_base.CatchId.astype(int)

    input_src_base = input_src_base.merge(output_flows[['ManningN','HydroID']],left_on='CatchId',right_on='HydroID')

    input_src_base = input_src_base.rename(columns=lambda x: x.strip(" "))
    input_src_base = input_src_base.apply(pd.to_numeric,**{'errors' : 'coerce'})
    input_src_base['TopWidth (m)'] = input_src_base['SurfaceArea (m2)']/input_src_base['LENGTHKM']/1000
    input_src_base['WettedPerimeter (m)'] = input_src_base['BedArea (m2)']/input_src_base['LENGTHKM']/1000
    input_src_base['WetArea (m2)'] = input_src_base['Volume (m3)']/input_src_base['LENGTHKM']/1000
    input_src_base['HydraulicRadius (m)'] = input_src_base['WetArea (m2)']/input_src_base['WettedPerimeter (m)']
    input_src_base['HydraulicRadius (m)'].fillna(0, inplace=True)
    input_src_base['Discharge (m3s-1)'] = input_src_base['WetArea (m2)']* \
    pow(input_src_base['HydraulicRadius (m)'],2.0/3)* \
    pow(input_src_base['SLOPE'],0.5)/input_src_base['ManningN']

    # set nans to 0
    input_src_base.loc[input_src_base['Stage']==0,['Discharge (m3s-1)']] = 0

    output_src = input_src_base.drop(columns=['CatchId'])
    if output_src.HydroID.dtype != 'int': output_src.HydroID = output_src.HydroID.astype(int)

    # update rating curves
    print("Update rating curves for short reaches.")
    for index, row in sml_segs.iterrows():
        if pd.notnull(row[1]):
          new_values = input_src_base.loc[input_src_base['HydroID'] == row[1]][['Stage', 'Discharge (m3s-1)']]
          for pf_index, pf_row in new_values.iterrows():
              input_src_base.loc[(input_src_base['HydroID']== row[0]) & (input_src_base['Stage']== pf_row[0]),['Discharge (m3s-1)']] = pf_row[1]
        else:
          print("Reach HydroID: {} does not have upstream reach - update skipped".format(row[1]))

    if extent == 'FR':
        output_src = output_src.merge(input_majorities[['HydroID','feature_id']],on='HydroID')
    elif extent == 'MS':
        output_src = output_src.merge(crosswalk[['HydroID','feature_id']],on='HydroID')

    output_crosswalk = output_src[['HydroID','feature_id']]
    output_crosswalk = output_crosswalk.drop_duplicates(ignore_index=True)

    # make hydroTable
    output_hydro_table = output_src.loc[:,['HydroID','feature_id','Stage','Discharge (m3s-1)']]
    output_hydro_table.rename(columns={'Stage' : 'stage','Discharge (m3s-1)':'discharge_cms'},inplace=True)
    if output_hydro_table.HydroID.dtype != 'str': output_hydro_table.HydroID = output_hydro_table.HydroID.astype(str)
    output_hydro_table['HydroID'] = output_hydro_table.HydroID.str.zfill(8)
    output_hydro_table['fossid'] = output_hydro_table.loc[:,'HydroID'].apply(lambda x : str(x)[0:4])
    if input_huc.fossid.dtype != 'str': input_huc.fossid = input_huc.fossid.astype(str)
    output_hydro_table = output_hydro_table.merge(input_huc.loc[:,['fossid','HUC8']],how='left',on='fossid')
    # if output_hydro_table.HydroID.dtype != 'int': output_hydro_table.HydroID = output_hydro_table.HydroID.astype(int)
    if output_flows.HydroID.dtype != 'str': output_flows.HydroID = output_flows.HydroID.astype(str)
    output_flows['HydroID'] = output_flows.HydroID.str.zfill(8)
    output_hydro_table = output_hydro_table.merge(output_flows.loc[:,['HydroID','LakeID']],how='left',on='HydroID')
    output_hydro_table['LakeID'] = output_hydro_table['LakeID'].astype(int)
    output_hydro_table = output_hydro_table.rename(columns={'HUC8':'HUC'})
    output_hydro_table.drop(columns='fossid',inplace=True)
    if output_hydro_table.feature_id.dtype != 'str': output_hydro_table.feature_id = output_hydro_table.feature_id.astype(str)

    # write out based on mode
    if calibration_mode == True:
        output_hydro_table.to_csv(output_hydro_table_fileName,index=False)
    else:
        # make src json
        output_src_json = dict()
        hydroID_list = unique(output_src['HydroID'])

        for hid in hydroID_list:
            indices_of_hid = output_src['HydroID'] == hid
            stage_list = output_src['Stage'][indices_of_hid].astype(float)
            q_list = output_src['Discharge (m3s-1)'][indices_of_hid].astype(float)

            stage_list = stage_list.tolist()
            q_list = q_list.tolist()

            output_src_json[str(hid)] = { 'q_list' : q_list , 'stage_list' : stage_list }

        # write out
        output_catchments.to_file(output_catchments_fileName, driver="GPKG",index=False)
        output_flows.to_file(output_flows_fileName, driver="GPKG", index=False)
        output_src.to_csv(output_src_fileName,index=False)
        output_crosswalk.to_csv(output_crosswalk_fileName,index=False)
        output_hydro_table.to_csv(output_hydro_table_fileName,index=False)

        with open(output_src_json_fileName,'w') as f:
            json.dump(output_src_json,f,sort_keys=True,indent=2)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Subset vector layers')
    parser.add_argument('-d','--input-catchments-fileName', help='DEM derived catchments', required=True)
    parser.add_argument('-a','--input-flows-fileName', help='DEM derived streams', required=True)
    parser.add_argument('-s','--input-srcbase-fileName', help='Base synthetic rating curve table', required=True)
    parser.add_argument('-l','--output-catchments-fileName', help='Subset crosswalked catchments', required=True)
    parser.add_argument('-f','--output-flows-fileName', help='Subset crosswalked  streams', required=True)
    parser.add_argument('-r','--output-src-fileName', help='Output crosswalked synthetic rating curve table', required=True)
    parser.add_argument('-j','--output-src-json-fileName',help='Output synthetic rating curve json',required=True)
    parser.add_argument('-x','--output-crosswalk-fileName',help='Crosswalk table',required=True)
    parser.add_argument('-t','--output-hydro-table-fileName',help='Hydrotable',required=True)
    parser.add_argument('-w','--input-huc-fileName',help='HUC8 boundary',required=True)
    parser.add_argument('-b','--input-nwmflows-fileName',help='Subest NWM burnlines',required=True)
    parser.add_argument('-y','--input-nwmcatras-fileName',help='NWM catchment raster',required=False)
    parser.add_argument('-m','--mannings-n',help='Mannings n. Accepts single parameter set or list of parameter set in calibration mode. Currently input as csv.',required=True)
    parser.add_argument('-z','--input-nwmcat-fileName',help='NWM catchment polygon',required=True)
    parser.add_argument('-p','--extent',help='MS or FR extent',required=True)
    parser.add_argument('-c','--calibration-mode',help='Mannings calibration flag',required=False,action='store_true')

    args = vars(parser.parse_args())

    input_catchments_fileName = args['input_catchments_fileName']
    input_flows_fileName = args['input_flows_fileName']
    input_srcbase_fileName = args['input_srcbase_fileName']
    output_catchments_fileName = args['output_catchments_fileName']
    output_flows_fileName = args['output_flows_fileName']
    output_src_fileName = args['output_src_fileName']
    output_src_json_fileName = args['output_src_json_fileName']
    output_crosswalk_fileName = args['output_crosswalk_fileName']
    output_hydro_table_fileName = args['output_hydro_table_fileName']
    input_huc_fileName = args['input_huc_fileName']
    input_nwmflows_fileName = args['input_nwmflows_fileName']
    input_nwmcatras_fileName = args['input_nwmcatras_fileName']
    mannings_n = args['mannings_n']
    input_nwmcat_fileName = args['input_nwmcat_fileName']
    extent = args['extent']
    calibration_mode = args['calibration_mode']

    add_crosswalk(input_catchments_fileName,input_flows_fileName,input_srcbase_fileName,output_catchments_fileName,output_flows_fileName,output_src_fileName,output_src_json_fileName,output_crosswalk_fileName,output_hydro_table_fileName,input_huc_fileName,input_nwmflows_fileName,input_nwmcatras_fileName,mannings_n,input_nwmcat_fileName,extent,calibration_mode)
