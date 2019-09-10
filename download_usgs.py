import logging as log
import os
import numpy as np

from stompy.io.local import usgs_sfbay
from stompy import utils

cache_dir='cache'

def download_usgs(start,end,output,z_min=None, z_max=None, field='Temperature'):
    if not os.path.exists(cache_dir):
        log.info("Creating local cache directory: %s"%cache_dir)
        os.makedirs(cache_dir)
        
    ds=usgs_sfbay.usgs_sfbay_dataset(start,end,
                                     cache_dir=cache_dir)
    ds=ds.copy()
    
    depth_sel=np.ones(ds.dims['depth'],np.bool8)
    if z_min is not None:
        depth_sel = depth_sel & (-ds.depth.values>=z_min)
    if z_max is not None:
        depth_sel = depth_sel & (-ds.depth.values<=z_max)
    ds_zavg=ds.isel(depth=depth_sel).mean(dim='depth')
    df=ds_zavg[field].to_dataframe()

    # limit to just the columns we care about.
    df=df.reset_index().loc[:, ['time','latitude','longitude',field]]

    # Drop missing
    valid=(~df['time'].isnull())&(~df[field].isnull())
    df_valid=df[ valid ]
    df_valid.to_csv(output,index=False)
    
    return df_valid

import argparse
if __name__=="__main__":
    parser=argparse.ArgumentParser(description='Download USGS transect data, format for extrapolation.')

    # these will default to the period of the data.
    parser.add_argument("-s", "--start", help="Starting date", required=True)
    parser.add_argument("-e", "--end", help="Ending date", required=True)
    parser.add_argument("-o", "--output", help="Path to CSV file for output", default=None)
    args=parser.parse_args()

    start=np.datetime64(args.start)
    end  =np.datetime64(args.end)

    if args.output is None:
        args.output="usgs_sfbay-%s-%s.csv"%( utils.to_datetime(start).strftime('%Y%m%d'),
                                             utils.to_datetime(end).strftime('%Y%m%d') )
        log.info("Output file: %s"%args.output)

    download_usgs(start=start,end=end,output=args.output)


# For direct testing     
#   download_usgs(start=np.datetime64("2012-09-01"),
#                 end= np.datetime64("2013-11-01"),
#                 output="test/test-usgs.csv")


