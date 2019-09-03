import logging

logging.basicConfig(level=logging.INFO)
log=logging.getLogger('extrapolate_to_grid')

import argparse
import numpy as np
import pandas as pd

from stompy import utils
from stompy.spatial import interp_4d, proj_utils
from stompy.grid import unstructured_grid

##

class ExtrapolateToGrid(object):
    alpha=1e-5 # smoothness in space
    beta=0.5   # smoothness in time
    time_scale=2.0 # time scale for temporal smoothing
    n_layers=1 # increase for 3D output. Assumes all layers dense -- no sparse z-layers.

    dt=np.timedelta64(86400,'s')
    t_ref=None

    plot_mode=None
    
    def __init__(self,grid,data,start,end,output,t_ref,**kw):
        self.grid=grid
        self.data=data
        self.start=start
        self.end=end
        self.output=output
        self.t_ref=t_ref

        # optional arguments
        utils.set_keywords(self,kw)

        # Do the deed
        # DWAQ segment files are written as binary, with each time step as:
        # <4-byte little-endian time in seconds since reference time>
        # <4-byte little-endian floating> * N_segments
        #    data is in segment order, all surface segments, then all segments on the
        #    next layer down, and so on.
        with open(output,'wb') as fp:
            for t in self.time_steps():
                log.info("Processing %s"%( utils.to_datetime(t).strftime('%Y-%m-%d %H:%M')))
                t_sec=(t-self.t_ref)/np.timedelta64(1,'s')
                t_sec=np.array([t_sec],'<i4')
                data2d=self.field_2d(t)
                assert np.all(np.isfinite(data2d)),"Error -- getting some non-finite values"
                data3d=self.extrude_to_3d(data2d)
                fp.write( t_sec.tobytes() )
                fp.write(data3d.astype('<f4'))
                self.post_frame(t,data2d,data3d)

    fig=None
    coll=None
    txt=None
    def post_frame(self,t,data2d,data3d):
        if self.plot_mode is None:
            return

        import matplotlib.pyplot as plt
        import stompy.plot.cmap as scmap
        if self.fig is None:
            self.fig=plt.figure()
            ax=self.fig.add_subplot(1,1,1)
            
            vmin=self.data['value'].min()
            vmax=self.data['value'].max()
            cmap=scmap.load_gradient('turbo.cpt')

            self.coll=self.grid.plot_cells(values=data2d,clim=[vmin,vmax],
                                           cmap=cmap,ax=ax)
            ax.axis('equal')
            plt.colorbar(self.coll,ax=ax)
            self.txt=ax.text(0.05,0.05,"text",transform=ax.transAxes)
        else:
            self.coll.set_array(data2d)
        self.txt.set_text( str(t) )
        self.fig.canvas.draw()
        plt.pause(0.01)
                
    def time_steps(self):
        """
        Return sequence of datetime64 according to start,end and dt.
        """
        return np.arange(self.start, self.end + self.dt, self.dt)
    
    def field_2d(self,t_interp):
        """
        Extrapolate the data to a single point in time given by t_interp.
        Returns an array of the extrapolated values on the 2D grid.
        """
        # interpolate in time
        time_int=interp_4d.interp_to_time_per_ll(self.data,t_interp,
                                                 lat_col='y',lon_col='x')

        # calculate weight that accounts for spacing in time
        time_int['weight']=(1+time_int.time_offset/self.time_scale)**(-self.beta)

        return interp_4d.weighted_grid_extrapolation(self.grid,time_int,alpha=self.alpha)

    def field_3d(self,tstamp):
        return self.extrude_to_3d(self.field_2d(tstamp))

    def extrude_to_3d(self,data):
        """
        data: 1-D array of values, one value per 2D cell (i.e. element)
        returns: 1-D array of values, one value per 3D cell (i.e. segment).
        Assumes dense ordering of segments, and that self.n_layers is set
        correctly.
        """
        return np.tile(data,self.n_layers)

##
# That was called like this from the waq_scenario code:
    def create_temperature_field(self):
        real_func=temp_seg_func.TempSegFunc(self.hydro)
        def seg_func(t_sec): # time conversion done in Scenario since we have the scu
            tstamp=self.time0 + self.scu*t_sec
            return real_func(tstamp)

        # pull times from hydro
        dt=7*86400 # days
        t_secs=np.arange( hydro.t_secs[0]-dt,
                          hydro.t_secs[-1]+2*dt,
                          dt)
        return waq_scenario.ParameterSpatioTemporal(func_t=seg_func,
                                                    times=t_secs)

##


if __name__=="__main__":
    parser=argparse.ArgumentParser(description='Extrapolate point data in space and time.')

    parser.add_argument("-g", "--grid",help="Path to DWAQ grid geometry netcdf.",default=None,required=True)
    
    parser.add_argument("-p", "--projection",help="Map projection for the grid.",default="EPSG:26910")

    # these will default to the period of the data.
    parser.add_argument("-r", "--reference", help="Reference date for DWAQ run (YYYY-MM-DDTHH:MM)", default=None, required=True)
    parser.add_argument("-s", "--start", help="Date of start of output (YYYY-MM-DDTTHH:MM)",default=None)
    parser.add_argument("-e", "--end", help="Date of end of output (YYYY-MM-DDTHH:MM)",default=None)
    
    parser.add_argument("-o", "--output", help="Path to DWAQ segment file for output", default="output.seg")

    parser.add_argument("-d", "--data",help="Input data",nargs='+')

    parser.add_argument("-i", "--interval",help="Time step in output, suffix 's' for seconds, 'D' for days", default='1D')

    parser.add_argument("-l", "--layers", help="Set number of layers for output", default=1)

    parser.add_argument("-m", "--plot", help="Try to display plots of each time step as they are processed",
                        action='store_true')

    parser.add_argument("-a","--alpha", help="Degree of smoothing. Smaller is smoother and more stable.", default=1e-5,
                        type=float)
    
    args=parser.parse_args()

    # LOAD THE GRID
    # for DWAQ flowgeom.nc, does this need to be tweaked?
    if '_net.nc' in args.grid:
        log.info("Reading grid from %s, assuming DFM format"%args.grid)
        grid=unstructured_grid.UnstructuredGrid.read_dfm(args.grid)
    else:
        log.info("Reading grid from %s, assuming UGRID format"%args.grid)
        grid=unstructured_grid.UnstructuredGrid.read_dfm(args.grid)

    # ASSEMBLE INPUT DATA
    dfs=[]
    
    for data_file in args.data:
        if ':' in data_file:
            data_file,field = data_file.split(':')
        else:
            field=None
            
        log.info("Reading data file from %s"%data_file)
        df=pd.read_csv(data_file,parse_dates=['time'])

        if field is None:
            data_cols=[ col
                        for col in df.columns
                        if col not in ['time','longitude','latitude'] ]
            if len(data_cols)!=1:
                raise Exception("Cannot determine data column.  Candidates were: %s"%
                                (' '.join(data_cols)))
            field=data_cols[0]

        # standardize column order
        df=df.loc[:,[ 'time','latitude','longitude',field] ]
        # and standardize field name to something generic
        df=df.rename(columns={field:'value'})
        dfs.append(df)

    combined=pd.concat(dfs)

    # Filter out bad data points
    bad=False
    for field in ['time','latitude','longitude','value']:
        bad=bad|df[field].isnull()
    n_bad=bad.sum()
    if n_bad:
        log.warning("%d records dropped because of missing data"%n_bad)
        df=df[ ~bad ].copy()
    
    dt=np.timedelta64(args.interval[:-1],args.interval[-1])
    
    if args.start is None:
        args.start=utils.floor_dt64(combined.time.min(),dt)
        log.info("Start date defaults to start of data: %s"%args.start)
    else:
        args.start=np.datetime64(args.start)
        
    if args.end is None:
        args.end=utils.ceil_dt64(combined.time.max(),dt)
        log.info("End date defaults to end of data: %s"%args.end)
    else:
        args.end=np.datetime64(args.end)

    args.reference=np.datetime64(args.reference)
        
    # REPROJECT
    log.info("Reprojecting to %s"%args.projection)
    mapper=proj_utils.mapper('WGS84',args.projection)
    ll=df.loc[:,['longitude','latitude']].values
    xy=mapper(ll)
    df['x']=xy[:,0]
    df['y']=xy[:,1]

    if args.plot:
        plot_mode='each'
    else:
        plot_mode=None
        
    ExtrapolateToGrid(grid,data=df,
                      start=args.start,end=args.end,dt=dt,
                      n_layers=args.layers,
                      t_ref=args.reference,
                      plot_mode=plot_mode,
                      alpha=args.alpha,
                      output=args.output)
