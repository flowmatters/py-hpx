import pandas as pd
import numpy as np
import xarray as xr
import sys
import os
from glob import glob
import re
from datetime import datetime

HP1_DIR=os.path.join(os.getcwd(),'HP1')
HP1_EXE='hp1.exe'
OUTPUTS={
    'Volume':{
        'time':True,
        'file':'T_LEVEL',
        'column':'Volume'
        
    }
}
DEFAULT_SIM_START=datetime(2000,1,1)
DEFAULT_SIM_END=datetime(2000,1,10)

M_TO_MM = 1e-3

ATMOS_HEADER='''Pcp_File_Version=4
*** BLOCK I: ATMOSPHERIC INFORMATION  **********************************
   MaxAL                    (MaxAL = number of atmospheric data-records)
%s
 DailyVar  SinusVar  lLay  lBCCycles lInterc lDummy  lDummy  lDummy  lDummy  lDummy
       f       f       f       f       f       f       f       f       f       f
 hCritS                 (max. allowed pressure head at the soil surface)
      0
'''
ATMOS_FOOTER="end*** END OF INPUT FILE 'ATMOSPH.IN' **********************************"

DEFAULT_CSV_OPTIONS={
    'sep':r'\s+'
}

CSV_OPTIONS={
    'A_LEVEL':{
        'header':0,
        'skiprows':[0,1,3,4],
        'skipfooter':1
    },
    '^solute':{
        'header':0,
        'skiprows':[0,1,3],
        'skipfooter':1
    },
    'BALANCE':{
        'skiprows':[0,1,2,3,4,5,6,7,8,10],
        'skipfooter':3
    },
    'NOD_INF':{
        'skiprows':[0,1,2,3,4,5,6,7,8,9,11,12],
        'skipfooter':1
    },
    'PROFILE':{
        'skiprows':[0,1,2,3,4,5,6,8],
        'skipfooter':1
    },
    'T_LEVEL':{
        'skiprows':[0,1,2,3,4,5,7,8],
        'skipfooter':1
    },
    'RUN_INF':{
        'skiprows':[0,1,2,3,4,5,6,8],
        'skipfooter':1
    },
    'OBS_NODE':{
        'skiprows':10,
        'skipfooter':1
    }
}

def write_atmosph_in(df,fn=None):
    WIDTH=12
    HEADER=('' +''.join([col.rjust(WIDTH) for col in df.columns]) + '\n')[1:]
    float_format=lambda x: ('' if np.isnan(x) else ('%f'%x)).rjust(WIDTH-1)
    format_first = {
        'tAtm':lambda x: ('' if np.isnan(x) else ('%f'%x)).rjust(WIDTH)
    }
    to_str_options = {
    #    'formatters':formats,
        'formatters':format_first,
        'float_format':float_format,
        'col_space':WIDTH-1,
        'header':False,
        'index':False,
        'na_rep':''
    }
    PREFIX=0
    result = (ATMOS_HEADER%str(len(df)).rjust(7)) + HEADER + (' '*PREFIX) + df.to_string(**to_str_options) + '\n' + ATMOS_FOOTER
    if fn is None:
        return result

    if os.path.isdir(fn):
        fn = os.path.join(fn,'ATMOSPH.IN')

    f = open(fn,'w')
    f.write(result)
    f.close()

def read_atmosph_in(fn):
    if os.path.isdir(fn):
        fn = os.path.join(fn,'ATMOSPH.IN')
    return pd.read_table(fn,skiprows=8,delimiter=r'\s+',skipfooter=1,engine='python')

def run_hydrus(sim):
    cwd = os.getcwd()
    try:
        sim_abs = os.path.abspath(sim)
        os.chdir(sim)
        assert os.path.exists(HP1_EXE)

        #print(sim_abs)
        f = open('path.dat','w')
        f.write(sim_abs)
        f.close()
        cmd_line = f'{HP1_EXE} {os.getcwd()}'
        print(f'Running with command: {cmd_line}')
        res = os.system(cmd_line)
        assert res==0
    finally:
        os.chdir(cwd)

    return read_outputs(sim)

def read_output_file(fn):
    base = os.path.basename(fn)
    options = DEFAULT_CSV_OPTIONS.copy()
    for k,v in CSV_OPTIONS.items():
        if re.match(k,base) is not None:
            options.update(**v)
            break
    
    if 'skipfooter' in options:
        options['engine']='python'

    return pd.read_table(fn,**options)
    
def read_outputs(model):
    model = os.path.abspath(model)
    output_filenames = list(glob(os.path.join(model,'*.out')))
    result = {}
    for fn in output_filenames:
        base = os.path.basename(fn)
        try:
            result[base[:-4]] = read_output_file(fn)
        except:
            pass
            #print('Error reading %s'%base)


    return result

def extract_silo(var,time_period,lat,lon,silo_path):
    ds = xr.open_mfdataset(glob(os.path.join(silo_path,'*.%s.nc'%var)))
    res = ds[var].loc[time_period,(lat-1e-6):(lat+1e-6),(lon-1e-6):(lon+1e-6)][:,0,0].to_series().dropna()
    res * M_TO_MM
    return res

def extract_climate(lat,lon,start,end,silo_path):
    time_period = slice(start,end)
    return extract_silo('daily_rain',time_period,lat,lon,silo_path),extract_silo('et_morton_wet',time_period,lat,lon,silo_path),

def run_for_location(lat,lon):
    rainfall, pet = extract_climate(lat,lon)
    if not len(rainfall) or not len(pet):
        print('Missing climate data at %f,%f'%(lat,lon))
        return None

    atmosph = read_atmosph_in(ORIG_MODEL)
    atmosph = atmosph[:len(rainfall)]
    atmosph.Prec = np.array(rainfall)
    atmosph.rSoil = np.array(pet)
    write_atmosph_in(atmosph,MODEL)
    outputs = run_hydrus(MODEL)
    return outputs
    # all in a temp directory
    # Extract climate (rainfall, PET) for cell
    # Put into atmosph.in
    
    # Extract soil properties (TODO)
    # Put into X.dat

    # Run model

    # res = extract results
    # return res

def extract_results(path):
    pass
    # Load X.out
    # return...

def ping(s):
    print(s,end='')
    sys.stdout.flush()

def min_threshold_filter(threshold):
    return lambda v,lat,lng: v>=threshold

def silo_climate_source(variable,silo_path):
    return lambda time_period,lat,lng: extract_silo(variable,time_period,lat,lng,silo_path)

class SpatialHPxRun(object):
    def __init__(self,
                 domain,
                 cell_filter=lambda v,lat,lng: True,
                 climate_sources={},):
        self.domain = domain
        self.cell_filter = cell_filter
        self.climate_sources = climate_sources

    def run(self,sim_start=DEFAULT_SIM_START,sim_end=DEFAULT_SIM_END):
        # mask = xr.open_rasterio(path)[0,:,:]
        coords_to_run = np.where(mask>threshold)
        dims = self.domain.dims

        self.coord_df = pd.DataFrame(dict(zip(dims,coords)))
        self.coord_df['lat'] = mask.indexes['y'][coord_df['y']]
        self.coord_df.lat = self.coord_df.lat.round(2)
        self.coord_df['lng'] = mask.indexes['x'][coord_df['x']]
        self.coord_df.lng = self.coord_df.lng.round(2)

        date_range = pd.date_range(sim_start,sim_end)

        def make_array(spec):
            if spec.get('time',False):
                shp = (len(date_range),mask.shape[0],mask.shape[1])
                return xr.DataArray(np.nan*np.ones(shp,dtype='f'),(date_range,mask.indexes['y'],mask.indexes['x']),('time','lat','lng'))
            return xr.DataArray(np.zeros(mask.shape,dtype='f'),(mask.indexes['y'],mask.indexes['x']),('lat','lng'))

        result_arrays = {o:make_array(spec) for o,spec in OUTPUTS.items()}
        results = xr.Dataset(result_arrays,{'time':date_range,'lat':mask.indexes['y'],'lng':mask.indexes['x']})
        for i,row in coord_df.iterrows():
            if (i % 100)==0: 
                ping('\n*, %d %s'%(i,str(datetime.now())))
            else:
                ping('.')
    #        print(i,row.y,row.x,row.lat,row.lng)
            res = run_for_location(row.lat,row.lng)
            if res is None: continue
            for o,spec in OUTPUTS.items():
                if spec.get('time',False):
                    results[o][:,int(row.y),int(row.x)] = np.nan
                    values = np.array(res[spec['file']][spec['column']])
                    results[o][:len(values),int(row.y),int(row.x)] = values
                else:
                    raise Exception('Not supported')
    #                results[o][int(row.y),int(row.x)] = 1
    #         if i > 3:
    #             break
            #break
        # initialise outputs / spatial domain
        # for each cell over threshold:
        #   res = run_for_location(lat,lon)
        #   for each output variable:
        #       # load_results
        return results
