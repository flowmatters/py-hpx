import pandas as pd
import numpy as np
import xarray as xr
import sys
import os
from glob import glob
import re
from datetime import datetime
import tempfile
import shutil

HP1_DIR=os.path.join(os.getcwd(),'HP1')
HP1_EXE=os.path.join(HP1_DIR,'hp1.exe')

DEFAULT_OUTPUTS={
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
        'skiprows':[0,1,2,3,4,6,7,8],
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

def hpx_exe(fn=None):
    global HP1_EXE
    if fn is not None:
        HP1_EXE=fn
    return HP1_EXE

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

def substitute_templates(dest,parameters={},default=-9999):
    templates = glob(os.path.join(dest,'*.tpl'))
    for tpl in templates:
        dest = tpl[:-4]
        with(open(tpl,'r')) as fp:
            tpl_txt = fp.read()
        for k,v in parameters.items():
            key = f'${k}$'
            tpl_txt = tpl_txt.replace(key,v.rjust(len(key)))
        with(open(dest,'w')) as fp:
            fp.write(tpl_txt)

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
        cmd_line = f'{HP1_EXE} .'
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
    output_filenames = list(set(glob(os.path.join(model,'*.out'))).union(set(glob(os.path.join(model,'*.OUT')))))
    result = {}
    for fn in output_filenames:
        base = os.path.basename(fn)
        try:
            result[base[:-4]] = read_output_file(fn)
        except:
            pass
            #print('Error reading %s'%base)


    return result

def extract_results(path):
    pass
    # Load X.out
    # return...

def ping(s):
    print(s,end='')
    sys.stdout.flush()

def min_threshold_filter(threshold):
    return lambda v,lat,lng: v>=threshold

def grid_parameters(grid,precision=3,format_str=None):
    if format_str is None:
        format_str = '%.'+str(precision)+'f'
    return lambda lat,lng: format_str%grid.sel(y=lat,x=lng,method='nearest',tolerance=1e-4)

class SpatialHPxRun(object):
    def __init__(self,
                 model,
                 domain,
                 cell_filter=lambda v,lat,lng: True,
                 climate_sources={},
                 parameters={},
                 outputs=DEFAULT_OUTPUTS):
        self.model = model
        self.domain = domain
        self.cell_filter = cell_filter
        self.climate_sources = climate_sources
        self.parameters = parameters
        self.outputs=outputs

    def _find_run_coords(self):
        dim_translate = {
            'y':'lat',
            'x':'lng'
        }
        all_coords = np.where(np.logical_not(np.isnan(self.domain)))
        dims = self.domain.dims
        coordinate_table = dict(zip(dims,all_coords))
        df = pd.DataFrame(coordinate_table)
        df['v'] = self.domain.isel(x=xr.DataArray(np.array(df.x),dims='z'),y=xr.DataArray(np.array(df.y),dims='z'))
        for d in dims:
            coord = dim_translate.get(d,d)
            df[coord] = self.domain.coords[d][np.array(df[d])]
        df = df[df.apply(lambda row: self.cell_filter(row['v'],row['lat'],row['lng']),axis=1)]
        return df.reset_index()

    def _make_temp_model(self):
        res = tempfile.mktemp(prefix='hpx_',dir='.')
        shutil.copytree(self.model,res)
        return res

    def _remove_temp_model(self,tmp_model):
        assert tmp_model != self.model
        shutil.rmtree(tmp_model)

    def _run_for_location(self,time_period,lat,lng):
        climate_inputs = {k:getter(time_period,lat,lng) for k,getter in self.climate_sources.items()}
        empty_climate=False
        for k,vals in climate_inputs.items():
            if not len(vals):
                print(f'Missing {k} data at {lat},{lng}')
                empty_climate=True
        if empty_climate:
            print(f'Skipping {lat},{lng}')
            return None
        
        atmosph = read_atmosph_in(self.model)
        atmosph = atmosph[:len(time_period)]
        for k,v in climate_inputs.items():
            atmosph[k] = np.array(v)

        parameters = {k:getter(lat,lng) for k,getter in self.parameters.items()}
        parameters['HPXDIR'] = os.path.abspath(os.path.join(HP1_DIR,'..'))
        tmp_model = self._make_temp_model()
        try:
            write_atmosph_in(atmosph,tmp_model)
            substitute_templates(tmp_model,parameters)
            outputs = run_hydrus(tmp_model)
            return outputs
        finally:
            self._remove_temp_model(tmp_model)
        return None

    def run(self,sim_start=DEFAULT_SIM_START,sim_end=DEFAULT_SIM_END):
        # mask = xr.open_rasterio(path)[0,:,:]
        coords_to_run = self._find_run_coords()
        dims = self.domain.dims
        date_range = pd.date_range(sim_start,sim_end)

        ping(f'Running {len(coords_to_run)} cells\n')

        def make_array(spec):
            if spec.get('time',False):
                shp = (len(date_range),self.domain.shape[0],self.domain.shape[1])
                return xr.DataArray(np.nan*np.ones(shp,dtype='f'),(date_range,self.domain.indexes['y'],self.domain.indexes['x']),('time','lat','lng'))
            return xr.DataArray(np.zeros(self.domain.shape,dtype='f'),(self.domain.indexes['y'],self.domain.indexes['x']),('lat','lng'))

        result_arrays = {o:make_array(spec) for o,spec in self.outputs.items()}
        results = xr.Dataset(result_arrays,{'time':date_range,'lat':self.domain.indexes['y'],'lng':self.domain.indexes['x']})
        for i,row in coords_to_run.iterrows():
            if (i % 100)==0: 
                ping('\n*, %d %s'%(i,str(datetime.now())))
            else:
                ping('.')
    #        print(i,row.y,row.x,row.lat,row.lng)
            res = self._run_for_location(date_range,row.lat,row.lng)
            if res is None: continue
            for o,spec in self.outputs.items():
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
