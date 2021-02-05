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
from .gw import set_gw_depth

HP1_DIR=os.path.join(os.getcwd(),'HP1')
HP1_EXE=os.path.join(HP1_DIR,'hp1.exe')

DEFAULT_OUTPUTS={
    'Volume':{
        'time':True,
        'file':'T_LEVEL',
        'column':'Volume'
        
    },
    'pH':{
        'dims':[
            'time',
            'dist_x'
        ],
        'file':'nod_inf_chem',
        'column':'pH'
    },
    'pe':{
        'dims':[
            'time',
            'dist_x'
        ],
        'file':'nod_inf_chem',
        'column':'pe'
    },
    'k_pyrite':{
        'dims':[
            'time',
            'dist_x'
        ],
        'file':'nod_inf_chem',
        'column':'k_pyrite'
    },
    'dk_pyrite':{
        'dims':[
            'time',
            'dist_x'
        ],
        'file':'nod_inf_chem',
        'column':'dk_pyrite'
    },
    'mass_H2O':{
        'dims':[
            'time',
            'dist_x'
        ],
        'file':'nod_inf_chem',
        'column':'mass_H2O'
    },
    'm_O2':{
        'dims':[
            'time',
            'dist_x'
        ],
        'file':'nod_inf_chem',
        'column':'m_O2'
    }
    	          	    # mass_H2O	       Fe(2)	       Fe(3)	       S(-2)	        S(6)	          Na	           K	          Mg	          Ca	        C(4)	          Cl	          Al	        m_O2	    si_O2(g)	   si_CO2(g)	si_Ferrihydrite	  si_Fe(OH)2	si_H-Jarosite	si_Na-Jarosite	si_K-Jarosite	  si_calcite	 si_gibbsite	   si_gypsum	   si_pyrite	    pressure	   total_mol	      volume	     g_O2(g)	    g_CO2(g)	    k_pyrite	   	PressureHead	     siFerri	         gCl	          gK	         gCa	         gMg	         gNa	        gSO4	         gAl	         gS2	        gFe2	        gFe3	'
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

def min_threshold_filter(threshold,source):
    return lambda lat,lng: source(lat,lng)>=threshold

def max_threshold_filter(threshold,source):
    return lambda lat,lng: source(lat,lng)<=threshold

def and_filter(*filters):
    return lambda lat,lng: np.all([f(lat,lng) for f in filters])

def grid_parameters(grid,precision=3,format_str=None,numeric=False):
    lat_ys={}

    y_step = float(grid.y[5]-grid.y[0])/5
    y0 = grid.y[0]
    def y_i(lat):
        y = lat_ys.get(lat,None)
        if y is None:
            y = lat_ys[lat] = int((lat-y0)/y_step)
        return y

    lng_xs={}
    x_step = float(grid.x[5]-grid.x[0])/5
    x0 = grid.x[0]
    def x_i(lng):
        x = lng_xs.get(lng,None)
        if x is None:
            x = lng_xs[lng] = int((lng-x0)/x_step)

        return x

    if numeric:
        return lambda lat,lng:float(grid[y_i(lat),x_i(lng)])
#        return lambda lat,lng:float(grid.sel(y=lat,x=lng,method='nearest',tolerance=1e-4))

    if format_str is None:
        format_str = '%.'+str(precision)+'f'
#     return lambda lat,lng: format_str%grid.sel(y=lat,x=lng,method='nearest',tolerance=1e-4)
    return lambda lat,lng: format_str%grid[y_i(lat),x_i(lng)]

# def grid_parameters(grid,precision=3,format_str=None,numeric=False):
#     y_step = float(grid.y[5]-grid.y[0])/5
#     x_step = float(grid.x[5]-grid.x[0])/5
#     y0 = grid.y[0]
#     x0 = grid.x[0]
#     y_i = lambda lat: int(y0 + lat*y_step)
#     x_i = lambda lng: int(x0 + lng*x_step)

#     if numeric:
#         return lambda lat,lng:float(grid.sel(y=lat,x=lng,method='nearest',tolerance=1e-4))

#     if format_str is None:
#         format_str = '%.'+str(precision)+'f'
#     return lambda lat,lng: format_str%grid.sel(y=lat,x=lng,method='nearest',tolerance=1e-4)

class SpatialHPxRun(object):
    def __init__(self,
                 model,
                 domain,
                 cell_filter=lambda v,lat,lng: True,
                 climate_sources={},
                 parameters={},
                 gw_depth=lambda lat,lng: None,
                 outputs=DEFAULT_OUTPUTS):
        self.model = model
        self.domain = domain
        self.cell_filter = cell_filter
        self.climate_sources = climate_sources
        self.parameters = parameters
        self.gw_depth = gw_depth
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
        rows_to_keep = [self.cell_filter(lat,lng) for lat,lng in zip(df['lat'],df['lng'])]
        df = df[rows_to_keep]
#00        df = df[df.apply(lambda row: self.cell_filter(row['lat'],row['lng']),axis=1)]
        return df.reset_index()

    def _make_temp_model(self):
        res = tempfile.mktemp(prefix='hpx_',dir='.')
        shutil.copytree(self.model,res)
        return res

    def _remove_temp_model(self,tmp_model):
        assert tmp_model != self.model
        shutil.rmtree(tmp_model)

    def _run_for_location(self,time_period,lat,lng,dry_run=False):
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
        nan_param = False
        for k,v in parameters.items():
            if isinstance(v,str):
                continue

            if np.isnan(v):
                nan_param = True
                print(f'{k} is nan at {lat},{lng}. Skipping cell')
        if nan_param:
            return None

        if not 'ROOTDIR' in parameters:
            parameters['ROOTDIR'] = os.path.abspath(os.path.join(HP1_DIR,'..'))
        tmp_model = self._make_temp_model()
        try:
            write_atmosph_in(atmosph,tmp_model)
            substitute_templates(tmp_model,parameters)

            gw_depth = self.gw_depth(lat,lng)
            if gw_depth is not None:
                set_gw_depth(gw_depth,os.path.join(tmp_model,'Profile.dat'))

            if dry_run:
                return None

            outputs = run_hydrus(tmp_model)
            return outputs
        finally:
            self._remove_temp_model(tmp_model)
        return None

    def run_single(self,sim_start=DEFAULT_SIM_START,sim_end=DEFAULT_SIM_END,cell_number=0,dry_run=False):
        coords_to_run = self._find_run_coords()
        date_range = pd.date_range(sim_start,sim_end)

        row = coords_to_run.iloc[cell_number]
        res = self._run_for_location(date_range,row.lat,row.lng,dry_run)
        return res, {'lat':row.lat,'lng':row.lng}

    def extract_results(self,res):
        if res is None: return None
        processed = {}
        for o,spec in self.outputs.items():
            print(f'Processing output: {o}')
            if spec.get('time',False):
                try:
                    processed[o] = np.array(res[spec['file']][spec['column']])
                except:
                    print(res[spec['file']].columns)
                    raise
            elif 'dims' in spec:
                try:
                    data = res[spec['file']]
                    processed[o] = {
                        'dims':{d:np.array(data[d]) for d in spec['dims']},
                        'values':np.array(data[spec['column']])
                    }
                except:
                    print(res[spec['file']].columns)
                    raise
            else:
                raise Exception('Not supported')

        return processed

    def init_results(self,example,date_range):
        # What about other dimensions?
        all_dims = {
            'time':date_range,
            'lat':self.domain.indexes['y'],
            'lng':self.domain.indexes['x']
        }

        def make_array(lbl,spec):
            if spec.get('time',False):
                shp = (len(date_range),self.domain.shape[0],self.domain.shape[1])
                return xr.DataArray(np.nan*np.ones(shp,dtype='f'),(date_range,self.domain.indexes['y'],self.domain.indexes['x']),('time','lat','lng'))
            dims = spec.get('dims',[])
            if len(dims):
                dim_vals = [list(sorted(set(example[spec['column']]['dims'][d]))) for d in dims]
                stored_dims = []
                for i,d in enumerate(dims):
                    if d in all_dims:
                        if len(dim_vals[i]) == len(all_dims[d]):
                            stored = d
                        else:
                            stored = f'{d}_{spec["column"]}'
                    else:
                        stored = d

                    if not stored in all_dims:
                        all_dims[stored] = dim_vals[i]
                    stored_dims.append(stored)
                dims = stored_dims
                dims += ['lat','lng']
                dim_vals += [
                    self.domain.indexes['y'],
                    self.domain.indexes['x']                   
                ]
                shp = [len(dv) for dv in dim_vals]
                return xr.DataArray(np.nan*np.ones(tuple(shp),dtype='f'),tuple(dim_vals),tuple(dims))

            return xr.DataArray(np.nan*np.ones(self.domain.shape,dtype='f'),(self.domain.indexes['y'],self.domain.indexes['x']),('lat','lng'))

        result_arrays = {o:make_array(o,spec) for o,spec in self.outputs.items()}
        results = xr.Dataset(result_arrays,all_dims)
        return results

    def run(self,sim_start=DEFAULT_SIM_START,sim_end=DEFAULT_SIM_END,dry_run=False):
        # mask = xr.open_rasterio(path)[0,:,:]
        coords_to_run = self._find_run_coords()
        # dims = self.domain.dims
        date_range = pd.date_range(sim_start,sim_end)
        ping(f'Running {len(coords_to_run)} cells\n')

        results = None

        for i,row in coords_to_run.iterrows():
            if (i % 100)==0: 
                ping('\n*, %d %s'%(i,str(datetime.now())))
            else:
                ping('.')
    #        print(i,row.y,row.x,row.lat,row.lng)
            res = self._run_for_location(date_range,row.lat,row.lng,dry_run)

            if dry_run:
                continue

            if results is None:
                results = self.init_results(res,date_range)

            processed = self.extract_results(res)
            if processed is None: continue

            for o,spec in self.outputs.items():
                if spec.get('time',False):
                    results[o][:,int(row.y),int(row.x)] = np.nan
                    try:
                        values = processed[o]
                    except:
                        print(res[spec['file']].columns)
                        raise
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
