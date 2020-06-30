import xarray as xr
import string
from glob import glob
import os

MM_DAY_TO_M_YR = 1e-3 * 365.25

SILO_CLIMATE_VARIABLES={
    'Prec':'daily_rain',
    'rSoil':'et_morton_wet'
}

CC_LAT_CELL_SIZE=1.25
CC_LON_CELL_SIZE=1.875

CLIMATE_CHANGE_SCALING_VARIABLES={
    'Prec':'pr',
    'rSoil':'wvap'
}
CLIMATE_CHANGE_SCALING_FILENAMES={
    'Prec':'pr',
    'rSoil':'wvap-morton'
}
CC_FN_TEMPLATE = string.Template('../spatial-inputs/climate-change/${model}/${scenario}/${variable}/${variable}_Amon_${model}_${scenario}_r1i1p1_${change_timeperiod}--change-wrt-1986-2005-seasavg-clim_native_AUS.nc')
CC_MONTHS=[
    'january',
    'february',
    'march',
    'april',
    'may',
    'june',
    'july',
    'august',
    'september',
    'october',
    'november',
    'december'
]


def extract_silo(var,time_period,lat,lon,silo_path):
    ds = xr.open_mfdataset(glob(os.path.join(silo_path,'*.%s.nc'%var)))
    res = ds[var].loc[time_period,(lat-1e-6):(lat+1e-6),(lon-1e-6):(lon+1e-6)][:,0,0].to_series().dropna()
    # res * M_TO_MM
    return res

def extract_climate(lat,lon,start,end,silo_path):
    time_period = slice(start,end)
    return extract_silo('daily_rain',time_period,lat,lon,silo_path),extract_silo('et_morton_wet',time_period,lat,lon,silo_path),

def silo_climate_source(variable,silo_path,scale=1):
    return lambda time_period,lat,lng: extract_silo(variable,time_period,lat,lng,silo_path) * scale

def scaled_climate_source(hp_var,gcm,rcp,change_timeperiod,historical_data):
    gcm_var = CLIMATE_CHANGE_SCALING_VARIABLES[hp_var]
    gcm_fn = CLIMATE_CHANGE_SCALING_FILENAMES[hp_var]
    fn = CC_FN_TEMPLATE.substitute(variable=gcm_fn,scenario=rcp,model=gcm,change_timeperiod=change_timeperiod)
    scaling_data = xr.open_dataset(fn,decode_times=False)
    monthly_factors = [scaling_data[f'{gcm_var}_{mth}'][0] for mth in CC_MONTHS]
    dims = ['month','lat','lon']
    coords = {
        'month':range(1,13),
        'lat':monthly_factors[0]['lat'][...],
        'lon':monthly_factors[0]['lon'][...]
    }
    monthly_scaling = (xr.DataArray(monthly_factors,dims=tuple(dims),coords=coords) * 0.01) + 1.0

    def generate(time_period,lat,lng):
        hist = historical_data(time_period,lat,lng)
        lat_delta=(CC_LAT_CELL_SIZE-1e-2)/2
        lon_delta=(CC_LON_CELL_SIZE-1e-2)/2
        scale_series = monthly_scaling.loc[time_period.month,(lat-lat_delta):(lat+lat_delta),(lng-lon_delta):(lng+lon_delta)]
        assert scale_series.shape[1]==1
        assert scale_series.shape[2]==1
        scale_series = scale_series[:,0,0]
        return hist * scale_series
    
    return generate

def silo_climate(silo_path,scale=MM_DAY_TO_M_YR):
  return {k:silo_climate_source(v,silo_path,scale) for k,v in SILO_CLIMATE_VARIABLES.items()}

def scaled_climate(gcm,rcp,change_timeperiod,historical_data):
  return {k:scaled_climate_source(k,gcm,rcp,change_timeperiod,v) for k,v in historical_data.items()}


