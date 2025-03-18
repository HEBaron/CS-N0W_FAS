# Import the required packages
import numpy as np
import pandas as pd
import datetime as dt
import xarray as xr

# variables and filepaths to input
scen_nos = 3 # number of future scenarios
first_year = 2020 # first year of projections
start_year = 2025 # start year for which you want to produce future abstractions
end_year = 2026 # end year for which you want to produce future abstractions
netcdf_output = False
abstractions_filepath = 'file_path/England_Monthly_Abstractions_1km_GW_SW_TW_199901_201412.csv' # baseline abstractions
scaling_filepath = 'file_path/scaling.csv' # future scenario scaling csv file
grouping_filepath = 'file_path/grouping.csv' # grouping csv file (for collecting abstractions together into coarser water use sectors)
savepath='file_path/future_demand' # filepath for outputs

# assign WRZ (rows will have a null value if they're in a DCWW WRZ)
def get_spatial_match(X,Y,spatial_match,spatial_grid):
    sm = spatial_grid[spatial_match].sel(Northing=Y,Easting=X).values
    if np.isnan(sm):
        sm = np.nan
    else:
        sm = int(sm)
    return sm

# to process the scaling and grouping dataframes
class scaling_dataframe(pd.DataFrame):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
   
    def cols_scale(cls):
        return ['Primary Code','Secondary Code', 'Use Code']
    def cols_group(cls):
        return ['Primary Code','Secondary Code', 'Use Code','Source Code']
    # make the scaling dataframe more machine readable     
    def listify(self,value):
        if isinstance(value, str) and ';' in value:
            return value.split(';')
        else:
            return [value]        

    # return only the national ("global") scalings 
    def national(self,start_year,first_year):
        national_scale = self[~self[str(first_year)].astype(str).str.contains('\.csv')]
        national_scale.loc[:,str(start_year)] = national_scale[str(start_year)].astype(float)
        return national_scale

    # return only the spatial scalings 
    def spatial(self, scen_nos, start_year, first_year):
        spatial_filter = self[self[str(first_year)].astype(str).str.contains('\.csv')]
        spatial_scale = {}
        for x in range(0,spatial_filter.shape[0],scen_nos):
            scens = {}
            for y in range(scen_nos):
                scens.update({y: pd.read_csv(spatial_filter[str(first_year)].iloc[y+x])})
            spatial_scale.update({spatial_filter['Name'].iloc[x]: scens})
        return spatial_filter, spatial_scale

# to extract a workable dataframe from the CEDA format csv
class metadata_df():
    def __init__(self, filepath, metadata_lineno, col_names_lineno):
        self.filepath = filepath
        self.metadata_lineno = metadata_lineno
        self.col_names_lineno = col_names_lineno
        self.column_names = []
        self.df = None
        # Read metadata and initialize column names
        self._read_metadata()        
        # Read the actual data
        self._read_data()
        
    def _read_metadata(self):
        #Extracts column names from the metadata section.
        with open(self.filepath, 'r') as file:
            metadata_lines = [next(file).strip() for _ in range(self.col_names_lineno)]

        # Parse metadata to extract column names
        for i in range(self.metadata_lineno, self.col_names_lineno - 1, 2):
            long_name_line = metadata_lines[i]
            type_line = metadata_lines[i + 1]
            long_name = long_name_line.split(',')[2].strip()  # Extract the long_name (3rd column)
            self.column_names.append(long_name)
    
    #Read the actual data starting from col_names_lineno & skipping final line
    def _read_data(self):
        with open(self.filepath, 'r') as file:
            total_lines = sum(1 for _ in file)
        # Read the data, skipping metadata and the final line
        self.df = pd.read_csv(self.filepath, skiprows=self.col_names_lineno, nrows=total_lines - self.col_names_lineno - 2)  
        # Assign the parsed column names to the dataframe
        self.df.columns = self.column_names

    #Returns the processed dataframe.
    def get_dataframe(self):
        return self.df 
        
    #Returns the mean monthly abstractions for each row across the year range=(ave_year_start,ave_year_end).
    def average_abs(self,ave_year_start=2010,ave_year_end=2015):
        pattern = 'Monthly Abstraction (m3/month;'
        selected_columns = [col for col in self.df.columns if pattern in col and str(ave_year_start) <= col.split('-')[1] <= str(ave_year_end)]
        
        #Group by month (e.g., 'Jan', 'Feb', ...) and calculate the mean over the selected years
        monthly_averages = {}
        months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        
        for month in months:
            # Extract columns that match the month across year range
            monthly_columns = [col for col in selected_columns if col.startswith(f"{pattern} {month}-")]
            # Calculate the mean for these columns across rows
            if monthly_columns:
                monthly_averages[f'{month} (m3/month)'] = self.df[monthly_columns].mean(axis=1)

                
        non_monthly_columns = [col for col in self.df.columns if pattern not in col]
        #Create a new DataFrame with the monthly averages
        df_monthly_avg = pd.DataFrame(monthly_averages)
        df_final = pd.concat([ self.df[non_monthly_columns], df_monthly_avg], axis=1)
        return df_final  
        

        
# to apply the scaling (and grouping) to the abstractions
class abs_dataframe(pd.DataFrame):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
    #pick out abstraction columns 
    def unit_columns(self):
        pattern =  r'\(m3\/month\)'  # Regular expression pattern for columns with unit in the header
        unit_columns = self.filter(regex=pattern).columns.tolist()
        return unit_columns
        
    #using the filters described in scaling, create a filter condition     
    def sector_filter(self,scaling,cols):
        all_filters = None
        for index, row in scaling.iterrows():
            filter_condition = None
            for col in cols:
                if row[col][0] == "ALL":
                    continue
                else:
                    if filter_condition is None:
                        filter_condition = (self[col].isin(row[col]))
                    else:
                        filter_condition &= (self[col].isin(row[col]))
            if all_filters is None:
                all_filters = filter_condition
            else:
                all_filters |= filter_condition
        return all_filters

    #if the scaling factors are national, apply them to the various filters for each scenario
    def update_national(self, scaling, cols, year, scen_nos, scen):
        for x in range(0,scaling.shape[0],scen_nos):
            filter_condition = self.sector_filter(scaling.iloc[x:x+scen_nos],cols)
            factor = scaling[str(year)].iloc[x+scen]
            self.loc[filter_condition, self.unit_columns()] *= factor

    #if the scaling factors are spatial, apply them to the various filters for each scenario, matching on the columns provided
    def update_spatial(self, scaling, cols, factors, year, scen_nos, scen):
        for x in range(0,scaling.shape[0],scen_nos):
            scale_key = spatial_filter['Name'].iloc[x]
            factor = factors[scale_key][scen]
            match_columns = factor.columns[0]
            null_row = {col: (np.nan if col in match_columns else 1) for col in factor.columns} # add row with scaling factor 1 for abstractions with no WRZ in the spatial scaling files
            factor = pd.concat([factor,pd.DataFrame([null_row])],ignore_index=True)
            filter_condition = self.sector_filter(scaling.iloc[x:x+scen_nos],cols)
            merged_df = pd.merge(self,factor,on=match_columns,how='left')  
            merged_df.loc[filter_condition, self.unit_columns()] = merged_df.loc[filter_condition, self.unit_columns()].multiply(merged_df.loc[filter_condition, str(year)], axis=0)
            merged_df.drop([col for col in factor.columns if col not in match_columns], axis=1, inplace=True)
            self.__dict__.update(merged_df.__dict__)

    #output as csv
    def output_csv(self,year,scen_name,savepath):
        self.to_csv(f'{savepath}/abstractions_{year}_{scen_name}.csv',index=False)
        
    # group the abstractions according to the grouping file and output as netcdfs
    def group(self, grouping, cols, year, scen_name, scen_nos, mask,savepath):
        time_range = pd.date_range(start=dt.datetime(year, 1, 1), end=dt.datetime(year, 12, 1), freq='MS')
        dims = ('Northing','Easting','time')
        coords = {'Northing': mask['Northing'].values,'Easting': mask['Easting'].values,'time':  time_range}  
        all_filters = None
        for source in ['SW','GW']:
            ds = xr.Dataset(coords=coords)
            encoding = {}
            for x in range(0,grouping.shape[0],scen_nos):
                var = grouping['Name'].iloc[x]
                units = grouping['Units'].iloc[x]
                grid = np.zeros((len(mask['Northing']), len(mask['Easting']), len(time_range)))
                filter_condition = self.sector_filter(grouping.iloc[x:x+scen_nos],cols)
                filter_condition &= (self['Source Code'].isin([source]))                
                if all_filters is None:
                    all_filters = filter_condition
                else:
                    all_filters |= filter_condition
                agg_dict = {col: 'sum' for col in self.unit_columns()}
                grouped = self.loc[filter_condition].groupby(['OSGB Easting (m)', 'OSGB Northing (m)']).agg(agg_dict).reset_index()
                if grouped.shape[0]==0: continue
                for month_index, month in enumerate(self.unit_columns()):
                    for _, row in grouped.iterrows():
                        x_index = np.where(mask['Easting'] == row['OSGB Easting (m)'])[0][0]
                        y_index = np.where(mask['Northing'] == row['OSGB Northing (m)'])[0][0]
                        grid[y_index, x_index, month_index] = row[month] 
                grid[mask] = np.nan   
                encoding.update({var: {'_FillValue': np.nan}})
                ds[var] = xr.DataArray(grid, dims=dims, coords=coords)
                ds[var].attrs["units"] = units
                ds[var].attrs["standard_name"] = var
            ds.to_netcdf(f'{savepath}/abstractions_{source}_{year}_{scen_name}.nc', encoding=encoding)


# Extract dataframe of abstractions from CEDA csv
metadata_lineno = 21  # Line where metadata ends
col_names_lineno = 418  # Line where column name data ends
df_obj = metadata_df(abstractions_filepath, metadata_lineno, col_names_lineno)
abstractions = df_obj.average_abs()
abstractions['Use Code'] = abstractions['Use Code'].astype(str)

#read in and process the scaling factors 
scaling = pd.read_csv(scaling_filepath)
scaling = scaling_dataframe(scaling)
scaling[scaling.cols_scale()] = scaling[scaling.cols_scale()].map(scaling.listify) 
national_scale = scaling.national(start_year,first_year)
spatial_filter, spatial_scale = scaling.spatial(scen_nos, start_year,first_year)  
#read scenario names (used to name output files)
scen_names = scaling.Scenario[:scen_nos].to_list()

# create column(s) to match spatial scaling on (here we use WRZ - rows will have a null value if they're in a DCWW WRZ)
spatial_mappings = spatial_filter[~spatial_filter['Spatial Map'].isna()]
for index, row in spatial_mappings.iterrows():
    df = pd.read_csv(row[str(first_year)])
    spatial_match = df.columns[0]
    if spatial_match not in abstractions.columns:
        spatial_grid_filepath = row['Spatial Map']
        spatial_grid = xr.open_dataset(spatial_grid_filepath)
        abstractions[spatial_match] = abstractions.apply(lambda row: get_spatial_match(row['OSGB Easting (m)'],row['OSGB Northing (m)'],spatial_match,spatial_grid),axis=1)

if netcdf_output:
    #read in and process the grouping information (required to output netcdfs)
    grouping = pd.read_csv(grouping_filepath)
    grouping = scaling_dataframe(grouping)
    grouping[grouping.cols_group()] = grouping[grouping.cols_group()].map(grouping.listify)   
    # use your chosen mask to define the limits, resolution, and mask for the output netcdfs (if required)
    mask = xr.open_dataset(grouping['Mask'].iloc[0])
    mask =  mask['mask'].isnull()

# example lines for using the grouping functionality for the baseline abstractions
#my_df = abs_dataframe(abstractions.copy())
#my_df.group(grouping,grouping.cols_group(),2015,'baseline',scen_nos,mask,savepath)

#looping through year and scenarios, scale and output the future abstractions
for year in range(start_year,end_year+1):
    for scen in range(scen_nos):
        my_df = abs_dataframe(abstractions.copy())
        my_df.update_spatial(spatial_filter,scaling.cols_scale(),spatial_scale,year,scen_nos,scen)
        my_df.update_national(national_scale,scaling.cols_scale(),year,scen_nos,scen)
        my_df.output_csv(year,scen_names[scen],savepath)
        if netcdf_output:
            my_df.group(grouping,grouping.cols_group(),year,scen_names[scen],scen_nos,mask,savepath) 