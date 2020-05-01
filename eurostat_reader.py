import zipfile
import glob
import os
from functools import reduce

import pandas as pd
import numpy as np


def eurostat_reader(file_path, na_proportion=0.8,
                    cols_to_drop={"all": {'UNIT'},
                                  "arrivals_at_tourist_accommodation_establishments.zip": {'C_RESID'},
                                  "average_length_of_stay_at_hospitals.zip": {'SEX', 'AGE', 'ICD10'},
                                  "crude_death_rate.zip": {'SEX', 'AGE'},
                                  "disposable_income_per_inhabitant.zip": {'DIRECT'},
                                  "percentage_education_attainment.zip": {'SEX', 'AGE'},
                                  "population_numbers.zip": {'SEX', 'AGE'},
                                  "hospital_discharges_per_100k_inhabitants.zip": {'AGE', 'INDIC_HE', 'SEX'}}):

    with zipfile.ZipFile(file_path, 'r') as zip_file:
        for file in zip_file.namelist():
            if 'Data' in file:
                with zip_file.open(file) as data_file:
                    df = pd.read_csv(data_file, encoding="ISO-8859-1")

    base_file_path = os.path.basename(file_path)
    # Drop desired columns
    if base_file_path in cols_to_drop.keys():
        if (len(cols_to_drop[base_file_path].intersection(df.columns)) > 0):
            df = df.drop(columns=cols_to_drop[base_file_path].intersection(df.columns))
    df = df.drop(columns=cols_to_drop["all"])

    # Clean Value column and convert to float
    df["Value"] = df["Value"].apply(
        lambda x: x.replace(",", "")).replace({":": None}, regex=False).astype(np.float32)

    # Check if there is no data for certain years and then drop those years
    na_check = df.set_index('TIME')['Value'].isna().all(level=0)
    if len(na_check.index[na_check]) > 0:
        df = df.set_index('TIME').drop(list(na_check.index[na_check])).reset_index()

    # Find the column containing the relevant values to spread on
    value_col = [col for col in df.columns if col not in {"TIME", "GEO", "Value"}]

    if len(value_col) > 1:
        raise ValueError(f"Too many columns available to spread on for file '{file_path}', "
                         f"namely {value_col}. Check the data and add columns to remove to "
                         "cols_to_drop.")

    # Pivot from long to wide
    if len(value_col) == 1:
        df = df.pivot_table(index=['TIME', 'GEO'], columns=value_col[0],
                            values='Value').reset_index()
        try:
            del df.columns.name
        except AttributeError:
            pass

    # Select only the latest data for which all data is known
    test = (df.groupby('TIME').apply(lambda x: len(x) - x.isna().sum()) > 0).all(axis=1)
    max_year = test[test].index.max()
    df = df.set_index("TIME").loc[max_year].reset_index()

    # Drop fully NA rows and columns
    df = df.dropna(how='all').dropna(axis=1, how='all')
    na_props = df.isna().sum().divide(df.apply(len))

    if len(na_props[na_props > na_proportion].index) > 0:
        print(f"Over {na_proportion*100} percent NAs in the column(s) "
              f"{na_props[na_props > na_proportion].index} in file {base_file_path}. "
              "Dropping these columns")

        df = df.drop(columns=na_props[na_props > na_proportion].index)

    if len(value_col) == 0:
        df = df.rename(columns={'Value': base_file_path.replace('.zip', '')})

    return df


dfs = []
keep_time=False

for file in glob.glob("data/eurostat/*.zip"):
    if 'by_rail_by_loading_unloading_region' in file:
        continue
    if 'TOCLEAN' in file:
        # TODO: Skip these files for now. These are per NUTS3 region and need to be aggregated.
        continue

    df = eurostat_reader(file)

    if keep_time:
        dfs.append(df)
    else:
        dfs.append(df.drop(columns=["TIME"]))

# Merge (outer join) the list of DataFrames into one big DataFrame
if keep_time:
    df_merged = reduce(lambda x, y: pd.merge(x, y, on=['TIME', 'GEO'], how='outer'),
                       dfs).sort_values(by=['TIME', 'GEO'])
else:
    df_merged = reduce(lambda x, y: pd.merge(x, y, on=['GEO'], how='outer'),
                       dfs).sort_values(by=['GEO'])

# Drop extra regions
extra_region = (df_merged["GEO"] ==
                "Extra-Regio NUTS 1") | (df_merged["GEO"] == "Extra-Regio NUTS 2")
df_merged = df_merged.drop(extra_region[extra_region].index)

# Propagate constant values over years for each region - Total area
if keep_time:
    df_merged = df_merged.set_index('GEO')
    constant_columns = ['Total area']

    for region in df_merged.index.unique():
        for col in constant_columns:
            value = df_merged.loc[region, col][df_merged.loc[region, col].notna()]

            try:
                df_merged.loc[region, col] = [value for i in range(len(df_merged.loc[region, col]))]
            except ValueError:
                print(f"Region: {region}, Value: {value}")

    df_merged = df_merged.reset_index()

dict_rename = {
    'GEO': 'region',
    'Freight and mail loaded': 'airFreightLoaded',
    'Freight and mail unloaded': 'airFreightUnloaded',
    'Passengers carried (arrival)': 'airPassengersArrived',
    'Passengers carried (departures)': 'airPassengersDeparted',
    'Total area': 'area',
    ('Hotels; holiday and other short-stay accommodation; camping grounds, recreational vehicle '
     'parks and trailer parks'): 'touristArrivals',
    'In-patient average length of stay (in days)': 'hospitalStay',
    'broadband_access': 'broadbandAccess',
    'Malignant neoplasms (C00-C97)': 'deathRateCancer',
    'Diabetes mellitus_x': 'deathRateDiabetes',
    'Ischaemic heart diseases': 'deathRateChd',
    'Influenza (including swine flu)': 'deathRateInfluenza',
    'Pneumonia_x': 'deathRatePneumonia',
    'Disposable income, net': 'disposableIncome',
    'gdp_per_inhabitant': 'gdpPerInhabitant',
    'Medical doctors': 'medicalDoctors',
    'Nurses and midwives': 'nursesMidwives',
    'Available beds in hospitals (HP.1)': 'availableBeds',
    'Curative care beds in hospitals (HP.1)': 'curativeCareBeds',
    'Long-term care beds in hospitals (HP.1)': 'longtermCareBeds',
    'Other beds in hospitals (HP.1)': 'otherBeds',
    'Psychiatric care beds in hospitals (HP.1)': 'psychiatricCareBeds',
    'Rehabilitative care beds in hospitals (HP.1)': 'rehabilitativeCareBeds',
    'Diabetes mellitus_y': 'dischargeRateDiabetes',
    'Neoplasms': 'dischargeRateCancer',
    'Hypertensive diseases': 'dischargeRateHypertension',
    'Diseases of the respiratory system (J00-J99)': 'dischargeRateRespiratory',
    'Other ischaemic heart disease': 'dischargeRateChd',
    'Pneumonia_y': 'dischargeRatePneumonia',
    'Tuberculosis': 'dischargeRateTB',
    'Internet use: interaction with public authorities (last 12 months)': 'internetContactAuthorities',
    'Electrified railway lines': 'lengthElectrifiedRailway',
    'Motorways': 'lengthMotorways',
    'Navigable canals': 'lengthCanals',
    'Navigable rivers': 'lengthRivers',
    'Other roads': 'lengthOtherRoads',
    'Railway lines with double and more tracks': 'lengthLargeRailway',
    'Total railway lines': 'lengthRailway',
    'Freight loaded': 'maritimeFreightLoaded',
    'Freight unloaded': 'maritimeFreightUnloaded',
    'Passengers disembarked': 'maritimePassengersDisembarked',
    'Passengers embarked': 'maritimePassengersEmbarked',
    'Median age of population': 'medianAge',
    'Less than primary, primary and lower secondary education (levels 0-2)': 'lowerEducation',
    'Upper secondary and post-secondary non-tertiary education (levels 3 and 4)': 'higherEducation',
    'Tertiary education (levels 5-8)': 'tertiaryEducation',
    'population_density': 'populationDensity',
    'population_numbers': 'populationNumbers',
    'risk_of_poverty_or_social_exclusion': 'riskOfPovertyOrSocialExclusion'
}

if keep_time:
    dict_rename['TIME'] = 'time'

df_merged = df_merged.rename(columns=dict_rename)
df_merged.to_csv("data/merged_eurostat.csv", index=False)

# Inter-city railroad connections
# file_path = ("data/eurostat/passengers_by_rail_by_loading_unloading_region.zip")

# with zipfile.ZipFile(file_path, 'r') as zip_file:
#     for file in zip_file.namelist():
#         if 'Data' in file:
#             with zip_file.open(file) as data_file:
#                 df = pd.read_csv(data_file, encoding="ISO-8859-1")

# # GEO is constant, simply "Italy"
# df = df.drop(columns=['UNIT', 'GEO'])

# # Clean Value column and convert to float
# df["Value"] = df["Value"].apply(
#     lambda x: x.replace(",", "")).replace({":": None}, regex=False).astype(np.float32)

# # Check if there is no data for certain years and then drop those years
# na_check = df.set_index('TIME')['Value'].isna().all(level=0)
# if len(na_check.index[na_check]) > 0:
#     df = df.set_index('TIME').drop(list(na_check.index[na_check])).reset_index()

# # Pivot from long to wide
# df = df.pivot_table('Value', index=['TIME', 'C_LOAD'], columns='C_UNLOAD').reset_index()
# del df.columns.name

# df.to_csv("data/interregion_railroad_travel.csv", index=False)
