
import pandas as pd
from pathlib import Path

def read_nen_weather(start_date, end_date):  # date in YYYYMD format
    """Reads the weather data from the NEN 5060 standard.

    Reads the CSV files provided from the NEN 5060 standard year. This function reads the weather data from
    a chosen start date to the chosen end date and adds the year, month and day columns together
    into the YYYYMD format.

    Args:
        start_date (int):   the starting point for which the weather data will be gathered in the array
        end_date (int):     the end point for which the weather data will be gathered in the array

    Returns:
        Array containing 26 columns with weather data, most relevant YYYYMMDD, HH, T:

        - YYYYMD (int): Date (YYYY=year,M=month,D=day).
        - H (int): Time in hours, the hourly division 05 runs from 04.00 UT to 5.00 UT.
        - T (float) Temperature at 1.50 m at the time of observation (°C).
    """

    directory_prefix = ''  # WeatherData\ for main directory
    nen_weather_data = pd.read_csv(directory_prefix + 'NEN5060-A2a.csv',
                                   header=5, sep=r'\s*,\s*', engine='python')

    # Read the data and convert into an array from the start till end date
    nen_weather_data["YYYYMD"] = (nen_weather_data['Y'].astype(str) + nen_weather_data['M'].astype(str)
                                  + nen_weather_data['D'].astype(str))
    nen_weather = nen_weather_data.set_index("YYYYMD", drop=False)
    nen_weather_range = nen_weather.loc[str(start_date):str(end_date), :]

    return nen_weather_range


def nen5060_to_dataframe(xl_tab_name: str="nen5060 - energie") -> pd.DataFrame :
    """ conversion from NEN5060 spreadsheet tab into Dataframe.

    Args:
        xl_tab_name: (str) tabname from NEN5060 spreadsheet ("nen5060 - energie", "ontwerp 1%" or "ontwerp 5%")

    Returns:
        pandas Dataframe with contents of NEN5060 tabsheet

    """
    # print(Path.cwd())
    data_dir = Path.cwd() / 'NEN_data'
    output_dir = Path.cwd()/'working'/'submit'
    NENdata_path = data_dir / 'NEN5060-2018.xlsx'
    print(NENdata_path)
    xls = pd.ExcelFile(NENdata_path)
    print(xls.sheet_names)  # Check sheet names

    # select sheet "nen5060 - energie" by NEN default
    df5060 = pd.DataFrame()
    # df5060 = pd.read_excel(xls, 'nen5060 - energie')  # this file is part of NEN 5060 20018
    # NEN5060-2018.xlsx has two lines with column headers
    # first line is column name, second line is measurement unit
    df5060 = pd.read_excel(xls, xl_tab_name, header=[0,1])  # this file is part of NEN 5060 20018
    ind = df5060.index
    print(ind.values)
    print(df5060.head())
    print(df5060.columns)

    return df5060 # pandas Dataframe

