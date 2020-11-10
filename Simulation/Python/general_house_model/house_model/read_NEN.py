import numpy as np
from pathlib import Path
import os
import pandas as pd

data_dir = Path.cwd() /'data'
#data_dir = Path.cwd() /'data'
output_dir = Path.cwd()/'working'/'submit'
NENdata_path = data_dir/'NEN5060-2018.xlsx'
print(NENdata_path)
xls = pd.ExcelFile(NENdata_path)    