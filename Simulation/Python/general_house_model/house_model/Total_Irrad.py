#def irrad_any(rground,iday,LST,t):
'''
(scalar) gamma = azimuth angle of the surface,
 east:gamma = -90, west:gamma = 90
 south:gamma = 0, north:gamma = 180
(scalar) beta = inclination angle of the surface,
 horizontal: beta=0, vertical: beta=90
'''

import	numpy 		as np
import	qsun  
import	pandas 		as pd
# import	read_NEN
from utils import nen5060_to_dataframe


NUM = nen5060_to_dataframe(xl_tab_name="nen5060 - energie")

#Convert data frame to array
to_array=NUM.to_numpy()
# to_array=np.delete(to_array, 0, 0)  # removed patch

"""    
dom=to_array[:,2]          # day of month
hod=to_array[:,3]          # hour of day
qglob_hor=to_array[:,4]
qdiff_hor=to_array[:,5]
qdir_hor=to_array[:,6]
qdir_nor=to_array[:,7]
Toutdoor=to_array[:,8]/10
phioutdoor=to_array[:,9]
xoutdoor=to_array[:,10]/10
pdamp=to_array[:,11]                            # wrong assignment
vwind=to_array[:,12]/10                         # at 10 m height
dirwind=to_array[:,13]
cloud=to_array[:,14]/10
rain=to_array[:,15]/10                          # wrong assignment
"""

# assignment to numpy arrays directly from dataframe columns
dom = NUM.loc[:,'DAY(datum)'].values                           # day of month
hod = NUM.loc[:,'HOUR(uur)'].values                            # hour of day
qglob_hor = NUM.loc[:,'globale_zonnestraling'].values          # global irradiation
qdiff_hor = NUM.loc[:, 'diffuse_zonnestraling'].values         # diffuse irradiation
qdir_hor = NUM.loc[:, 'directe_zonnestraling'].values          # direct irradiation
qdir_nor = NUM.loc[:, 'directe_normale_zonnestraling'].values  # DNI
Toutdoor = NUM.loc[:, 'temperatuur'].values / 10.0             # temperature
phioutdoor = NUM.loc[:, 'relatieve_vochtigheid'].values        # %RH
xoutdoor = NUM.loc[:, 'absolute_vochtigheid'].values / 10.0    # AH

rain = NUM.loc[:, 'neerslaghoeveelheid'].values / 10.0         # precipitation
vwind = NUM.loc[:, 'windsnelheid'].values / 10.0               # wind speed
dirwind = NUM.loc[:, 'windrichting'].values                    # wind dir
cloud = NUM.loc[:, 'bewolkingsgraad'].values / 8.0             # cloud coverage
sunduration = NUM.loc[:, 'zonneschijnduur'].values / 10.0       # daily solar duration
pdamp = NUM.loc[:, 'dampspanning'].values                      # vapour pressure

#______________##__________________________


'''
(scalar) gamma = azimuth angle of the surface,
    east:gamma = -90 (270), west:gamma = 90
    south:gamma = 0, north:gamma = 180
    (scalar) beta = inclination angle of the surface,
    horizontal: beta=0, vertical: beta=90
'''

#time : hour of the year
t=(np.array(list(range(1,8761)))-1)*3600

#ground reflection is ignored
rground=0 

#day of the year 
iday = 1+np.floor(t/(24*3600))
#local time in hour : from 0 to 23:00
LST  = np.floor((t/3600) % 24)

#Define an empty matrix
E= np.zeros((8760,9,4))


'# -90 (E), -45 (SE), 0 (S), 45 (SW), 90 (W), 135 (NW), 180 (N), 225 (NE)'

k=-1

for j in range(9):
    if k<7:
        
        gamma=45*(k-1)    #gamma -90 (E), -45 (SE), 0 (S), 45 (SW), 90 (W), 135 (NW), 180 (N), 225 (NE)
        beta=90
        
    else:
            
        gamma=90
        beta=0
    
    k=k+1
    
    for i in range(8760):
        temp = qsun.Irrad(qdiff_hor[i],qdir_nor[i],gamma,beta,rground,iday[i],LST[i])
        #E[i][j]=qsun(t[i],qdiff_hor[i],qdir_nor[i],gamma,beta,rground)
        E[:,j][i]=temp
    #n=n+1
#myarray = np.asarray(E)
myarray = E  

#add time and Total radiation together

qsunE=np.vstack((t,myarray[:,0,2]))
qsunSE=np.vstack((t,myarray[:,1,2]))
qsunS=np.vstack((t,myarray[:,2,2]))
qsunSW=np.vstack((t,myarray[:,3,2]))
qsunW=np.vstack((t,myarray[:,4,2]))
qsunNW=np.vstack((t,myarray[:,5,2]))
qsunN=np.vstack((t,myarray[:,6,2]))
qsunNE=np.vstack((t,myarray[:,7,2]))
qsunhor=np.vstack((t,myarray[:,8,2]))