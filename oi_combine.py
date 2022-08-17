#my_new_list = np.core.defchararray.add(tab_dist['prefix'][i], files)
from readcol import *

cd = True  ### Must always be FALSE when combined for the reconstruction (use false for CDIFF)
if cd == True:
    import oi_merge as oi_merge
else:
    import oi_merge2 as oi_merges#
import importlib
importlib.reload(oi_merge)

#### This routine combines several OIFITS files (snap-shots) into a unique file
data = 'calibrated_abdor_380.txt' ###OIFITS filenames to be combined
[files] = readcol(data, twod=False)
merged = oi_merge.oi_merge(files)
merged.write('COMB_JWST_SAM_abdor.fits') ### Output OIFITS filename
