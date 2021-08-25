#Generate 1000 samples of surrogate data using the function surrogate_dyn_corr
#and varying the initial random state from 0 to 999

import numpy as np
import time
import multiprocessing
from scipy.spatial.distance import squareform
import some_functions
from some_functions import surrogate_dyn_corr , phase_shift, dyn_corr_isc

X=np.load("../Data/ts_data_regions.npy" ) 

def surrogate_for_pool(rand):
	print("Random state %s" % rand)
	surrogate_dyn_corr(X,slide=18,cut=100,name_file='../Data/surrogate/sliding_all',rand_name=True, rand_state=rand)
	
if __name__ == '__main__':
	start_time = time.time()
	pool = multiprocessing.Pool()
	pool.map(surrogate_for_pool, range(900,1000))
	end_time = time.time()
	print("--- %s seconds ---" % (end_time- start_time))
