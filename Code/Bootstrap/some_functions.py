import numpy as np
from scipy.spatial.distance import squareform

def dyn_corr_isc(D, slide,verbose=True,summary_statistic=None):
    """Intersubject correlation
    For each moment in time compute the Pearson correlation between 
    the time series of length slide for each pair of voxels or ROI's, 
    defining a matrix correlation. If summary_statistic is None, return 
    t values (total time - slide) x R*(R-1)/2 values of Pearson correlation
    values for each subject. Alternatively, supply either
    np.mean or np.median to compute summary statistic (Fisher Z will
    be applied and inverted if using mean).
        
    The implementation is based on Scha professor's code
    
    Parameters
    ----------
    data : list or ndarray
        fMRI data for which to compute ISC
        
    pairwise : bool, default: False
        Whether to use pairwise (True) or leave-one-out (False) approach
        
    summary_statistic : None
        Return all ISCs or collapse using np.mean or np.median
    Returns
    -------
    iscs : subjects or pairs by pair of voxels ndarray
        ISC for each subject or pair (or summary statistic) per pair of voxels
    """
    
    # Infer subjects, TRs, voxels and print for user to check
    n_subjects = D.shape[2]
    n_TRs = D.shape[0]
    n_voxels = D.shape[1]
    
    
    print(f"Assuming {n_subjects} subjects with {n_TRs} time points "
          f"and {n_voxels} voxel(s) or ROI(s).\n")
          #f"Will compute sliding window analysis with a window length of -{neg_win} and +{pos_win} samples.")
    
    t_wind=n_TRs+1-slide
    corr_din=np.zeros((t_wind, int(n_voxels*( n_voxels-1)/2),n_subjects))        
    
    for s in range(n_subjects):
        if verbose:
            progress = 100 * (s/n_subjects)
            sys.stdout.write("\r%d%%" % progress)
            sys.stdout.flush() 
        for t in range(t_wind):
            D_new=np.corrcoef(D[t:t+slide,:,s].T)-np.identity(n_voxels)
            corr_din[t,:,s]=squareform(D_new,checks=False)
    
    if summary_statistic == np.mean:
        corr_din = np.tanh(summary_statistic(np.arctanh(corr_din), axis=2))
    elif summary_statistic == np.median:    
        corr_din = summary_statistic(corr_din, axis=2)
    elif not summary_statistic:
        pass
    else:
        raise ValueError("Unrecognized summary_statistic! Use None, np.median, or np.mean.")
    
    return corr_din

# ## Test statistic

#input: 74x193
def phase_shift(X,rand_state=None):
    r,t = X.shape
    Y=np.fft.rfft(X)
    amp=np.abs(Y)
    angle = np.angle(Y)

    prng= np.random.RandomState(rand_state)

    random_vec= prng.rand(1,angle.shape[1])*2*np.math.pi 
    angle2 = angle + random_vec
    amp2= amp*(np.cos(angle2)+1j*np.sin(angle2))

    Y2= np.fft.irfft(amp2,n=t) 
    
    return(Y2)

def surrogate_dyn_corr(data, slide, cut=100,name_file='',rand_name=True, rand_state=None):
	"""
	Phase randomization for dyn_corr

	For each pair of ROIs compute a null distribution of 
	dyn_corr where time series for each subject are phase 
	randomized prior to computing dyn_corr. 

	Input:- data: Matrix of data. Size = (time_points, regions, and n_subjects)
	  - slide: Size of window for corr_dyn (size in second= slide*TR)
	- - n_shifts: Number of copies of surrogate data. Default set to 1000
	  - cut: Cut proces for cut subjects everytime for generate corr_dyn 
	  (to not crash my computer). Default set to 100
	  - file: If not empty, the oupt will be stored in this file (omit.npy). Otherwise return the result
	  - rand_state: rand_state  
	At the end we generate n_shifts copies for dyn_corr output
	"""
	n_TRs = data.shape[0]
	n_regions = data.shape[1]
	n_subjects = data.shape[2]

	surrogate=[]
	for i in range(n_subjects):
		X=data[:,:,i].T
		Y=phase_shift(X,rand_state*n_subjects + i)
		surrogate.append(Y.T)

	s_data=np.stack(surrogate,axis=2)    

	r=n_subjects % cut
	if r == 0:
		k = int(n_subjects/cut)
	else:
		k = int(n_subjects/cut) + 1

	slidings=[]    
	for i in range(k):
		sliding=dyn_corr_isc(s_data[:,:,i*cut:(i+1)*cut],slide=slide,verbose=False,summary_statistic=np.mean)
		slidings.append(sliding)

	sliding_all=np.mean(np.stack(slidings,axis=2),axis=2)

	if len(name_file) == 0:
		return(sliding_all)
	else:
		if rand_name:
			if rand_state<10:
				name=name_file+'_00'+str(rand_state)+'.npy'	
			elif rand_state<100:
				name=name_file+'_0'+str(rand_state)+'.npy'
			else:
				name=name_file+'_'+str(rand_state)+'.npy'
		else:
			name=name_file

	np.save(name,sliding_all)
