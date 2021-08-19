# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.6.0
#   kernelspec:
#     display_name: mapper
#     language: python
#     name: mapper
# ---

# ## Labeling MNI
# <ol>
#     <li> ball_set
#     <li> coord2coord
#     <li> mni2region 
#     <li> mni2cl_region
#     <li> set_dictionaries_rois

import json
#------------------ set_dictionaries_rois --------------
#Set the dictionaries and the list of attributes for 
#function r_info
#Input: - file_info: file with json dictionaries
#       - file_nodes: file with json  attributes
#Output: - r_info: with those file sets 
#----------------- r_info -----------------
#Gives information of nodes of the given Atlas (Shen for me)
#Input: i: number of the region
#Output result: dictionarie with the MNI coordinates,
#      lobe, network and Brodmann area of region i

import itertools as it


#------------------ ball_set -----------------------
def ball_set(center,r):
    """ Generates the ball of integer coordinates of a ball with 
    integer center and radius r
    Input: - center: center of ball
           - r: radius of ball    
    Output:- set of coordinates into the ball
    """
    center=np.array(center)
    l=np.arange(-(r*r),(r*r)+1,1)
    lista=[]
    for element in it.product(l,l,l):
        norm=np.linalg.norm(element)
        if (norm*norm) <= (r*r):
            lista=lista+[center+element]
    return(lista)


#----------------------coord2coord-----------------------------
def coord2coord(coord,aff):
    """
    This function gives one coordinate and returns the coordinate 
    after applying the affine matrix
    Input: -coord: coordinate
           -aff: Affine matrix(4x4)
    Output -new_coord: new coordinate
    """
    import nibabel as nib
    coord=np.array(coord)
    
    vround = np.vectorize(my_round)
    vint=np.vectorize(int)
    
    new_coord=nib.affines.apply_affine(aff,coord) #coord to new_coord
    new_coord=vint(vround(new_coord))
    
    return(new_coord)


#--------------------- mni2region -----------------------
def mni2region(coord,nii):
    """
    This function gives the index of the region in the atlas
    for the given MNI coordinate 
    Input: - coord: MNI coordinate
           - nii: Atlas as nib object 
    Output: -index of region
    """
    coord=np.array(coord)
    atlas=nii.get_fdata()
    aff=nii.affine
    vox=coord2coord(coord,np.linalg.inv(aff))
    
    r=atlas[vox[0],vox[1],vox[2]]
    return(int(r))


#--------------------- mni2cl_region -----------------------
def mni2cl_region(coord,atlas):
    """
    This function gives the index of the region in the atlas
    for the given MNI coordinate (or the most closests)
    Input: - coord: MNI coordinate
           - atlas: Atlas 
    Output: -index of region
    """
    r=mni2region(coord,atlas)
    if r==0:
        R=0
        neighs=set()
        while True:
            R=R+1
            print(R)
            coord=np.array(coord)
            neigh2={tuple(row) for row in ball_set(coord,np.sqrt(R)) } #set of ball_set
            neigh=neigh2-neighs #minus the last ball_set
            neighs=neigh2
            poss_r=[mni2region(coords,atlas) for coords in neigh]
            poss_r=set(poss_r) - {0} #take away 0 values
            if len(poss_r)>0:
                return (list(poss_r))
                break
    else:
        return([r])


def set_dictionaries_rois(file_info,file_nodes):
    """ 
    Retrieve info for given region (by number)
    Input: - file_info: json file with attributes
           - file_node: The number each attribute for each node
    Output: -r_info function:
             For input i (region i+1) returns MNI coordinates, 
             lobe, Networks 1 and 2 and Broadmann area. 
    """
    with open(file_info, "r") as read_file:
        dict_info = json.load(read_file)
    read_file.close()
    with open(file_nodes, "r") as read_file:
        dict_nodes = json.load(read_file)
    read_file.close()
    
    def r_info(i):
        info=dict_nodes["rois"][i]
    
        l=str(info['attr'][0])
        n1=str(info['attr'][2])
        b=str(info['attr'][3])
        n2=str(info['attr'][4])
        
        result=dict()
        result['MNI']=[info['x'],info['y'],info['z']]
        result['Lobe']= dict_info['gui_Lobes'][l]
        result['Network']=dict_info['gui_Networks'][n1]
        result['Network2']=dict_info['gui_Networks2'][n2]
        result['Brodmann Area']=dict_info['gui_BrodLabels'][b]
        
        return(result)
    return(r_info)


# ## For regions
# <ol> 
#     <li> create limits
#     <li> regions_from_json
#     <li> sort_by_lobes
#     
#     

def create_limits(lengths):
    """
    Create the indices based on lengths
    """
    l=len(lengths)+1
    lim=list(np.zeros(l))
    for i in range(l):
        lim[i]=int(np.sum(lengths[:i]))
    return(lim)


def regions_from_json(file,keys=None,verbose=True):
    """
    Return for each key the list of indices and the total size
    Input:  -file: file with information
            -keys: list of keys
            -verbose: If true print key and size
    """
    with open(file, "r") as read_file:
        info = json.load(read_file)
    lista=[]
    tam=[]
    if keys is None:
        keys=info.keys()
    for key in keys:
        if verbose:
            print(key)
            print(len(info[key]))
        tam=tam+[len(info[key])]
        lista=lista+info[key]
    lista=np.array(lista)-1
    return({'array':lista,'length':tam})


def sort_by_lobes(vector):
    lobes=np.array([22, 33, 37, 50, 71, 82, 99, 119,128,133,
                    157,167,170,184,202,216,235,256,264,268])-1
    lobes2=np.zeros(len(vector))
    for i in range(len(vector)):
        lobe=np.sum(lobes < vector[i]) % 10
        lobes2[i]=lobe
    my_order=np.argsort(lobes2)
    return(my_order)


# ## Dynamics
# <ol>
#     <li> find_alpha_by_density
#     <li> full_correlation 
#     <li> my_avg
#     <li> network_efficiency
#     <li> network_weight

def find_alpha_by_density(array,limit=50,alpha=0,step=0.01):
    total = 100
    while total > 0:
        t,r= array.shape
        density=np.zeros(t)
        for i in range(t):
            density[i]=np.sum(array[i,:]>alpha)*100/r
        total=np.sum(density>limit)
        alpha=alpha+step
        
    return(alpha-step)   


import math
from scipy.stats import zscore
def full_correlation(block,array,inter=True,summary_statistic=None):
    (t,r1,r2) = block.shape
    if t != len(array):
        raise ValueError("array has not same dimensions as block")
    if inter:
        corr=[]
        for i in range(r1):
            for j in range(r2):
                x= np.corrcoef(zscore(block[:,i,j]),zscore(array))[0,1]
                if math.isnan(x):
                    corr.append(0)
                else:
                    corr.append(x)        
    else:
        corr=[]
        for i in range(r1):
            for j in range(i+1,r1):
                x= np.corrcoef(zscore(block[:,i,j]),zscore(array))[0,1]
                if math.isnan(x):
                    corr.append(0)
                else:
                    corr.append(x) 

    if summary_statistic == np.mean:
        corr = np.tanh(summary_statistic(corr))
    elif summary_statistic == np.median:    
        corr = summary_statistic(corr)
    elif not summary_statistic:
        corr
    else:
        raise ValueError("Unrecognized summary_statistic! Use None, np.median, or np.mean.")
    
    return(corr)


def my_avg(x,blocks):
    y=np.copy(x)
    r=blocks.shape[0]
    for i in range(r):
        a=blocks[i,0]-1
        b=blocks[i,1]-1
        for j in np.arange(a,b+1):
            y[j]=np.mean(x[a:(b+1)])
    return(y)


def network_efficiency(block,inter=True,positive=True):
    """
    Measure of strength of connections within and between networks.
    The weight between networks is based on 
    "Dynamics of Intersubject Brain Networks during Anxious Anticipation"
    Najafi M., Kinnison J. and Pessoa L. 2017.
    Input: - block: connections between N1 and N2 networks.
           The dimension of the block (r1,r2) correspond the sizes of the networks
           - alpha: set the threshold for the weights
           - inter: If True the block represents two differents networks
           otherwise is the same network
           -positive: If True only consider positives weights otherwise
           only the negatives are considered
    Output: Strength of connections between networks 
    """
    (r1,r2) = block.shape
    if positive:
        new_block=block*(block>0)
    else:
        new_block=block*(block<0)
    array=np.ndarray.flatten(new_block)
    array=array[array != 0]
    w=np.sum(1/array)
    if inter:
        return(w/(r1*r2))
    else:
        if r1==r2:
            return(w/(r1*(r1-1)))
        else:
            warnings.warn('is not a square matrix for intra weights')


def network_weight(block,alpha=0,inter=True,positive=True):
    """
    Measure of strength of connections within and between networks.
    The weight between networks is based on 
    "Dynamics of Intersubject Brain Networks during Anxious Anticipation"
    Najafi M., Kinnison J. and Pessoa L. 2017.
    Input: - block: connections between N1 and N2 networks.
           The dimension of the block (r1,r2) correspond the sizes of the networks
           - alpha: set the threshold for the weights
           - inter: If True the block represents two differents networks
           otherwise is the same network
           -positive: If True only consider positives weights otherwise
           only the negatives are considered
    Output: Strength of connections between networks 
    """
    (r1,r2) = block.shape
    if positive:
        new_block=block*(block>alpha)
    else:
        new_block=block*(block<-alpha)
    w=np.sum(new_block)
    if inter:
        return(w/(r1*r2))
    else:
        if r1==r2:
            return(w/(r1*(r1-1)))
        else:
            warnings.warn('is not a square matrix for intra weights')


# ##  Distance
# <ol>
#     <li> dyn_corr_isc

from scipy.spatial.distance import squareform
from scipy.stats import pearsonr
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

# ## Filter
# <ol>
# <li> affinity_matrix
# <li> check_connected
# <li> cmdscale
# <li> DFS
# <li> my_isomap

from scipy.sparse import csr_matrix
def affinity_matrix(D,n_neighbors=5,ignore=0,heat_kernel=False,sigma=1,verbose=False):
    """
    Define the affinity martix for a given matrix distance
    Input:  -D: Distance matrix
            -n_neighbors: Number of considered closest neighbors for each node
            -ignore: Number of ignored closest neighbors (n_neighbors is kept)
            -heat_kernel: if False (default) edge weights are set to 1 
            otherwise are set to exp(-D[i,j]**2/sig) 
            -sig : normalize factor for heat kernel
            -verbose: if True print each node ad neighbors
    Output: -aff: Affinitiy matrix (in csr format)
    """
    m=D.shape[0]
    n=D.shape[1]
    if(m != n):
        raise ValueError("D is not a square matrix")
    else:
        aff=np.zeros((m,m))
        for i in range(m):
            diss=D[i,:] #Distance with all 
            idx = np.argsort(-diss)[::-1] #Sort indices ascending
            idx = idx[1+ignore:n_neighbors+1+ignore] #Take n neighbors (ignore itself)
            if (heat_kernel):
                if verbose:
                    print("i=%s, %s" %(i,idx))
                for j in idx:
                    aff[i,j]=np.exp(-D[i,j]**2/sigma)
                    aff[j,i]=aff[i,j]
            else:
                if (verbose == 1):
                    print("i=%s, %s" %(i,idx))
                aff[i,idx]=1 #Set ones to neighbors
                aff[idx,i]=1
    aff=csr_matrix(aff)
    return aff


def check_connected(aff):
    """
    Checks if the affinity matrix aff represents a connected graph
    Input: -aff: Affinitu matrix
    Output: True or False 
    """
    total=aff.shape[0]
    visited=[]
    DFS(aff,visited,0)
    return(total==len(visited))


def cmdscale(D,n_components=2):
    """                                                                                       
    Classical multidimensional scaling (MDS)                                                  
    From https://github.com/meccaLeccaHi/mech_turk/blob/master/analyze/cmdscale.py                                                                                    
    Parameters                                                                                
    ----------                                                                                
    D : (n, n) array                                                                          
        Symmetric distance matrix.                                                            
                                                                                               
    Returns                                                                                   
    -------                                                                                   
    Y : (n, p) array                                                                          
        Configuration matrix. Each column represents a dimension. Only the                    
        p dimensions corresponding to positive eigenvalues of B are returned.                 
        Note that each dimension is only determined up to an overall sign,                    
        corresponding to a reflection.                                                        
                                                                                               
    e : (n,) array                                                                            
        Eigenvalues of B.                                                                     
                                                                                               
    """
    # Number of points                                                                        
    n = len(D)
 
    # Centering matrix                                                                        
    H = np.eye(n) - np.ones((n, n))/n
 
    # YY^T                                                                                    
    B = -H.dot(D**2).dot(H)/2
 
    # Diagonalize                                                                             
    evals, evecs = np.linalg.eigh(B)
 
    #Takes only positive-eigenvalued components
    w, = np.where(evals > 0)
    evals=evals[w]
    evecs=evecs[:,w]
    idx   = np.argsort(evals)[::-1]
    if (n_components<= len(idx)):
        idx=idx[:n_components]
    else:
        warnings.warn('Number of components %s is greater than positive eigenvalues %s' %(n_components,len(idx)))
        
    # Sort by eigenvalue in descending order                                                  
    evals = evals[idx]
    V = evecs[:,idx]
 
    # Compute the coordinates using positive-eigenvalued components only                      
    
    L  = np.diag(np.sqrt(evals))
    Y  = V.dot(L)
 
    return Y, evals


def DFS(aff,visited,node):
    """
    Depth-first search algorithm.
    Input:  -aff: Affinity matrix
            -visited: list to be populated with visited nodes
            -node: initial node
    Returns a populated list of visited nodes from initial node (node)
    """
    if node not in visited:
        #print(node)
        visited.append(node) 
        vecinos=np.where(aff[node,:]!=0)[0]
        #print(vecinos)
        for veci in vecinos:
            DFS(aff,visited,veci)


# +
from sklearn.neighbors import NearestNeighbors, kneighbors_graph
from sklearn.manifold import MDS
from sklearn.utils.graph_shortest_path import graph_shortest_path

def my_isomap (X,n_neighbors=5, n_components=2, eigen_solver='auto', 
               eps=0.001, max_iter=300, path_method='auto',random_states=[None],
               neighbors_algorithm='auto', n_jobs=None, solver="cmds",
               metric='minkowski', p=2, metric_params=None):
        
    nbrs = NearestNeighbors(n_neighbors= n_neighbors,
                                      algorithm=neighbors_algorithm,
                                      metric=metric, p=p,
                                      metric_params=metric_params,
                                      n_jobs=n_jobs)
    
    nbrs.fit(X)
    kng = kneighbors_graph(nbrs, n_neighbors,
                           metric= metric, p=p,
                           metric_params=metric_params,
                           mode='distance', n_jobs=n_jobs)
    kng=(kng + kng.T)/2
    
    dist_matrix = graph_shortest_path(kng,method=path_method,
                                      directed=False)
    print(dist_matrix.shape)
    
    if solver == 'cmds':
        embedding = cmdscale(dist_matrix,n_components=p)[0]
    elif solver == 'mmds':
        st=1e10
        for r_state in random_states:
            mds=MDS(n_components=p,metric=True,dissimilarity='precomputed',max_iter=max_iter,random_state=r_state,eps=eps)
            new_emb = mds.fit_transform(dist_matrix)
            if (mds.stress_< st):
                emb=new_emb
                st=mds.stress_
                print("st=%s, r=%s" %(st,r_state))
        embedding=emb
    else:
        raise Exception("Sorry, solver must be cmds or mmds")
    
    return(embedding)


# -

# ## Test statistic
# <ol>
#     <li> phase_shift
#     <li> surrogate_dyn_corr

#input: 74x193
def phase_shift(X):
    r,t = X.shape
    Y=np.fft.rfft(X)
    amp=np.abs(Y)
    angle = np.angle(Y)

    prng= np.random.RandomState(None)

    random_vec= prng.rand(1,angle.shape[1])*2*np.math.pi 
    angle2 = angle + random_vec
    amp2= amp*(np.cos(angle2)+1j*np.sin(angle2))

    Y2= np.fft.irfft(amp2,n=t) 
    
    return(Y2)


import time
def surrogate_dyn_corr(data, slide, n_shifts=1000, cut=100,path=''):
    """Phase randomization for dyn_corr
    
    For each pair of ROIs compute a null distribution of 
    dyn_corr where time serioes for each subjest are phase 
    randomized prior to computing dyn_corr. 
    
    Input:- data: Matrix of data. Size = (time_points, regions, and n_subjects)
          - slide: Size of window for corr_dyn (size in second= slide*TR)
        - - n_shifts: Number of copies of surrogate data. Default set to 1000
          - cut: Cut proces for cut subjects everytime for generate corr_dyn 
          (to not crash my computer). Default set to 100
          - path: Path of folder where surrogate copies will be stored.
            If not specified 
    At the end we generate n_shifts copies for dyn_corr output
    """
    
    n_TRs = data.shape[0]
    n_regions = data.shape[1]
    n_subjects = data.shape[2]
    
    for j in range(n_shifts):
        surrogate=[]
        for i in range(n_subjects):
            X=data[:,:,i].T
            Y=phase_shift(X)
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
            time.sleep(10)

        sliding_all=np.mean(np.stack(slidings,axis=2),axis=2)
        
        if i<10:
            name=path+'/'+'sliding_all_'+'00'+str(j)+'.npy'
        elif i<100:
            name=path+'/'+'sliding_all_'+'0'+str(j)+'.npy'
        else:
            name=path+'/'+'sliding_all_'+str(j)+'.npy'
            
        np.save(name,sliding_all)


# ## Anchoring
# <ol>
#     <li> bapply
#     <li> compare_2samples
#     <li> critical_kolmogorov

from itertools import combinations_with_replacement,product
#================================================================
#bapply applies a matrix function on blocks of the matrix M
#Input:   -FUN: Function to be applied
#         -M: Matrix
#         -X,Y: List or matrix of intervals (n x 2). A element 
#         of the list or row would be c(1,5) representing the interval between indices 1 and 5
#         -comb: Use all combinations of intervals in X and Y, Otherwise the intervals will be matched
#         one to one
# Output: List of results of FUN applied to defined blocks
def bapply(FUN,M,X,Y=None,comb=False):
    def block_from_index(A,B,index):
        j1=index[0]
        j2=index[1]
        int1=A[j1,:]
        int2=B[j2,:]
        block=M[int1[0]:(int1[1]+1),int2[0]:(int2[1]+1)]
        return(block)
    
    if isinstance(X, list):
        X=np.array(X)
        
    if(X.shape[1]!=2):
        warnings.warn('X must have 2 columns')
    
    if(Y != None):
        if isinstance(Y, list):
            Y=np.array(Y)
        
        if(Y.shape[1]!=2):
            warnings.warn('Y must have 2 columns')
        rY=Y.shape[0]
    
    rX=X.shape[0]
    if (comb): #comb is True
        if (Y == None):
            index=list(combinations_with_replacement(range(rX),2))
            Y=X
        
        else:
            index=list(product(range(rX),range(rX)))
        
    else:
        index=np.reshape([np.arange(rX),np.arange(rX)],(2,rX)).T
        if(Y != None):
            if (rX != rY):
                warnings.warn("X and Y must have same dimensions")
        else:
            Y=X
        
    result=list()
    for i in range(len(index)):
        block=block_from_index(X,Y,index[i])
        result.append(FUN(block))
    
    return(result)   


from scipy.stats import ks_2samp,median_test, kruskal
def compare_2samples(x1,x2,stat_fun,params=None):
    '''
    Compute the stat_fun statisctic on 2 samples.
    Input: -x1, x2: Lists of arrays . The contained arrays must have
            the same dimensions m x n
           -stat_fun: Statistic to compare 1-D samples
           -params: Especifically params list for stat_fun (later)
    Output: -stat: Test statistic of size m x n
            -p: The p-value of test of size m x n
    '''
    n1=len(x1)
    n2=len(x2)
    print("Samples of", n1, "and", n2, "elements")
    dim=np.shape(x1[1])
    
    if len(dim) == 1:
        t=dim[0]
        stat=np.zeros(t)
        p=np.zeros(t)
        
        for i in range(t):
            sample1=np.zeros(n1)
            sample2=np.zeros(n2)
            
            for j in range(n1):
                sample1[j]=x1[j][i]
                
            for j in range(n2):
                sample2[j]=x2[j][i]
    
            result=stat_fun(sample1,sample2)
            stat[i]=result[0]
            p[i]=result[1]
        
    return([stat,p])


def critical_kolmogorov(alpha,n1,n2):
    c_alpha = {
      0.10: 1.22,
      0.05: 1.36,
      0.025: 1.48,
      0.01: 1.63,
      0.005: 1.73,
      0.001: 1.95
    } 
    c=c_alpha[alpha]
    D=c*np.sqrt((n1+n2)/(n1*n2))
    return(D)  


# ## Plot
# <ol>
#     <li> colorline
#     <li> make_segments
#     <li> truncate_colormaps

# +
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.collections as mcoll
import matplotlib.path as mpath

def colorline(
    x, y, z=None, cmap=plt.get_cmap('copper'), norm=plt.Normalize(0.0, 1.0),
        linewidth=3, alpha=1.0,zorder=1):
    """
    http://nbviewer.ipython.org/github/dpsanders/matplotlib-examples/blob/master/colorline.ipynb
    http://matplotlib.org/examples/pylab_examples/multicolored_line.html
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    """

    # Default colors equally spaced on [0,1]:
    if z is None:
        z = np.linspace(0.0, 1.0, len(x))

    # Special case if a single number:
    if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
        z = np.array([z])

    z = np.asarray(z)

    segments = make_segments(x, y)
    lc = mcoll.LineCollection(segments, array=z, cmap=cmap, norm=norm,
                              linewidth=linewidth, alpha=alpha,zorder=zorder)

    ax = plt.gca()
    ax.add_collection(lc)

    return lc


# -

def make_segments(x, y):
    """
    Create list of line segments from x and y coordinates, in the correct format
    for LineCollection: an array of the form numlines x (points per line) x 2 (x
    and y) array
    """

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments


import matplotlib.colors as colors
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


# ### Matlab thing

import sys
import scipy.io
import numpy as np

#----------------------------  find_structure ---------------------------
#This function converts MNI coordinate to a description of brain structure in aal
#Input: - mni : the coordinates (MNI) of some points, in mm.  It is Mx3 matrix
#where each row is the coordinate for one point.
#        -DB (optional): The database. If is omit, make sure TDdatabase.mat is in the 
#        same folder
#Output: -one_line_result: A list of M elements, each describing each point.
#        -table_result:  A  MxN matrix being N the size of the database (DB)
#
def find_structure(mni,DB=None):
    #Vectorize this functions
    vstr=np.vectorize(str)
    vround = np.vectorize(my_round)
    vint=np.vectorize(int)
    
    if DB is None:
        mat=scipy.io.loadmat('../Data/TDdatabase.mat')  
            
    mni=np.array(mni)        
    
    #round coordinates
    mni=vround(mni/2)*2 
    
    T=np.array([[2 ,0 ,0 ,-92],[0,2,0,-128],
                    [0,0,2,-74],[0,0,0,1]])
    
    index=mni2cor(mni,T)
    M=np.shape(index)[0]
    
    #-1 by python indexation
    index=vint(index) - 1 
    
    N=np.shape(mat['DB'])[1]
    table_result=np.zeros((M,N))
    table_result=table_result.tolist() #instead of [i,j] use [i][j]
    
    one_line_result=[""] * M
    
    for i in range(M):
        for j in range(N):
            #mat['DB'][0,j][0,0][0] is the j-th 3D-matrix 
            graylevel=mat['DB'][0,j][0,0][0][index[i,0],index[i,1],index[i,2]] 
            if graylevel == 0:
                 label = 'undefined'
            else:
                if j < (N-1):
                    tmp = ''
                else:
                    tmp =' (aal)' 
                    
                #mat['DB'][0,j][0,0][1]  is the list with regions
                label=mat['DB'][0,j][0,0][1][0,(graylevel-1)][0] + tmp
            
            table_result[i][j]=label
            one_line_result[i] = one_line_result[i] + ' // ' + label
    return(one_line_result,table_result)


import warnings

#---------------------------- mni2cor --------------------------------
def mni2cor(mni,T=np.array([[-4 ,0 ,0 ,84],[0,4,0,-116],
                    [0,0,4,-56],[0,0,0,1]])):
    """
    Convert mni coordinate to matrix coordinate
    Input: - mni : the coordinates (MNI) of some points, in mm.  It is Mx3 matrix
           where each row is the coordinate for one point.
           -T (optional): transform matrix coordinate is the returned coordinate in matrix.
    Output: -coords : Coordinate matrix
    """
    
    mni=np.array(mni)
    
    if len(np.shape(mni))==1:
        mni=mni.reshape((1, len(mni)))
    
    if np.shape(mni)[1] != 3:
        warnings.warn('are not 3-length coordinates')
        return(np.array([]))
        
    a=np.hstack((mni,np.ones((np.shape(mni)[0],1))))
    b=np.transpose(np.linalg.inv(T))
    coords=a.dot(b)
    coords=coords[:,0:3]
        
    vround = np.vectorize(my_round)
    coords = vround(coords)
    return(coords)


import numpy as np
#----------------- my_round -------------------
#Integer rounding that behaves like round from MATLAB
def my_round(x):
    r=x-np.floor(x)
    if(r == 0.5):
        if x<0: return(x-0.5)
        else: return(x+0.5)
    else:
        return(round(x))
