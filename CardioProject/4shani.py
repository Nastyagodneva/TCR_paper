import numpy as np
import scipy
from scipy.cluster.hierarchy import linkage, leaves_list, cut_tree


def get_var_order(sim):
    """
    generate var order based on hierarchical clustering of the given similarity matrix
    """
    
    dist_condensed = -(sim[np.triu_indices(sim.shape[0], 1)]-1)
    dist_condensed = np.clip(dist_condensed,0,2,dist_condensed)
    lnk = linkage(dist_condensed, method='average', metric='euclidean')
    leaf_order = leaves_list(lnk)
    return leaf_order


# say tmp is some distance matrix
tmporder = get_var_order(-tmp) # take - bc get_var_order expects similarity matrices
tmp = tmp[np.ix_(tmporder, tmporder)]
fig, axes = plt.subplots(figsize=(10,10))
axes.imshow(tmp, interpolation=None)
