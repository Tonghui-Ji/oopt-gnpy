import numpy as np
from gnpy.core.utils import lin2db,db2lin,snr2ber,ber2snr
from scipy.constants import pi,h,c
from gnpy.core.parameters import Parameters
from scipy.special import erfc
from scipy.stats import unitary_group
import matplotlib.pyplot as plt

def h_mc(D:int=4)->np.array:
    """
    generate random unitary matrix

    Parameters
    ----------
    D : int scalar
        number of modes.
    Returns
    -------
    h_rsop : complex-valued ndarray
        mode coupling matrix.

    References
    -------
    [1] https://arxiv.org/abs/math-ph/0609050  

    """
    unitary_group.rvs(D)
    


def h_mdl(D:int=4, mdl_pp_dB:float=1, type:str='active')->np.array:
    """
    generate mdl transfer matrix
    While the peak-to-peak MDG is relevant for predicting the transmission performance of weakly coupled MDM links.
    For a single mdl element, p-p and rms are both required to describle the property.
    Denote the distrubition of mdl value is uniform.
    A fixed variance uniform distribution is not easy to genearte, so we generate h_mdl using p-p and return rms value.

    Parameters
    ----------
    D : int scalar
        number of modes.
    mdl_pp_dB : float scalar
        p-p mdl value in dB.
    mdl_rms_dB : float scalar
        rms mdl value in dB.
    alpha : float scalar
        rotation angle of transfer matrix.
    beta : float scalar
        phase difference between two pols.
    type : str
        type of pdl model
    Returns
    -------
    h_mdl : complex-valued ndarray
        transfer matrix.

    References
    -------
    [1] Impact of Polarization- and Mode-Dependent Gain on the Capacity of Ultra-Long-Haul Systems.    
    """
    
    # TODO: In some papers, mdl_db = ln(mdl_lin), which I think is not consistent with the definitio of pdl. 
    #       So here I temporarily choose mdl_db = 10*log10(mdl_lin).

    
    if type.lower() == 'active':
        mdl_tmp_dB = np.sort(np.random.uniform(0, mdl_pp_dB/2, int((D-2)/2)))
        mdl_array_dB = np.concatenate((
            np.array([-mdl_pp_dB/2]), 
            -mdl_tmp_dB[::-1],
            mdl_tmp_dB,
            np.array([mdl_pp_dB/2])
        ))
        h_mdl = np.diag(db2lin(mdl_array_dB))
        mdl_rms_dB = np.std(mdl_array_dB)
        pass
    elif type.lower() == 'passive':
        pass

    else:
        raise TypeError('Unsuporrted pdl type')
    
    h_mc = unitary_group.rvs(D)
    h_mdl = np.dot(h_mc, h_mdl)
    h_mc = unitary_group.rvs(D)
    h_mdl = np.dot(h_mdl,h_mc)

    return h_mdl, mdl_rms_dB

def prod_mat(mat_list:list)->np.array:
    """
    product of matrix 

    Parameters
    ----------
    mat_list : array list
        matrix list to be producted.
    Returns
    -------
    mat_prod : complex-valued array
        product result.

    References
    -------    
    """    
    mat_prod = np.identity(mat_list[0].shape[0])

    for tmp_mat in mat_list:
        mat_prod = np.dot(tmp_mat,mat_prod)

    return mat_prod


if __name__ == '__main__': 
    N = 10000
    D = 2
    mdl_tot_pp_dB = np.zeros(N)
    mdl_tot_rms_coupled_dB = np.zeros(N)
    for j in range(N):
        mdl_tot_var_uncoupled_dB = 0
        h_mdl_list = []
        for i in range(1000):
            h_mdl_tmp,mdl_tmp_rms_dB = h_mdl(D=D, mdl_pp_dB=1)
            mdl_tot_var_uncoupled_dB += mdl_tmp_rms_dB**2 
            h_mdl_list.append(h_mdl_tmp)
        mdl_tot_rms_uncoupled_dB = np.sqrt(mdl_tot_var_uncoupled_dB)
        h_mdl_tot = prod_mat(h_mdl_list)
        U,S,V = np.linalg.svd(h_mdl_tot)
        mdl_tot_pp_dB[j] = 2*lin2db(max(S)/min(S))
        mdl_tot_rms_coupled_dB[j] = np.std(2*lin2db(S))
        mdl_tot_rms_coupled_esti_dB = mdl_tot_rms_uncoupled_dB*np.sqrt(1+mdl_tot_var_uncoupled_dB/(12*(1-1/D**2)))

        pass

    plt.hist(mdl_tot_pp_dB, bins=10, edgecolor='black')
    plt.show()