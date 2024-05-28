import numpy as np
from gnpy.core.utils import lin2db,db2lin,snr2ber,ber2snr
from scipy.constants import pi,h,c
from gnpy.core.parameters import Parameters
from scipy.special import erfc

def h_pdl(D:int=4, mdl_pp_dB:float=0, mdl_rms_dB:float=0, type:str='active')->np.array:
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
    mat_prod = np.array([[1,0],[0,1]])

    for tmp_mat in mat_list:
        mat_prod = np.dot(tmp_mat,mat_prod)

    return mat_prod

def pdl_pen(params:Parameters)->float:
    """
    Calculate snr pen induced by pdl based on mmse criterion.

    Parameters
    ----------
    params : Parameters object
        Simulation Parameters.
    Returns
    -------
    snr_pen : real scalar
        Snr penatly induce by pdl.

    References
    -------
    [1] Analysis of Impact of Polarization Dependent Loss in Point to Multi-Point Subsea Communication Systems.    
    [2] QAM BER for AWGN channel: https://www.etti.unibw.de/labalive/experiment/qam/
    """
    link_config = params.link_params.link_config
    snr_trx_dB = params.snr_trx_dB
    Rs = params.Rs
    Fc = params.Fc

    oa_params_list = getattr(params,"oa_params_list",None)
    wss_params_list = getattr(params,"wss_params_list",None)
    fiu_params_list = getattr(params,"fiu_params_list",None)
    fiber_params_list = getattr(params,"fiber_params_list",None)
    voa_params_list = getattr(params,"voa_params_list",None)
    alpha_list = getattr(params,"alpha_list",None)
    beta_list = getattr(params,"beta_list",None)

    sig_power_dBm = params.sig_power_dBm
    sig_power_lin = db2lin(sig_power_dBm)/1e3

    n_power_lin = np.zeros(len(link_config)+2)

    n_power_lin[0] = sig_power_lin / db2lin(snr_trx_dB) / 2
    h_pdl_list = []


    for i,comp in enumerate(link_config):
        if comp.lower() == 'wss':
            wss_params = wss_params_list.pop()
            loss_lin = db2lin(getattr(wss_params,"loss_dB",8))
            pdl_dB = getattr(wss_params,"pdl_dB",0.5)
            sig_power_lin = sig_power_lin / loss_lin
            n_power_lin = n_power_lin / loss_lin
            h_pdl_list.append(h_pdl(pdl_dB=pdl_dB,alpha=alpha_list[i],beta=beta_list[i]))
            pass
        elif comp.lower() == 'oa':
            oa_params = oa_params_list.pop()
            gain_lin = db2lin(getattr(oa_params,"gain_dB",16))
            nf_lin = db2lin(getattr(oa_params,"nf_dB",4.5))
            pdl_dB = getattr(oa_params,"pdl_dB",0.3)
            sig_power_lin = sig_power_lin * gain_lin
            n_power_lin = n_power_lin * gain_lin
            n_power_lin[i+1] = h*Fc*Rs*(gain_lin*nf_lin-1)
            h_pdl_list.append(h_pdl(pdl_dB=pdl_dB,alpha=alpha_list[i],beta=beta_list[i]))
            pass
        elif comp.lower() == 'fiber':
            fiber_params = fiber_params_list.pop()
            loss_lin = db2lin(getattr(fiber_params,"loss_dB",16))
            pdl_dB = getattr(fiber_params,"pdl_dB",0)
            sig_power_lin = sig_power_lin / loss_lin
            n_power_lin = n_power_lin / loss_lin     
            h_pdl_list.append(h_pdl(pdl_dB=pdl_dB,alpha=alpha_list[i],beta=beta_list[i]))   
            pass
        elif comp.lower() =='fiu':
            fiu_params = fiu_params_list.pop()
            loss_lin = db2lin(getattr(fiu_params,"loss_dB",1))
            pdl_dB = getattr(fiu_params,"pdl_dB",0.1)
            sig_power_lin = sig_power_lin / loss_lin
            n_power_lin = n_power_lin / loss_lin  
            h_pdl_list.append(h_pdl(pdl_dB=pdl_dB,alpha=alpha_list[i],beta=beta_list[i])) 
            pass
        elif comp.lower() == 'voa':
            voa_params = voa_params_list.pop()
            loss_lin = db2lin(getattr(voa_params,"loss_dB",1))
            pdl_dB = getattr(voa_params,"pdl_dB",0.1)
            sig_power_lin = sig_power_lin / loss_lin
            n_power_lin = n_power_lin / loss_lin  
            h_pdl_list.append(h_pdl(pdl_dB=pdl_dB,alpha=alpha_list[i],beta=beta_list[i])) 
        else:
            raise TypeError('Unsupported component type')
    n_power_lin[-1] = sig_power_lin / db2lin(snr_trx_dB) / 2    
    h_pdl_list.append(h_pdl(pdl_dB=0,alpha=0,beta=0)) # 和收端ASE匹配，便于计算

    H_ch = prod_mat(h_pdl_list)
    alpha = n_power_lin/sum(n_power_lin)
    
    S = np.mat([[0,0],[0,0]],dtype='complex')
    for i in range(len(n_power_lin)):
        # Note: T and S is matrix, just for calculation convenient
        T = np.mat(prod_mat(h_pdl_list[i:]))
        S += alpha[i]*(T*T.H)

    W = S.I
    H_eq = W*H_ch

    snr_ase_lin = sig_power_lin/sum(n_power_lin)
    MMSE = np.abs((H_eq*H_eq.H + 1/snr_ase_lin*np.eye(2)).I)
    snr_esti_dB_x = lin2db(sig_power_lin/MMSE[0,0]/sum(n_power_lin)-1)
    snr_esti_dB_y = lin2db(sig_power_lin/MMSE[1,1]/sum(n_power_lin)-1)

    # convert SNR to ber, calculate average ber, conver average ber to average snr
    M = 16
    ber_x = snr2ber(M,snr_esti_dB_x)
    ber_y = snr2ber(M,snr_esti_dB_y)
    snr_pdl_dB = ber2snr(M,(ber_x + ber_y) / 2)

    snr_ase_dB = lin2db(sig_power_lin/sum(n_power_lin))
    snr_pen_dB = snr_ase_dB - snr_pdl_dB

    return snr_pen_dB


if __name__ == '__main__':
    h_pdl(D=4, mdl_pp_dB=3, mdl_rms_dB=1)