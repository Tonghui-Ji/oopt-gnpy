import numpy as np
from gnpy.core.utils import lin2db,db2lin,snr2ber,ber2snr
from scipy.constants import pi,h,c
from gnpy.core.parameters import Parameters
from scipy.special import erfc

def h_pdl(pdl_dB:float=0,alpha:float=0,beta:float=0,type:str='active')->np.array:
    """
    generate pdl transfer matrix

    Parameters
    ----------
    pdl_dB : float scalar
        pdl value in dB.
    alpha : float scalar
        rotation angle of transfer matrix.
    beta : float scalar
        phase difference between two pols.
    type : str
        type of pdl model
    Returns
    -------
    h_pdl : complex-valued ndarray
        transfer matrix.

    References
    -------
    [1] PDL in Optical Links: A Model Analysis and a Demonstration of a PDL-Resilient Modulation.    
    """
    def h_rsop(alpha,beta):
        h_rsop = np.array([[np.cos(alpha)*np.exp(-1j*beta),-np.sin(alpha)],[np.sin(alpha),np.cos(alpha)*np.exp(1j*beta)]])

        return h_rsop
    

    pdl_lin = 10**(pdl_dB*0.1)
    if type.lower() == 'active':
        gamma = (1 - pdl_lin) / (1 + pdl_lin)
        h_pdl = np.array([[1+gamma,0],[0,1-gamma]])
    elif type.lower() == 'passive':
        epsilon = 1 / pdl_lin
        h_pdl = np.array([[1,0],[0,np.sqrt(epsilon)]])
    else:
        raise TypeError('Unsuporrted pdl type')
    
    h_pdl = np.dot(h_rsop(-alpha,-beta),np.dot(h_pdl,h_rsop(alpha,beta)))


    return h_pdl

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
    params = Parameters()
    params.snr_trx_dB = 17
    params.sig_power_dB = 3
    params.Rs = 92e9
    params.Fc = 193e12

    link_params = Parameters()
    link_config = ["wss","wss","oa","voa","fiber","fiu","oa"]
    link_params.link_config = link_config
    link_params.num_oa = link_config.count("oa")
    link_params.num_wss = link_config.count("wss")
    link_params.num_fiber = link_config.count("fiber")
    link_params.num_fiu = link_config.count("fiu")
    link_params.num_voa = link_config.count("voa")
    params.link_params = link_params

    wss_params_list = []
    oa_params_list = []
    fiu_params_list = []
    fiber_params_list = []
    voa_params_list = []
    for i in range(link_params.num_wss):
        wss_params = Parameters()
        wss_params.loss_dB = 8
        wss_params.type = 'lcos_model'
        wss_params.bw_otf = 10e9
        wss_params.bw = 100e9
        wss_params.pdl_dB = 0.55
        wss_params_list.append(wss_params)

    for i in range(link_params.num_oa):
        oa_params = Parameters()
        oa_params.nf_dB = 4.5
        oa_params.gain_dB = 15
        oa_params.pdl_dB = 0.1
        oa_params_list.append(oa_params)

    for i in range(link_params.num_fiber):
        fiber_params = Parameters()
        fiber_params.loss_dB = 16
        fiber_params_list.append(fiber_params)

    for i in range(link_params.num_fiu):
        fiu_params = Parameters()
        fiu_params.loss_dB = 0.7
        fiu_params.pdl_dB = 0.1
        fiu_params_list.append(fiu_params)

    for i in range(link_params.num_voa):
        voa_params = Parameters()
        voa_params.loss_dB = 1
        voa_params.pdl_dB = 0.1
        voa_params_list.append(voa_params)

    params.oa_params_list = oa_params_list[::-1]
    params.wss_params_list = wss_params_list[::-1]
    params.fiu_params_list = fiu_params_list[::-1]
    params.fiber_params_list = fiber_params_list[::-1]
    params.voa_params_list = voa_params_list[::-1]

    pdl_pen(params)