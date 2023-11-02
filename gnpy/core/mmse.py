import numpy as np
import scipy.special as sp
from scipy.constants import pi,h,c
import matplotlib.pyplot as plt
from support.utils import lin2dB,dB2lin,parameters

def h_wss_f(f_axis,wss_params):
    def super_guassian(f_axis,bw,n,m):
        sigma_sg = bw/(2*np.sqrt( 2*np.log(np.sqrt(10**(m/10)))**(1/n) ))
        h_sg_f = 1/sigma_sg/np.sqrt(2*np.pi)*np.exp(-(f_axis**2/2/sigma_sg**2)**n)
        h_sg_f = h_sg_f/np.max(np.abs(h_sg_f))

        return h_sg_f


    def lcos_model(f_axis,bw,bw_otf):
        sigma = bw_otf/2/np.sqrt(2*np.log(2))
        h_lcos_f = 0.5*sigma*np.sqrt(2*np.pi)*(sp.erf((bw/2-f_axis)/(np.sqrt(2)*sigma))-sp.erf((-bw/2-f_axis)/(np.sqrt(2)*sigma)))

        h_lcos_f = h_lcos_f/np.max(np.abs(h_lcos_f))

        return h_lcos_f
    
    if wss_params is None:
        wss_params = []

    # check input parameters
    type = getattr(wss_params, "type", "super_guassian")
    bw = getattr(wss_params, "bw", 0)
    if type == 'super_guassian':
        m = getattr(wss_params, "m", 6)
        n = getattr(wss_params, "n", 6)
        h_wss_f = super_guassian(f_axis,bw,n,m)
    elif type == 'lcos_model':
        bw_otf = getattr(wss_params, "bw_otf", 10e9)
        h_wss_f = lcos_model(f_axis,bw,bw_otf)
    elif type == 'ideal':
        h_wss_f = np.ones_like(f_axis)

    return h_wss_f

def rc_filter(N,alpha,Rs,Fs):
    f_axis = np.linspace(-Fs/2,Fs/2-Fs/N,N)
    Ts = 1/Rs
    h_f = Ts*np.cos(np.pi*Ts/2/alpha*(np.abs(f_axis)-(1-alpha)/(2*Ts)))
    h_f[np.where(np.abs(f_axis)<=(1-alpha)/(2*Ts))] = Ts
    h_f[np.where(np.abs(f_axis)>(1+alpha)/(2*Ts))] = 0

    return h_f


def wss_pen(params):
    N = params.N
    Rs = params.Rs
    Fs = params.Fs
    Fc = params.Fc
    link_config = params.link_params.link_config
    roll_off = params.roll_off
    snr_trx_dB = params.snr_trx_dB

    oa_params_list = params.oa_params_list
    wss_params_list = params.wss_params_list
    fiu_params_list = params.fiu_params_list
    fiber_params_list = params.fiber_params_list
    voa_params_list = params.voa_params_list

    sig_power_dB = params.sig_power_dB
    sig_power_lin = dB2lin(sig_power_dB)

    f_axis = np.linspace(-Fs/2,Fs/2-Fs/N,N)

    h_rc_f = rc_filter(N,roll_off,Rs,Fs)


    h_f = np.ones((N,len(link_config)+2))
    n_power_lin = np.zeros(len(link_config)+2)

    n_power_lin[0] = sig_power_lin / dB2lin(snr_trx_dB) / 2


    for i,comp in enumerate(link_config):
        if comp == 'wss':
            wss_params = wss_params_list.pop()
            loss_lin = dB2lin(wss_params.loss_dB)
            sig_power_lin = sig_power_lin / loss_lin
            n_power_lin = n_power_lin / loss_lin
            h_f[:,i] = h_wss_f(f_axis,wss_params)
            pass
        if comp == 'oa':
            oa_params = oa_params_list.pop()
            gain_lin = dB2lin(oa_params.gain_dB)
            nf_lin = dB2lin(oa_params.nf_dB)
            sig_power_lin = sig_power_lin * gain_lin
            n_power_lin = n_power_lin * gain_lin
            n_power_lin[i+1] = h*Fc*Rs*(gain_lin*nf_lin-1)
            pass
        if comp == 'fiber':
            fiber_params = fiber_params_list.pop()
            loss_lin = dB2lin(fiber_params.loss_dB)
            sig_power_lin = sig_power_lin / loss_lin
            n_power_lin = n_power_lin / loss_lin        
            pass
        if comp =='fiu':
            fiu_params = fiu_params_list.pop()
            loss_lin = dB2lin(fiu_params.loss_dB)
            sig_power_lin = sig_power_lin / loss_lin
            n_power_lin = n_power_lin / loss_lin  
            pass
        if comp == 'voa':
            voa_params = voa_params_list.pop()
            loss_lin = dB2lin(voa_params.loss_dB)
            sig_power_lin = sig_power_lin / loss_lin
            n_power_lin = n_power_lin / loss_lin  
            pass
    n_power_lin[i+1] = sig_power_lin / dB2lin(snr_trx_dB) / 2

    C_f = np.prod(h_f,axis=1)
    S_f = np.zeros(N) 
    alpha = n_power_lin/sum(n_power_lin)

    for i in range(len(n_power_lin)):
        T_f = np.prod(h_f[:,i:],axis=1)
        S_f += alpha[i]*T_f

    W_f = 1/np.sqrt(S_f)
    H_eq_f = W_f*C_f
    fig,ax = plt.subplots(nrows=2,ncols=2)
    ax[0,0].plot(f_axis,np.abs(C_f))
    ax[0,0].set_title('C_f')
    ax[0,1].plot(f_axis,np.abs(S_f))
    ax[0,1].set_title('S_f')
    ax[1,0].plot(f_axis,np.abs(W_f))
    ax[1,0].set_title('W_f')
    ax[1,1].plot(f_axis,np.abs(H_eq_f))
    ax[1,1].set_title('H_eq_f')    
    plt.show()

    snr_esti_dB = lin2dB(sig_power_lin*np.mean(H_eq_f[int(N/2):int(3*N/2)])/sum(n_power_lin))
    pass



if __name__ == '__main__':
    params = parameters()
    params.N = 1024
    params.Rs = 92e9
    params.Fs = 2*params.Rs
    params.Fc = 193e12
    
    params.roll_off = 0.1 
    params.snr_trx_dB = 17
    params.sig_power_dB = 3

    link_params = parameters()
    link_config = ["wss","oa","voa","fiber","fiu","oa"]
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
        wss_params = parameters()
        wss_params.loss_dB = 7
        wss_params.type = 'lcos_model'
        wss_params.bw_otf = 10e9
        wss_params.bw = 100e9
        wss_params_list.append(wss_params)

    for i in range(link_params.num_oa):
        oa_params = parameters()
        oa_params.nf_dB = 4.5
        oa_params.gain_dB = 15
        oa_params_list.append(oa_params)

    for i in range(link_params.num_fiber):
        fiber_params = parameters()
        fiber_params.loss_dB = 16
        fiber_params_list.append(fiber_params)

    for i in range(link_params.num_fiu):
        fiu_params = parameters()
        fiu_params.loss_dB = 0.7
        fiu_params_list.append(fiu_params)

    for i in range(link_params.num_voa):
        voa_params = parameters()
        voa_params.loss_dB = 1
        voa_params_list.append(voa_params)

    params.oa_params_list = oa_params_list[::-1]
    params.wss_params_list = wss_params_list[::-1]
    params.fiu_params_list = fiu_params_list[::-1]
    params.fiber_params_list = fiber_params_list[::-1]
    params.voa_params_list = voa_params_list[::-1]

    wss_pen(params)


    pass