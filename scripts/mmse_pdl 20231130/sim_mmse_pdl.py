from gnpy.core.parameters import Parameters
from gnpy.core.utils import generate_random_numbers
import numpy as np
import matplotlib.pyplot as plt
from gnpy.core.elements import Fiber, Edfa, Roadm
from gnpy.core.mmse_pdl import pdl_pen
from support.multi_processing import MyProcess
import time

def sim_pdl_pen(*args):
    params = Parameters()
    params.snr_trx_dB = 17
    params.sig_power_dBm = 0
    params.Rs = 100e9
    params.Fc = 193e12

    link_params = Parameters()
    link_config = ["wss", "oa", "voa", "fiu", "fiber", "oa",  "voa", "fiu", "fiber", "oa", "voa", "fiu",  "fiber", "oa", "voa", "fiu",  "fiber", "wss"]*3
    wss_params_list = []
    oa_params_list = []
    fiber_params_list = []
    voa_params_list = []
    fiu_params_list = []

    for comp in link_config:
        if comp.lower() == 'wss':
            wss_params = Parameters()
            wss_params.loss_dB = 6
            wss_params.pdl_dB = 0.5
            wss_params_list.append(wss_params)
        elif comp.lower() == 'oa':
            oa_params = Parameters()
            oa_params.nf_dB = 5
            oa_params.gain_dB = 18           
            oa_params.pdl_dB  = 0.3
            oa_params_list.append(oa_params)
        elif comp.lower() == 'fiber':
            fiber_params = Parameters()
            fiber_params.loss_dB = 16
            fiber_params_list.append(fiber_params)
        elif comp.lower() == 'voa':
            voa_params = Parameters()
            voa_params.loss_dB = 1
            voa_params.pdl_dB  = 0.1
            voa_params_list.append(voa_params)
        elif comp.lower() == 'fiu':
            fiu_params = Parameters()
            fiu_params.loss_dB = 1
            fiu_params.pdl_dB  = 0.1
            fiu_params_list.append(fiu_params)

    link_params.link_config = link_config
    link_params.num_oa = link_config.count("oa")
    link_params.num_wss = link_config.count("wss")
    link_params.num_fiber = link_config.count("fiber")
    link_params.num_fiu = link_config.count("fiu")
    link_params.num_voa = link_config.count("voa")
    params.link_params = link_params
    params.oa_params_list = oa_params_list[::-1]
    params.wss_params_list = wss_params_list[::-1]
    params.fiber_params_list = fiber_params_list[::-1]
    params.voa_params_list = voa_params_list[::-1]
    params.fiu_params_list = fiu_params_list[::-1]
    params.alpha_list = generate_random_numbers('pdf',(lambda x:np.sin(2*x),0,np.pi/2),len(link_config))
    params.beta_list = generate_random_numbers('uniform',(0,np.pi*2),len(link_config))
    snr_pen_dB = pdl_pen(params)

    return snr_pen_dB

if __name__ == '__main__':
    N_sim = 10000
    # 多进程
    args_list = [0]*N_sim
    mp_start_time = time.time()
    mp = MyProcess(func=sim_pdl_pen,args_list=args_list)
    snr_pen_dB = mp.run_multiProcess()

    mp_end_time = time.time()
    print("多进程耗时: {:.2f}秒".format(mp_end_time - mp_start_time))
    hist, bin_edges = np.histogram(snr_pen_dB, bins=30, density=True)   
    plt.bar(bin_edges[:-1], 10*np.log10(hist), width=np.diff(bin_edges), edgecolor="k", alpha=0.7,label='simulation')
    plt.xlabel('PDL')
    plt.ylabel('Frequency')
    plt.title('Histogram')
    plt.legend()
    plt.show()