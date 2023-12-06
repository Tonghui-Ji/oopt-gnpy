'''
This file is based on elements, to create elements through simple parameters.
2023.12.06: 不熟悉GNPY中各个器件的实际工作原理，很难完成这部分工作，暂时先放着，该文件咱不可用，后面根据需要再写。
'''
from gnpy.core import elements

def fiber_from_params(uid='fiber_0', type='ssmf', length=80, length_unit='km',delta_p=0):
    el_config = {}
    el_config['uid'] = uid
    el_config['metadata'] = {}
    el_config['params'] = {'length':length,'length_units':length_unit,'con_in':None,'con_out':None} 
    el_config['operatioal'] = {}

    if type.lower() == 'ssmf':
        el_config['params']['loss_coef'] = 0.17
        el_config['params']['pmd_coef'] = 1.265e-12
    else:
        raise Warning('Unsupported fiber type')


    return elements.Fiber(**el_config)

def edfa_from_params(uid='edfa_0',type='advanced_mode', mode='apc', gain_target=16,delta_p=0):
    el_config = {}
    el_config['uid'] = uid
    el_config['metadata'] = {}
    el_config['operational'] = {'gain_target':18,'tilt_target':0}
    el_config['params'] = {
        'f_min': 191.35e12,
        'f_max': 196.1e12,
        'type_variety': '',
        'type_def': type,
        'gain_flatmax': None,
        'gain_min': None,
        'p_max': None,
        'nf_model': None,
        'dual_stage_model': None,
        'preamp_variety': None,
        'booster_variety': None,
        'nf_min': None,
        'nf_max': None,
        'nf_coef': None,
        'nf0': None,
        'nf_fit_coeff': None,
        'nf_ripple': 0,
        'dgt': None,
        'gain_ripple': 0,
        'tilt_ripple': 0,
        'f_ripple_ref': None,
        'out_voa_auto': False,
        'allowed_for_design': False,
        'raman': False,
        'pmd': 0,
        'pdl': 0,
        'advance_configurations_from_json': None
    }

    if mode.lower() == 'apc':
        el_config['operational']['delta_p'] = delta_p
    elif mode.lower() == 'agc':
        el_config['operational']['gain_target'] = gain_target
    else:
        raise Warning('Unsupported edfa mode')

    return elements.Fiber(**el_config)


if __name__ == '__main__':
    fiber = fiber_from_params(uid='fiber_0', type='ssmf', length=80, length_unit='km',delta_p=0)

    pass
