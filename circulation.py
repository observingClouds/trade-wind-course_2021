"""
Functions to calculate circulations
"""

def streamfunction_p(omega, p_levels, firstProfile=None, alpha=1):
    """
    Horizontal integrating the vertical velocity to get the 
    overall streamfunction as in Bretherton et al. (2005) Eq. page 4282)
    but in pressure coordinates:

    psi_i(p)  = psi_(i-1)(p) + omega_{i-1}(p)/g

    omega       :   vertical velocity  [Pa/s]
    p_levels    :   height             [Pa]
    firstProfile:   vertical velocity in first profile (standard: 0)
    """
    import numpy as np
    g        = 9.80665   #Gravitational acceleration [m/s]
    psi      = np.empty((np.shape(omega)[0]+1,np.shape(omega)[1]))
    if firstProfile is None:
        psi[0,:] = 0
    else:
        psi[0,:] = firstProfile
        print(psi[0,:])

    for l in range(len(p_levels)): #per height
        for i in range(1,len(omega)+1): #per rank
            psi[i,l] = psi[i-1,l]+omega[i-1,l]*alpha
    psi = psi/g # to convert to kg/m**2/s

    return psi

def wtg(T,Tpot,Qcool_day,dTpot_dp):
    """ Calculate vertical velocity by using 
    the weak temperature gradient approximation 

    Input:
    Tpot:  potential temperature
    Qcool: cooling [K/d]

    Output:
    vertical velocity in Pa/s

    omega = Qcool/S where S= T/Tpot * dTpot/dp
    """
    Qcool_sec = Qcool_day/(24*60*60.)
    return Qcool_sec/(T/Tpot*(dTpot_dp))

