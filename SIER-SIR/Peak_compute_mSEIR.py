#!/usr/bin/env python
# coding: utf-8

# In[48]:


import numpy as np
#from __future__ import print_function     #for python 2
def compute_peak_ak(D0,alpha,sigma,N,gamma,R0,mu):
    
    I0=gamma*D0/(1-alpha)
    Im=sigma/(gamma+sigma) *(1-(1+np.log(R0))/R0)*N
    Tm = 1/mu * np.log(Im/I0)
    PDC = (1-alpha)*gamma*Im
    
    return {'Tm(days)':Tm,'PDC':PDC}


def compute_peak_ad(D0,alpha,sigma,N,gamma,R0,mu):
    gamma_p = gamma
    gamma_e = gamma/(1-alpha/3)
    
    Ip0=D0/gamma_p
    Im=sigma/(gamma_e+sigma) *(1-(1+np.log(R0))/R0)*N
    Ipm = (1-alpha)*gamma_e*Im/gamma_p
    Tm = 1/mu * np.log(Im/Ip0)
    Dpm = Ipm*gamma_p#(1-alpha)*gamma*Im
    
    return {'Tm(days)':Tm,'PDC':Dpm}


# In[51]:


#Delhi: 

## using values from the mail

alpha,sigma,gamma,R0=0.67,1/2,1/5,1.33

D0,N,mu=115.6,1.9e7,0.05

print('AK:',compute_peak_ak(D0,alpha,sigma,N,gamma,R0,mu))
print('AD:',compute_peak_ad(D0,alpha,sigma,N,gamma,R0,mu))


# In[50]:


#India:
## using values from 22may20, Table 1, 1st set.

alpha,sigma,gamma,R0=0.67,1/2,1/5,1.33

D0,N,mu=900,1.3e9,0.05

print('AK:',compute_peak_ak(D0,alpha,sigma,N,gamma,R0,mu))
print('AD:',compute_peak_ad(D0,alpha,sigma,N,gamma,R0,mu))


# In[ ]:





# In[ ]:




