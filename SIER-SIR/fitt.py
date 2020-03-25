import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate, optimize


#load data
dat = np.loadtxt('c_confirmed.csv', delimiter=',')

xdata = np.arange(int(len(dat)/4))
ydata = dat[:int(len(dat)/4)]
zs = np.polyfit(xdata, np.log(ydata), 1)
ps = np.poly1d(zs)
print zs
plt.plot(xdata, ps(xdata), '.')
plt.plot(xdata, np.log(ydata))
"""
#inital condition
N = 10**8 #Population size
inf0 = ydata[0]
sus0 = N-inf0
Exp0 = 0
Rec0 = 0

def seir(y, x, beta, sigma, gamma):
    sus = -beta*y[0]*y[2]/N
    expo = (beta*y[0]*y[2]/N)-sigma*y[1]
    infected = sigma*y[1] - gamma*y[2]
    reco = gamma*y[2]
    return sus, expo, infected, reco

def fit_ode(x, beta, sigma, gamma):
    return integrate.odeint(seir, (sus0, Exp0, inf0, Rec0), x, args=(beta, sigma, gamma))[:,2]


popt, pcov = optimize.curve_fit(fit_ode, xdata, ydata)
fitted = fit_ode(xdata, *popt)

plt.plot(xdata, ydata, 'o', mec = 'r')
plt.plot(xdata, fitted)
"""
plt.title("Fit of SEIR model to global infected cases")
plt.ylabel("Population infected")
plt.xlabel("Days")
plt.show()
#print "Optimal parameters:"+str(popt)
