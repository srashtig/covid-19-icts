# masterplot by dates
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime
from matplotlib import gridspec
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
plt.style.use('seaborn-poster')
plt.rcParams['axes.facecolor']='lightgray'
plt.rcParams['savefig.facecolor']='white'

state = ['Karnataka', 'Madhya Pradesh', 'Rajasthan', 'West Bengal', 'Kerala', 'Manipur']


#-------------- DATA ---------------------------------
f = open("refined_df_confirmed.csv", "r")
dates = f.readline().rstrip()
i_dates = dates.split(',')[1:]

t_values = [datetime.datetime.strptime(d, '%d/%m/%y').date() for d in i_dates]
#print t_values
i=0
fig = plt.figure()
ax = plt.subplot(111)
for x in f:
    y = x.rstrip().split(',')
    if y[0] in state:
        N = []
        time = []
        date_time=[]
        print y[0] # <-- country name
        z = [float(a) for a in y[1:]]   # <-- no of affected people in that state from day 0
        ti = z.count(0)+1  # <-- day when first infection appeared
        print dates.split(',')[ti+1-1] # printing the first day of infection
        print dates.split(',')[ti+1-1].replace('/', '-')
        tf = len(z)        # <-- last day in record
        for t in range(0, tf):#ti+1, tf-1):
            if z[t]>=0 and z[t-1]>=0:
                N.append(z[t])
                time.append(t-ti)
                date_time.append(t_values[t])
        ax.semilogy(date_time, N, marker ='o', lw=3, label=y[0])
formatter = mdates.DateFormatter("%d/%m")
ax.xaxis.set_major_formatter(formatter)
locator = mdates.WeekdayLocator()
ax.xaxis.set_major_locator(locator)
locator = mdates.DayLocator()
ax.xaxis.set_minor_locator(AutoMinorLocator(7))#MultipleLocator(7))
ax.tick_params(axis ='both', which ='minor', labelsize=13, length=4, width=2)
ax.tick_params(axis ='both', which ='major', labelsize=13, length=7, width=2)
plt.xlabel('Dates', fontsize=15)
plt.ylabel('Confirmed cases', fontsize=15)
ax.grid(linestyle='--', color='white')
ax.legend(loc='upper center', facecolor='white', fontsize=13, bbox_to_anchor=(0.5, 1.15), fancybox=True, shadow=True, ncol=5)
#ax.xaxis.set_label_position('top')
#ax.xaxis.set_ticks_position('top')
fig.autofmt_xdate(rotation=70)
#plt.tight_layout()
plt.savefig('lockdown/India-master-plot.png')
#plt.show()
f.close()


