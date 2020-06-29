# number of confirmed, rate, and log rate, indicating lockdown
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime
from matplotlib import gridspec
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
plt.style.use('seaborn-poster')
plt.rcParams['axes.facecolor']='lightgray'
plt.rcParams['savefig.facecolor']='white'

lockdown = {'Karnataka':['2020-01-30', '2020-04-30']
        , 'Rajasthan': ['2020-01-30', '2020-04-30']
        , 'Kerala': ['2020-01-30', '2020-04-30']
        , 'Madhya Pradesh':['2020-01-30', '2020-04-30']
        , 'West Bengal':['2020-01-30', '2020-04-30']
        , 'Tripura':['2020-01-30', '2020-04-30']
        , 'Jammu and Kashmir':['2020-01-30', '2020-04-30'] }
state = lockdown.keys()


#-------------- DATA ---------------------------------
f = open("refined_df_confirmed.csv", "r")
dates = f.readline().rstrip()
i_dates = dates.split(',')[1:]

t_values = [datetime.datetime.strptime(d, '%d/%m/%y').date() for d in i_dates]
#print t_values
i=0
for x in f:
    y = x.rstrip().split(',')
    if y[0] in state:
        if not lockdown[y[0]][1]:
            print "Find lock down date of "+y[0]
        else:
            begining_date = datetime.datetime.strptime(lockdown[y[0]][0], '%Y-%m-%d').date()
            lockdown_date = datetime.datetime.strptime(lockdown[y[0]][1], '%Y-%m-%d').date()
            l_D = lockdown_date.strftime('%m/%d/%y')
        N = []
        dNdt = []
        mu = []
        time = []
        date_time=[]
        print y[0] # <-- country name
        z = [float(a) for a in y[1:]]   # <-- no of affected people in that state from day 0
        ti = z.count(0)+1  # <-- day when first infection appeared
        print dates.split(',')[ti+1-1] # printing the first day of infection
        print dates.split(',')[ti+1-1].replace('/', '-')
        tf = len(z)        # <-- last day in record
        for t in range(ti+1, tf):
            if z[t]>50 and z[t-1]>50:
                N.append(z[t])
                dNdt.append((z[t]-z[t-2])/2.0)
                mu.append((((np.log(z[t])-np.log(z[t-1]))/1.0)+((np.log(z[t])-np.log(z[t-2]))/2.0))/2.0) # left difference
                time.append(t-ti)
                date_time.append(t_values[t])
        fig = plt.figure()
        gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])
        ax = plt.subplot(gs[0])
        ax1 = plt.subplot(gs[1], sharex=ax)
        formatter = mdates.DateFormatter("%d/%m")
        ax.xaxis.set_major_formatter(formatter)
        locator = mdates.WeekdayLocator()
        ax.xaxis.set_major_locator(locator)
        locator = mdates.DayLocator()
        ax.xaxis.set_minor_locator(AutoMinorLocator(7))#MultipleLocator(7))
        ax.semilogy(date_time, N, 'k', marker ='o', lw=3, label='Confirmed cases')
        ax.semilogy(date_time, dNdt, 'r', marker='o', lw=3, label='New Cases / Day')
        ax1.plot(date_time, mu, 'b', marker='o', lw=3, label='Exponential growth rate')
        ax1.set_xlabel('Dates', fontsize=20)
        ax1.set_ylim(0, max(mu)+0.1)
        if not lockdown[y[0]][1]:
            print "plotting error"
            ax.set_title(y[0]+' Start '+str(begining_date)+', Lockdown :', fontsize=20)
        else:
            ax1.axvline(lockdown_date, ls='--', c='orange', lw=3)
            ax.axvline(lockdown_date, ls='--', c='orange', lw=3)
            ax.fill_between(date_time, 1.5*max(N), where=[date_time[i]-lockdown_date>datetime.timedelta(0) for i in range(len(date_time))] , facecolor='gray', interpolate=True)
            ax1.fill_between(date_time, 0, 1, where=[date_time[i]-lockdown_date>datetime.timedelta(-1) for i in range(len(date_time))] ,facecolor='gray', interpolate=True)
            ax.set_title(y[0]+' Start '+str(begining_date)+', Lockdown '+str(lockdown_date), fontsize=20)
        ax.tick_params(axis='both', which='minor', labelsize=13, length=4, width=2)
        ax1.tick_params(axis='both', which='minor', labelsize=13, length=4, width=2)
        ax.tick_params(axis='both', which='major', labelsize=13, length=7, width=2)
        ax1.tick_params(axis='both', which='major', labelsize=13, length=7, width=2)
        ax.grid(linestyle='--', color='white')
        ax1.grid(linestyle='--', color='white')
        ax.legend(loc=0, fontsize=15, shadow=True, facecolor='white')
        ax1.legend(loc=0, fontsize=15, shadow=True, facecolor='white')
        ax.set_ylim(top=2*max(N))
        fig.autofmt_xdate(rotation=70)
        [i.set_linewidth(1.5) for i in ax1.spines.itervalues()]
        plt.tight_layout()
        plt.xlim([min(date_time), max(date_time)])
        plt.subplots_adjust(hspace=.0)
        fig.savefig('lockdown/Ind_'+y[0]+'_lock.png')
        plt.close()
f.close()




