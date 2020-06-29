# number of confirmed, rate, and log rate, indicating lockdown
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime
from matplotlib import gridspec
plt.style.use('seaborn-poster')
plt.rcParams['axes.facecolor']='lightgray'
plt.rcParams['savefig.facecolor']='white'

lockdown = {'Italy':['2020-01-31', '2020-03-08']
        , 'India': ['2020-01-30', '2020-03-24']
        , 'Korea South': ['2020-01-22', '']
        , 'Germany':['2020-01-27', '2020-03-22']
        , 'US':['2020-01-22', '']
        , 'United Kingdom':['2020-01-31', '2020-03-23']
        , 'Iran':['2020-02-19', '2020-03-14']
        , 'China':['2020-01-22', '2020-01-23']
        , 'France':['2020-01-24', '2020-03-16']
        , 'Spain':['2020-02-01', ''] }
#country = ['Italy', 'India', 'Korea South', 'Germany', 'US', 'United Kingdom', 'Iran', 'China', 'France', 'Spain']
country = lockdown.keys()
province = ['', 'Hubei']


#-------------- DATA ---------------------------------
f = open("time_series_global-confirmed_edit.csv", "r")
dates = f.readline()
i_dates = dates.split(',')[4:-2]

t_values = [datetime.datetime.strptime(d, '%m/%d/%y').date() for d in i_dates]
#print t_values
i=0
for x in f:
    y = x[:-1].split(',')
    if y[0] in province and y[1] in country:
        if not lockdown[y[1]][1]:
            print "Find lock down date of "+y[1]
        else:
            begining_date = datetime.datetime.strptime(lockdown[y[1]][0], '%Y-%m-%d').date()
            lockdown_date = datetime.datetime.strptime(lockdown[y[1]][1], '%Y-%m-%d').date()
            l_D = lockdown_date.strftime('%m/%d/%y')
        N = []
        dNdt = []
        mu = []
        time = []
        date_time=[]
        print y[1] # <-- country name
        z = [float(a) for a in y[4:]]   # <-- no of affected people in that state from day 0
        ti = z.count(0)+1  # <-- day when first infection appeared
        print dates.split(',')[ti+4-1] # printing the first day of infection
        print dates.split(',')[ti+4-1].replace('/', '-')
        tf = len(z)        # <-- last day in record
        for t in range(ti+1, tf-3):
            if z[t]>20 and z[t-1]>20:
                N.append(z[t])
                dNdt.append((((z[t+1]-z[t-1])/2.0)+((z[t]-z[t-1]))+((z[t+2]-z[t-1])/3.0))/3.0)
                mu.append(( np.log(z[t+1])-np.log(z[t-1]) )/2.0) # central difference
                time.append(t-ti)
                date_time.append(t_values[t])
        fig = plt.figure()
        gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
        ax = plt.subplot(gs[0])
        ax1 = plt.subplot(gs[1], sharex=ax)
        formatter = mdates.DateFormatter("%m/%d")
        ax.xaxis.set_major_formatter(formatter)
        locator = mdates.WeekdayLocator()
        ax.xaxis.set_major_locator(locator)
        locator = mdates.DayLocator()
        ax.xaxis.set_minor_locator(locator)
        ax.semilogy(date_time, N, 'k', marker ='o', lw=3, label='Confirmed cases')
        ax.semilogy(date_time, dNdt, 'r', marker='o', lw=3, label='New Cases / Day')
        ax1.plot(date_time, mu, 'b', marker='o', lw=3, label='Change of exponential growth rate')
        ax1.set_xlabel('Dates', fontsize=20)
        if not lockdown[y[1]][1]:
            print "plotting error"
            ax.set_title(y[1]+':'+y[0]+' Start '+str(begining_date)+', Lockdown :', fontsize=20)
        else:
            ax1.axvline(lockdown_date, ls='--', c='orange', lw=3)
            ax.axvline(lockdown_date, ls='--', c='orange', lw=3)
            ax.fill_between(date_time, 1.5*max(N), where=[date_time[i]-lockdown_date>datetime.timedelta(0) for i in range(len(date_time))] , facecolor='gray', interpolate=True)
            ax1.fill_between(date_time, 0, 1, where=[date_time[i]-lockdown_date>datetime.timedelta(-1) for i in range(len(date_time))] ,facecolor='gray', interpolate=True)
            ax.set_title(y[1]+':'+y[0]+' Start '+str(begining_date)+', Lockdown '+str(lockdown_date), fontsize=20)
        #ax.tick_params(axis='both', which='major', labelsize=15)
        #ax1.tick_params(axis='both', which='major', labelsize=15)
        ax.grid(linestyle='--', color='white')
        ax1.grid(linestyle='--', color='white')

        ax.legend(loc=0, fontsize=15, shadow=True, facecolor='white')
        ax1.legend(loc=0, fontsize=15, shadow=True, facecolor='white')
        ax.set_ylim(top=2*max(N))
        fig.autofmt_xdate(rotation=70)
        #ax1.spines['top'].set_color('k')
        [i.set_linewidth(1.5) for i in ax1.spines.itervalues()]
        plt.tight_layout() 
        plt.xlim([min(date_time), max(date_time)])
        plt.subplots_adjust(hspace=.0)
        fig.savefig('lockdown/'+y[1]+'_lock.png')
        #plt.show()
        plt.close()
f.close()




