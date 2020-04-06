
# ## covid-19 Data exploration from data set available by John hopkins CSSE

# In[1]:


import pandas as pd
import pandas as pd
import numpy as np
import os
import matplotlib.pylab as plt
import scipy
#import plotly.express as px


# # importing time series datasets of confirmed, recovered and deaths

# In[2]:


url='https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/'
df_confirmed = pd.read_csv(url + 'time_series_covid19_confirmed_global.csv',error_bad_lines=False)  
df_recovered = pd.read_csv(url + 'time_series_covid19_recovered_global.csv',error_bad_lines=False) 
df_deaths = pd.read_csv(url + 'time_series_covid19_deaths_global.csv',error_bad_lines=False) 
df_confirmed.rename(columns={'Country/Region':'Country'}, inplace=True)
df_recovered.rename(columns={'Country/Region':'Country'}, inplace=True)
df_deaths.rename(columns={'Country/Region':'Country'}, inplace=True)

## date upto which analysis is to be considered and name the output directory
date='4/2/20'
odir = date.replace('/','_') + '/'
try:  
    os.mkdir(odir)  
except OSError as error:  
    print("Directory already exists")


# In[3]:


df_confirmed.head()


# In[4]:


df_confirmed=df_confirmed.groupby(["Country"]).sum()
df_confirmed=df_confirmed.drop(columns=['Lat','Long'])
df_recovered=df_recovered.groupby(["Country"]).sum()
df_recovered=df_recovered.drop(columns=['Lat','Long'])
df_deaths=df_deaths.groupby(["Country"]).sum()
df_deaths=df_deaths.drop(columns=['Lat','Long'])
df_confirmed.head()


## extracting top 60 infected countries

# In[5]:


top_affected_countries=df_confirmed.sort_values([date],ascending=False).index[:40].values #20/20
#top_affected_countries=np.delete(top_affected_countries,[62,59,58,57,55])

top_affected_countries


# In[6]:


np.savetxt(odir+'top_affected_countries.txt',top_affected_countries, delimiter=" ", fmt="%s")


# In[7]:


df_confirmed_top=df_confirmed.loc[top_affected_countries,:]
df_recovered_top=df_recovered.loc[top_affected_countries,:]
df_deaths_top=df_deaths.loc[top_affected_countries,:]
df_confirmed_top.to_csv(odir+'confirmed-data-top-affected-countries.csv')
df_confirmed_top.head()


# In[8]:


df_active_top=df_confirmed_top-df_recovered_top - df_deaths_top
df_active1_top=df_confirmed_top-df_recovered_top

#df_active_top.to_csv(odir+'active-data-top-affected-countries.csv')
#df_active1_top.to_csv(odir+'active1-data-top-affected-countries.csv')


# In[9]:


plt.figure(figsize=(10,10))
plt.title('Total Number of Cases',fontsize=15)
y_pos=np.arange(0,len(top_affected_countries))
p1 =plt.barh(y_pos,width=np.flip((df_confirmed_top[date].values)),align='center',label='Confirmed')
p1 =plt.barh(y_pos,width=np.flip((df_active_top[date].values)),align='center', label='Active')
plt.yticks(y_pos,np.flip(top_affected_countries),rotation=0,fontsize=13)
plt.ylim(0,len(top_affected_countries))
plt.legend(fontsize=14)
plt.xlabel('Log Scale',fontsize=12)
plt.xscale('log')
plt.grid()
#plt.xlim([9e4,0])
plt.tight_layout()
plt.savefig(odir+'Confirmed_active_hor.png')
#pd.plotting.table(data=df_confirmed_top['3/21/20'])


# In[10]:


confirmed_mat=df_confirmed_top.values
deaths_mat=df_deaths_top.values
recovered_mat=df_recovered_top.values
active_mat = df_active_top.values


# In[11]:


plt.figure(figsize=(10,7))
idx=np.where(top_affected_countries== 'India')[0][0]
[plt.plot(np.log10(confirmed_mat[i,confirmed_mat[i,:]>50]),'.-',label=top_affected_countries[i]) for i in range(0,40,3)]
plt.plot(np.log10(confirmed_mat[idx,confirmed_mat[idx,:]>50]),'k*-',label=top_affected_countries[idx])
plt.legend(fontsize=12)
plt.xlabel('No. of days(start from cases = 50)',fontsize=12)
plt.ylabel('Log(Cases)',fontsize=12)
#plt.yscale('log')
plt.tight_layout()
#plt.savefig(odir+'timeseries.png')


## slope of log(#cases) v/s days: spread rate in logarithmic phase

# In[12]:


def calc_slopes(confirmed_mat,thr=3,guess=1,extra=False):
    #slopes=np.zeros(50)
    l=len(top_affected_countries)
    fit_params = np.zeros([l,5])
    for i in range(l):
        idx=np.where(confirmed_mat[i,:]>=thr)[0][:21]
        #print(idx)
        y = np.log(confirmed_mat[i,idx.T])
        if len(y)>0:
            x = np.arange(len(y))
            z, res, _, _, _ = np.polyfit(x, y, guess,full=True)
            
            fit_params[i,0]=z[0]
            fit_params[i,1]=z[1]
            fit_params[i,2]=idx[0]
            fit_params[i,3]= np.log(confirmed_mat[i,idx[0]])
            fit_params[i,4]=res                        
    return fit_params
    


# In[13]:


fit_params_confirmed=calc_slopes(confirmed_mat,thr=50)
fit_params_active=calc_slopes(active_mat,thr=50)


# In[14]:


#make dataframe  
def make_df_cases(fit_params_confirmed):
    confirmed_fit_df= pd.DataFrame()
    confirmed_fit_df['Country']=top_affected_countries
    confirmed_fit_df['slope(log cases)'] = fit_params_confirmed[:,0]
    confirmed_fit_df['start day'] = fit_params_confirmed[:,2].astype(int)
    confirmed_fit_df['start date']=df_confirmed_top.columns[confirmed_fit_df['start day']]
    confirmed_fit_df['intersecpt(log cases)'] = fit_params_confirmed[:,1]
    confirmed_fit_df['MSE fit'] = fit_params_confirmed[:,3]
    confirmed_fit_df['total cases'] =df_confirmed_top[date].values
    confirmed_fit_df['total deaths']=df_deaths_top[date].values
    confirmed_fit_df['total recovered']=df_recovered_top[date].values
    confirmed_fit_df['mortality'] = confirmed_fit_df['total deaths']/confirmed_fit_df['total cases']
    confirmed_fit_df=confirmed_fit_df.sort_values('slope(log cases)')
    #confirmed_fit_df=confirmed_fit_df[confirmed_fit_df.Country != 'China']
    #confirmed_fit_df=confirmed_fit_df[confirmed_fit_df.Country != 'Korea, South']
    #confirmed_fit_df=confirmed_fit_df[confirmed_fit_df.Country != 'Cruise Ship']
    #confirmed_fit_df=confirmed_fit_df[confirmed_fit_df.Country != 'United Arab Emirates']

    return confirmed_fit_df


## sorted from highest to lowest spread rate of confirm cases

# In[15]:


confirmed_fit_df=make_df_cases(fit_params_confirmed)
#confirmed_fit_df.to_csv(r'47_countries_confirmed_fit.csv', index = False)
plt.figure(figsize=(10,10))
plt.title('COVID-19: Rate of Spreading ',fontsize=20)
plt.xlabel('Slope of Ln(Cases) v/s Days',fontsize=15)
y_pos=np.arange(0,len(confirmed_fit_df['Country'].values))
plt.barh(y_pos,width=(confirmed_fit_df['slope(log cases)'].values),align='center',color='orange')
plt.yticks(y_pos,(confirmed_fit_df['Country'].values),rotation=0,fontsize=13)
plt.ylim(0,len(confirmed_fit_df['Country'].values))
plt.tight_layout()
plt.grid()
plt.savefig(odir+'spread_rate_bar_hor.png')


# In[16]:


# plot settings
from matplotlib import rc
import matplotlib
#matplotlib.rc('text.latex', preamble='\\usepackage{txfonts}')


#rc('text', usetex=True)
#rc('font', family='serif')
#rc('font', serif='times')
rc('mathtext', default='sf')
rc("lines", markeredgewidth=1)
rc("lines", linewidth=2)
rc('axes', labelsize=16)  # 24
rc("axes", linewidth=0.5)  # 2)
rc('xtick', labelsize=14)
rc('ytick', labelsize=14)
rc('legend', fontsize=10)  # 16
rc('xtick.major', pad=6)  # 8)
rc('ytick.major', pad=6)  # 8)
rc('xtick.minor', size=5)  # 8)
rc('ytick.minor', size=5)  # 8)


# In[17]:


# import temperature data
temp_df=pd.read_csv('GlobalLandTemperaturesByCountry.csv')#('training_data_with_weather_info_week_1.csv')
#temp=temp_df.groupby(['Country/Region']).mean()['temp']
temp=temp_df.groupby(['Country']).mean()['AverageTemperature']
confirmed_fit_df['temp'] = temp[confirmed_fit_df['Country'].values].values
#confirmed_fit_df['temp'] = (temp[confirmed_fit_df['Country'].values].values -32) * 5./9
#confirmed_fit_df.to_csv(r'top_countries_confirmed_fit.csv', index = False)
fig,ax=plt.subplots(figsize=(10,10))
confirmed_fit_df.plot('temp', 'slope(log cases)', kind='scatter', ax=ax,color='red')
confirmed_fit_df[['temp','slope(log cases)','Country']].apply(lambda row: ax.text(*row,fontsize=14),axis=1);
plt.xlabel('Avg. Temperature (degrees Celsius)')
plt.ylabel('Rate of Spreading(slope of ln(cases))')
plt.grid()
#import scipy
#x,y=confirmed_fit_df['temp'].values, confirmed_fit_df['slope(log cases)'].values
#nas = np.logical_or(np.isnan(x), np.isnan(y))
#coeff,p_val=scipy.stats.pearsonr(x[~nas], y[~nas])
#plt.title('Pearson Coeff = %.2f' %(coeff)+ ' and  p= %.2f' %(p_val),fontsize=14)
plt.savefig(odir+'temp_vs_slope.pdf')


# In[18]:


fig,ax=plt.subplots(figsize=(10,10))
confirmed_fit_df.plot('temp', 'total cases', kind='scatter', ax=ax,color='red')
confirmed_fit_df[['temp','total cases','Country']].apply(lambda row: ax.text(*row,fontsize=14),axis=1);
plt.xlabel('Avg. Temperature (Celsius)')
plt.ylabel('Total Cases')
#plt.ylim(0,1e5)
plt.yscale('log')
plt.grid()
x,y=confirmed_fit_df['temp'].values, confirmed_fit_df['total cases'].values
#nas = np.logical_or(np.isnan(x), np.isnan(y))
#coeff,p_val=scipy.stats.pearsonr(x[~nas], y[~nas])
plt.tight_layout()
#plt.title('Pearson Coeff = %.2f' %(coeff)+ ' and  p= %.2f' %(p_val),fontsize=14)
plt.savefig(odir+'temp_vs_TotalCases_log.png')


# In[21]:


fig,ax=plt.subplots(figsize=(10,10))
confirmed_fit_df.plot('temp', 'total deaths', kind='scatter', ax=ax,color='red')
confirmed_fit_df[['temp','total deaths','Country']].apply(lambda row: ax.text(*row,fontsize=14),axis=1);
plt.xlabel('Avg. Temperature (Celsius)')
plt.ylabel('Total Deaths')
plt.yscale('log')
plt.grid()
x,y=confirmed_fit_df['temp'].values, confirmed_fit_df['total deaths'].values
#nas = np.logical_or(np.isnan(x), np.isnan(y))
#coeff,p_val=scipy.stats.pearsonr(x[~nas], y[~nas])
#plt.title('Pearson Coeff = %.2f' %(coeff)+ ' and  p= %.2f' %(p_val),fontsize=14)
plt.savefig(date.replace('/','_')+'temp_vs_TotalDeaths.png')


# In[22]:


fig,ax=plt.subplots(figsize=(10,10))
confirmed_fit_df.plot('temp', 'mortality', kind='scatter', ax=ax,color='red')
confirmed_fit_df[['temp','mortality','Country']].apply(lambda row: ax.text(*row,fontsize=14),axis=1);
plt.xlabel('Avg. Temperature (Celsius)')
plt.ylabel('Mortality (fraction of death cases from total)')
#plt.yscale('log')
plt.grid()
x,y=confirmed_fit_df['temp'].values, confirmed_fit_df['mortality'].values
#nas = np.logical_or(np.isnan(x), np.isnan(y))
#coeff,p_val=scipy.stats.pearsonr(x[~nas], y[~nas])
#plt.title('Pearson Coeff = %.2f' %(coeff)+ ' and  p= %.2f' %(p_val),fontsize=14)
plt.savefig(odir+'temp_vs_MortalityFraction.png')

