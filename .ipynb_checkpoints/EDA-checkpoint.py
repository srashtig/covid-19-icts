import pandas as pd
import numpy as np
from matplotlib import pyplot as plt



# In[2]:


df = pd.read_csv('covid_19_india.csv')
df.drop(columns=['Sno','ConfirmedIndianNational','ConfirmedForeignNational','Time'],inplace=True)


df=df.groupby(['State/UnionTerritory',"Date"]).head()
df=df.dropna()
States=np.unique(df['State/UnionTerritory'].values)
States=States[States!='Unassigned']
States=States[States!='Nagaland']
States=States[States!='Nagaland#']


# In[6]:


top5aff_states=df.groupby(['State/UnionTerritory']).max().sort_values(['Confirmed'],ascending=False)[:15].index.values
#print(top5aff_states)


# In[6]:


dates=df[df['State/UnionTerritory'] == 'Kerala']['Date'].values
refined_df_confirmed=pd.DataFrame(columns=dates,index=States)
refined_df_cured=pd.DataFrame(columns=dates,index=States)
refined_df_deaths=pd.DataFrame(columns=dates,index=States)


# In[7]:



for state in States:
    try:
        df1=df[df['State/UnionTerritory'] == state]
        rec_date_idx=np.where(dates==df1['Date'].values[0])[0][0]
        if rec_date_idx >0:
            df2=pd.DataFrame()
            df2['Date']=dates[:rec_date_idx]
            df2['Confirmed'] =  np.zeros(rec_date_idx)
            df2['Cured']=np.zeros(rec_date_idx)
            df2['Deaths']=np.zeros(rec_date_idx)
            df2['State/UnionTerritory']=state
            df2=df2.append(df1,ignore_index=True)
        else: 
            df2=df1
        refined_df_confirmed.loc[state]=df2['Confirmed'].values
        refined_df_cured.loc[state]=df2['Cured'].values

        refined_df_deaths.loc[state]=df2['Deaths'].values
        
    except:
        print(state + ' not recorded')

    #df2.to_csv(state+'.csv',index=False)
refined_df_confirmed.index.name='State/UT'
refined_df_cured.index.name='State/UT'
refined_df_deaths.index.name='State/UT'

refined_df_cured.to_csv('refined_df_cured.csv')
refined_df_deaths
    
    
refined_df_confirmed.to_csv('refined_df_confirmed.csv')
refined_df_cured.to_csv('refined_df_cured.csv')
refined_df_deaths.to_csv('refined_df_deaths.csv')






