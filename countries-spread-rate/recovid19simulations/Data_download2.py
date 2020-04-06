import urllib2

# https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data/csse_covid_19_time_series

for Str in ['confirmed', 'recovered', 'deaths']:
    filedata = urllib2.urlopen('https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_%s_global.csv'%(Str))
    datatowrite = filedata.read()
    with open('time_series_global-%s.csv'%(Str), 'w') as f:
        f.write(datatowrite)
    f.close()
    g = open('time_series_global-%s.csv'%(Str), 'r')
    h = open('time_series_global-%s_edit.csv'%(Str), 'w')
    for x in g:
        if "Korea, South" not in x:
           h.write(x)
        if "Korea, South" in x:
            y = x.replace('"Korea, South"' , 'South Korea')
            h.write(y)
    g.close()
    h.close()
