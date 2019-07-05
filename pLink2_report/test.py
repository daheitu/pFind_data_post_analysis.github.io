import plotly.plotly as py
import plotly.tools as tls
import plotly
import matplotlib.pyplot as plt
plotly.tools.set_credentials_file(username='daheitu', api_key='••••••••••')
bubbles_mpl = plt.figure()

# doubling the width of markers
x = [0,2,4,6,8,10]
y = [0]*len(x)
s = [20*4**n for n in range(len(x))]
plt.scatter(x,y,s=s)

#plotly_fig = tls.mpl_to_plotly(bubbles_mpl)
#py.iplot(plotly_fig, filename='mpl-bubbles')