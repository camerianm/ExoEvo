import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np

def pseudo_heat(df, xval, yval, colorcolumn, colortype):
	# Input: nparray = Numpy array; columnkeys = list of strings; colorcolumn = col name; colortype = continuous or discrete
	# df = pd.DataFrame(data=nparray, columns=columnkeys)
	if colortype == "continuous":
		plot = px.scatter(df, x=xval, y=yval, opacity = 0.5, color=colorcolumn, template = 'plotly_white+presentation',
			hover_data=['Ev', 'visc0', 'Urey', 'temp'], color_continuous_scale=px.colors.sequential.Rainbow) #diverging.Spectral[::-1])
		plot.layout.xaxis.title.text=xval
		plot.layout.yaxis.title.text=yval
		plot.update_traces(marker=dict(size=30.0))  #, color='rgba(1, 1, 1, 0.5)'))
	# if colortype == 'discrete':
	# 	df[colorcolumn] = df[colorcolumn].astype(str)
	# 	plot = px.line(df, x="time", y="temp", color=colorcolumn, template = 'plotly_white+presentation',
	# 		color_discrete_sequence=px.colors.qualitative.Vivid, hover_data=['Ev', 'visc0', 'log10visc', 'Ra'])
	# 	plot.layout.xaxis.title.text='Time, Ga' 
	# 	plot.layout.yaxis.title.text='Temp, K'

	plot.layout.font.family='Arial'
	ymin, ymax = 800, 2200 #min(df[yval]), max(df[yval])
	xmin, xmax = min(df[xval]), max(df[xval])
	plot.update_xaxes(showline=True, ticks="inside", linewidth=2, linecolor='black', mirror=True, range=[xmin, xmax]) #, range=[0, 4.5])
	plot.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True, range=[ymin, ymax]) #[11, 18]) #
	plot.update_layout(showlegend=False) 
	#yaxis_type='log')
	plot.show()
	return df


def evolution_colorcoded(df, colorcolumn, colortype):
	# Input: nparray = Numpy array; columnkeys = list of strings; colorcolumn = col name; colortype = continuous or discrete
	# df = pd.DataFrame(data=nparray, columns=columnkeys)
	if colortype == "continuous":
		plot = px.scatter(df, x="time", y="temp", opacity = 0.05, color=colorcolumn, template = 'plotly_white+presentation',
			hover_data=list(df.keys()), color_continuous_scale=px.colors.sequential.Bluered) #diverging.Spectral[::-1])
		plot.layout.xaxis.title.text='Time, Ga'
		plot.layout.yaxis.title.text='Temp, K'
		plot.update_traces(marker=dict(size=5.0, opacity=0.3))  #, color='rgba(1, 1, 1, 0.5)'))
	if colortype == 'discrete':
		df[colorcolumn] = df[colorcolumn].astype(str)
		plot = px.line(df, x="time", y="temp", color=colorcolumn, template = 'plotly_white+presentation',
			color_discrete_sequence=px.colors.qualitative.Vivid, hover_data=list(df.keys()))
		plot.layout.xaxis.title.text='Time, Ga'
		plot.layout.yaxis.title.text='Temp, K'
	if colortype == None:
		plot = px.line(df, x="time", y="temp", template = 'plotly_white+presentation',
			hover_data=list(df.keys()))
		plot.update_traces(line=dict(width=3.0, color='rgba(1, 1, 1, 0.1)'))
		plot.layout.xaxis.title.text='Time, Ga'
		plot.layout.yaxis.title.text='Temperature, K'
	plot.layout.font.family='Arial'
	plot.update_xaxes(showline=True, ticks="inside", linewidth=2, linecolor='black', mirror=True, range=[0, list(df['time'])[-1]])
	ymin, ymax = min(df['temp']), max(df['temp'])
	plot.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True, range=[ymin, ymax]) #[11, 18]) #
	plot.update_layout(showlegend=False) 
	#yaxis_type='log'
	plot.show()
	return plot

def plot_pd_mylimits(df, p):
	#input required: pandas data frame, plus dictionary formatted as follows (example case)
	# p1 = {'x': 'time', 'y': 'temp', 'showlegend': False, 'title': 'Thermal history sample case', 'xlim': (0.0, 4.55), 'ylim': (800,2000)}
	a = df.plot(x = p['x'], y = p['y'], title = p['title'], s=0.5, xlim=xlim, ylim=ylim,
		legend = p['showlegend'], kind = 'scatter', color='8C1D40', alpha=1.0, figsize=(5,5))
	fname = p['x']+p['y']+'.png'
	plt.savefig(fname)
	plt.show()
	return None

def plot_pd_autolimit(df, p):
    #input required: pandas data frame, plus dictionary formatted as follows (example case)
    # p1 = {'x': 'time', 'y': 'temp', 'showlegend': False, 'title': 'Thermal history sample case'}
    if np.log10(np.max(np.abs(df[p['y']]))) - np.log10(np.min(np.abs(df[p['y']]))) > 3: # if the y limits span >3 orders of magnitude,
        a = df.plot(x = p['x'], y = p['y'], title = p['title'], logy=True, legend = p['showlegend'], kind = 'scatter', color='8C1D40', alpha=0.8, s=0.3, figsize=(8, 8))
    else:
        a = df.plot(x = p['x'], y = p['y'], title = p['title'], logy=False, legend = p['showlegend'], kind = 'scatter', color='8C1D40', alpha=0.8, s=0.3, figsize=(8, 8))
    fname = p['x']+p['y']+'passfail.png'
    plt.savefig(fname)
    plt.show()
    return a


def plot_boolean(df, p):
	#input required: pandas data frame with a boolean column (here called 'passfail' - will change to slicing, etc. later
	# ...plus dictionary formatted as follows (example case)
	# p1 = {'x': 'time', 'y': 'temp', 'showlegend': False, 'title': 'Thermal history sample case', 'colorcolumn': 'passfail'}
    y_logscale, x_logscale = False, False #automatically assumes linear scale 
    plt.scatter(df[df[colorcolumn]==False][p['x']], df[df[colorcolumn]==False][p['y']],  c='#8C1D40', alpha=0.8, s=0.3)
    plt.scatter(df[df[colorcolumn]==True][p['x']], df[df[colorcolumn]==True][p['y']],  c='#FFC627',  alpha=0.8, s=0.3)
    if np.log10(np.max(np.abs(df[p['y']]))) - np.log10(np.min(np.abs(df[p['y']]))) > 3: # if the y limits span >3 orders of magnitude,
        if p['y'] != 'temp':
            plt.yscale('log')
    if np.log10(np.max(np.abs(df[p['x']]))) - np.log10(np.min(np.abs(df[p['x']]))) > 3: # if the x limits span >3 orders of magnitude,
        if p['x'] != 'time':
            plt.xscale('log')
    fname = p['x']+p['y']+colorcolumn+'.png'
    plt.savefig(fname)
    plt.show()
    return fname
