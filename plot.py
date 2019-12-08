import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd

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
			hover_data=['Ev', 'visc0', 'log10visc', 'Ra'], color_continuous_scale=px.colors.sequential.Bluered) #diverging.Spectral[::-1])
		plot.layout.xaxis.title.text='Time, Ga'
		plot.layout.yaxis.title.text='Temp, K'
		plot.update_traces(marker=dict(size=5.0, opacity=0.3))  #, color='rgba(1, 1, 1, 0.5)'))
	if colortype == 'discrete':
		df[colorcolumn] = df[colorcolumn].astype(str)
		plot = px.line(df, x="time", y="temp", color=colorcolumn, template = 'plotly_white+presentation',
			color_discrete_sequence=px.colors.qualitative.Vivid, hover_data=['planet', 'Ev', 'visc0', 'log10visc', 'Ra'])
		plot.layout.xaxis.title.text='Time, Ga'
		plot.layout.yaxis.title.text='Temp, K'
	if colortype == None:
		plot = px.line(df, x="time", y="temp", template = 'plotly_white+presentation',
			hover_data=['Ev', 'visc0', 'log10visc', 'Ra'])
                plot.update_traces(line=dict(width=3.0, color='rgba(1, 1, 1, 0.1)'))
		plot.layout.xaxis.title.text='Time, Ga'
                plot.layout.yaxis.title.text='Temperature, K'
	plot.layout.font.family='Arial'
        plot.update_xaxes(showline=True, ticks="inside", linewidth=2, linecolor='black', mirror=True, range=[0, 4.55])
        ymin, ymax = min(df['temp']), max(df['temp'])
	plot.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True, range=[ymin, ymax]) #[11, 18]) #
	plot.update_layout(showlegend=False) 
	#yaxis_type='log'
	plot.show()
	return plot

def plot_pd(df, p):
	a = df.plot(x = p['x'], y = p['y'], title = p['title'], s=0.5, xlim=(0, 4.5), ylim=(800,2200),
		legend = p['showlegend'], kind = 'scatter', color='blue', alpha=1.0, figsize=(5,5))
	plt.show()
	return None

def plot_pd_autolim(df, p):
        a = df.plot(x = p['x'], y = p['y'], title = p['title'], s=0.5, #xlim=(0, 4.5), ylim=(800,2200),
		legend = p['showlegend'], kind = 'scatter', color='red', alpha=1.0)
	plt.show()
        return None
