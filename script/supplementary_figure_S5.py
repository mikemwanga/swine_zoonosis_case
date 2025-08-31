
# PACKAGES version
#kaleido==0.1.0
#pandas == 2.2.3
#plotly == 6.2.0
#numpy== 1.26.0


#import necessary libraries
import pandas as pd
import os, sys, glob
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import numpy as np
# import kaleido

# DEFINE CONSTANTS
black_color = 'rgb(0,0,0)'
textcolor='black'
gridcolor=' #ccd1d1 '
fontfamily= "sans-serif"
barcolor="#afb2b3"
bgcolor='rgba(0, 0, 0, 0)'
margin=dict(l=10, r=10, t=40, b=5)
tickfont=12
leftmost_axes = []
cols = 3  

#SET PATHS
# Path to the directory containing coverage data files
# CURRDIR=os.getcwd() #CALLING THE SCRIPT AS A PYTHON MODULE HENCE USES CURRENT DIRECTORY
path='../data/'


# #FUNCTION TO PLOT DATA
def plot_coverage(data):
    """Plot the coverage data using Plotly Express line plot."""


    fig = px.line(data, x = 'Position',y='log_depth', facet_col='gene', color='sample', width=1200,height=900,
                facet_col_wrap=3,range_y=[0,5],
                category_orders={"gene": ["PB2", "PB1", "PA", "HA",'NP','NA','MP','NS']})
    fig.update_yaxes(gridcolor=gridcolor,title='Read Depth (log10)', linecolor=gridcolor, mirror=True)
    fig.update_xaxes(matches=None, title=None, tickfont=dict(size=tickfont,color=textcolor),showticklabels=True)
    fig.update_xaxes(linecolor=gridcolor, mirror=True)
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    
    fig.update_layout(margin=margin,plot_bgcolor=bgcolor, paper_bgcolor=bgcolor,
                    legend=dict(title=None,orientation='v', x=0.8,y=.03,bordercolor="black",borderwidth=0.5),
                    font_family=fontfamily,font_color=textcolor,font_size=tickfont)
    
    for i, ax in enumerate(fig.select_yaxes(), start=1):
        if (i - 1) % cols == 0:
            leftmost_axes.append(ax)

    fig.for_each_yaxis(lambda y: y.update(title='Read Depth (log10)') if y in leftmost_axes else y.update(title=''))
    fig.write_image(os.path.join('../figure','Supplimentary_Figure_S5.pdf'))


#PROCESS THE DATA
def process_data():
    """Process the coverage data from text files and prepare it for plotting."""

    files = glob.glob(os.path.join(path, '*.txt'))
    depth_data = []

    for file in files:
        filename = file.split('/')[-1].split('.')[0]
        depth= pd.read_table(file, usecols=['Reference_Name','Position', 'Coverage Depth'])
        depth['sample'] = filename
        depth['gene'] = depth['Reference_Name'].str.split('_').str[1]
        depth_data.append(depth)
    
    merged_data = pd.concat(depth_data)
    merged_data['log_depth'] = np.log10(merged_data['Coverage Depth'].replace(0, 1e-3))

    return merged_data


#MAIN FUNCTION CALL
plot_coverage(process_data())