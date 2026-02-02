import pandas as pd
import glob
import numpy as np
import plotly.express as px

#generate segment specific coverage plots.
path='../data/coverage_data/'

files=glob.glob(path+'*.txt')

coverage_dfs=[]
for file in files:
    df=pd.read_table(file,sep='\t', names=['Reference_Name','Position',
                                             'Consensus_Count'])
    
    #have coverage in log10 using numpy
    df['Coverage (log10)']=np.log10(df['Consensus_Count']+1).round(1)

    #split to get gene and samplename
    df['sample_name']=df['Reference_Name'].str.split('_').str[0]
    df['sample_name']=df['sample_name'].map({'A1':'A0001','A2':'A0002','A3':'A0003',
                                             'A4':'A0004','A5':'A0005','H1':'H0001'})
    df['gene']=df['Reference_Name'].str.split('_').str[1]
    #treat NA as segment character not missing data
    df['gene']=df['gene'].fillna('NA.')
    # #drop where sample_name is missing
    coverage_dfs.append(df)
    # print(df.head() )
coverage_df=pd.concat(coverage_dfs)
# print(coverage_df.head())

#generate faceted figures by segment 
bgcolor='rgba(0, 0, 0, 0)'
margin=dict(l=50, r=10, t=40, b=50)
gridcolor='#ccd1d1'
tickfont=22
textsize=16
textcolor='black'

color_discrete_map={'A0001':"#652e2f",'A0002':'#984547','A0003':'#a94d4f','A0004':'#c28283',
                'A0005':'#dcb7b8','H0001':"#0745f0"}
figure = px.line(coverage_df, x='Position', y='Coverage (log10)', 
                 color='sample_name',facet_col='gene', facet_col_wrap=4,
                 category_orders={'sample_name':['A0001','A0002','A0003','A0004','A0005','H0001']},
                 facet_row_spacing=0.2,color_discrete_map=color_discrete_map)
figure.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))

#specify line thickness
figure.update_traces(line=dict(width=3))

figure.update_yaxes(gridcolor=gridcolor, linecolor=gridcolor, mirror=True, range=[0,6],
                     zeroline=True,zerolinecolor=gridcolor,tickfont=dict(size=tickfont,color=textcolor), ticks='outside')
figure.update_xaxes(matches=None, title=None, 
                    tickfont=dict(size=tickfont,color=textcolor),showticklabels=True,
                     ticks='outside')
figure.update_xaxes(linecolor=gridcolor, mirror=True)


figure.update_layout(margin=margin,plot_bgcolor=bgcolor, paper_bgcolor=bgcolor,
                    legend=dict(title=None,orientation='h', 
                                x=0.3,y=-0.1,font=dict(size=18, color='black')),
                    font_family='Arial',font_color='black',font_size=textsize)

figure.write_image('../figure/Figure_S2.pdf', width=1200, height=700)
