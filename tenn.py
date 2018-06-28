import numpy as np

from bokeh.io import curdoc
from bokeh.layouts import row,column, widgetbox
from bokeh.models import ColumnDataSource
from bokeh.models.callbacks import CustomJS
from bokeh.models.widgets import Slider,PreText, TextInput,Dropdown
from bokeh.plotting import figure
import numpy as np
import tau_tenn

# Set up data
N = 101
x = np.linspace(0, 1, N)
zeros = np.zeros(N)
source1 = ColumnDataSource(data=dict(x=x, y=zeros))
source2 = ColumnDataSource(data=dict(x=x, y=zeros))
source4 = ColumnDataSource(data=dict(x=x, y=zeros))
source5 = ColumnDataSource(data=dict(x=x, y=zeros))

# Set up plot
plot = figure(plot_height=200, plot_width=400, title="q",
              tools="pan,reset,save,wheel_zoom",
              x_range=[0, 1], y_range=[0, 10])
              
plot2 = figure(plot_height=200, plot_width=400, title="ne (10^20 m^-3)",
              tools="pan,reset,save,wheel_zoom",
              x_range=[0, 1], y_range=[0, 4])    

plot.line('x', 'y', source=source1, line_width=3, line_alpha=0.6)
plot2.line('x', 'y', source=source2, line_width=3, line_alpha=0.6)
              
              
plot3 = figure(plot_height=400, plot_width=400,
              tools="pan,reset,save,wheel_zoom",
              x_range=[0, 10], y_range=[-5, 5],toolbar_location=None)            
              
plot4 = figure(plot_height=200, plot_width=400,title='error',
              tools="pan,reset,save,wheel_zoom",
              x_range=[0, 100], y_range=[-1., 1.],toolbar_location=None)
            
plot5 = figure(plot_height=200, plot_width=400,title='Te',
              tools="pan,reset,save,wheel_zoom",
              x_range=[0, 1], y_range=[0, 25],toolbar_location=None)

plot4.line('x', 'y', source=source4, line_width=3, line_alpha=0.6)
plot5.line('x', 'y', source=source5, line_width=3, line_alpha=0.6)

a = 1.7/3.0

theta = np.arange(0,2*np.pi+2*np.pi*0.01,2*np.pi*0.01)
x = np.arcsin(0.7)
R = 1.7+a*np.cos(theta+x*np.sin(theta))
Z = 1.8*a*np.sin(theta)
source3 = ColumnDataSource(data=dict(x=R, y=Z))             
plot3.line('x', 'y', source=source3, line_width=3, line_alpha=0.6)

tokamaks =  [("DIII_D", "d3d"),("JET","JET") ,("ITER steady-state", "ITERSS"),("ARC","ARC"),("CAT-DEMO","CAT-DEMO"),("ARIES ACT 2","ACT2")]
starting_tok = Dropdown(label="Starting Tokamak", button_type="warning", menu=tokamaks)

majorRadius = Slider(title="Major Radius", value=1.7, start=0.5, end=10.0,step=0.1, callback_policy='mouseup')
Bt = Slider(title="Bt", value=2.2, start=1.0, end=12.0,step=0.1, callback_policy='mouseup')
delta = Slider(title="triangularity", value=0.7, start=0.0, end=1.0,step=0.05,callback_policy='mouseup')
kappa = Slider(title="elongation", value=1.8, start=1.0, end=3.0,step=0.05,callback_policy='mouseup')
aspect = Slider(title="Aspect Ratio", value=3.0, start=2.0, end=5.0,step=0.05, callback_policy='mouseup')
plasmaCurrent = Slider(title="Ip", value=1.0, start=0.5, end=25.0,step=0.05, callback_policy='mouseup')
fgw_ped = Slider(title="Greenwald fraction at pedestal", value=0.6, start=0.2, end=1.5,step=0.02, callback_policy='mouseup')

Paux = Slider(title="Aux power", value=7.0, start=1.0, end=100.,step=0.5, callback_policy='mouseup')
Zeff = Slider(title="Zeff", value=2.0, start=1.0, end=5.0,step=0.05, callback_policy='mouseup')
q0 = Slider(title="On-axis q", value=2.0, start=0.5, end=5.0,step=0.1, callback_policy='mouseup')
error = Slider(title="Error required for Convergence", value=0.05, start=0.01, end=0.1,step=0.01, callback_policy='mouseup')
densityPeaking = Slider(title="Density peaking", value=1.5, start=1.0, end=3.0,step=0.05, callback_policy='mouseup')
tempShape = Slider(title="Temp Shape", value=1.5, start=0.9, end=3.,step=0.1, callback_policy='mouseup')
Tratio = Slider(title="Ti/Te", value=1.0, start=0.1, end=4.0,step=0.05, callback_policy='mouseup')
stats = PreText(text='OUTPUTS', width=500)

def update_data(attrname, old, new):

    x = np.linspace(0, 1, 101) 
    S = 0.5*(1.+kappa.value**2*(1+0.3*delta.value))*(1+1.5/(aspect.value)**2)

    q95 = 5.*majorRadius.value/(aspect.value)**2*Bt.value/plasmaCurrent.value*S
    y = q0.value +(q95-q0.value)*x*x
    source1.data = dict(x=x, y=y)

    a = float(majorRadius.value)/float(aspect.value)
    n_gw = plasmaCurrent.value/(np.pi*a*a)
    ne_ped = fgw_ped.value*n_gw
    ne_core = densityPeaking.value*ne_ped
    
    nval,tval=tau_tenn.Hmode_profiles(nped14 = ne_ped,ncore14 = ne_core)   
    source2.data = dict(x=x, y=nval)
    

    theta = np.arange(0.0,2*np.pi+0.01*2*np.pi,0.01*2*np.pi)
    x = np.arcsin(delta.value)

    R = majorRadius.value+a*np.cos(theta+x*np.sin(theta))
    Z = kappa.value*a*np.sin(theta)
    source3.data = dict(x=R,y=Z)

    inputs={}
    inputs['Paux']=Paux.value
    inputs['aspect']=aspect.value
    inputs['Bt']=Bt.value
    inputs['Ip']=plasmaCurrent.value
    inputs['R']=majorRadius.value
    inputs['kappa']=kappa.value
    inputs['delta']=delta.value
    inputs['density_peaking']=densityPeaking.value
    inputs['tempShape']=tempShape.value
    inputs['zeff']=Zeff.value
    inputs['q0']=q0.value
    inputs['fgw_ped']=fgw_ped.value
    inputs['error']=error.value
    inputs['Tratio']=Tratio.value
    
    outputs =  tau_tenn.tau_tenn(inputs,tglf_model='nn')
    nval,tval=tau_tenn.Hmode_profiles(nped14 = ne_ped,ncore14 = ne_core,tcorekEV =outputs['T0'],tpedkEV=outputs['teped'],
                                      texpin=tempShape.value,texpout=tempShape.value)

    source4.data = dict(x=range(0,len(outputs['error'])),y=outputs['error'])
    x = np.linspace(0,1,101)
    source5.data = dict(x=x,y=tval)
    
    stats.text = 'OUTPUTS\n'
    try:
        stats.text += 'Pfus = ' +str(outputs['Pfus'])+   '\n\n'
        stats.text += 'Qplasma = ' +str(outputs['Qplasma'])+   '\n\n'
        stats.text += 'H98y2  = ' + str(outputs['H98'])+  '\n\n'
        stats.text += 'fbs = ' +str(outputs['fbs'])+ '\n\n'
        stats.text += 'Te,ped(keV) = ' + str(outputs['teped'])+ '\n\n'
        stats.text += 'taue(s) = ' +str(outputs['taue'])+ '\n\n'
        stats.text += 'betaN = '+ str(outputs['betaN'])+'\n\n' 
    except:
        stats.text += 'YOU BROKE ME!'
    return

def starting_tokamaks(attrname, old, new):
   if starting_tok.value == 'ITERSS':
       Paux.value = 73.0
       aspect.value = 3.1
       Bt.value = 5.3
       plasmaCurrent.value = 9.0
       majorRadius.value = 6.2
       kappa.value = 1.8
       delta.value = 0.5
       densityPeaking.value = 1.5
       Zeff.value = 2.0
       q0.value = 2.0
       fgw_ped.value = 0.8
       Tratio.value = 1.0

   if starting_tok.value == 'ARC':
       Paux.value = 39.0
       aspect.value = 3.0
       Bt.value = 9.2
       plasmaCurrent.value = 7.8
       majorRadius.value = 3.3
       kappa.value = 1.85
       delta.value = 0.5
       densityPeaking.value = 2.0
       Zeff.value = 2.0
       q0.value = 3.0
       fgw_ped.value = 0.4
       Tratio.value = 1.0

   if starting_tok.value == 'JET':
       [Paux.value,aspect.value,Bt.value,plasmaCurrent.value,majorRadius.value,kappa.value,delta.value,densityPeaking.value,Zeff.value,q0.value,fgw_ped.value,Tratio.value]=[30.,3.0,3.0,3.0,3.0,1.8,0.5,1.6,2.0,2.0,0.6,1.0]
   if starting_tok.value == 'd3d':
       Paux.value = 10.0
       aspect.value = 2.95
       Bt.value = 2.2
       plasmaCurrent.value = 1.0
       majorRadius.value = 1.7
       kappa.value = 1.8
       delta.value = 0.55
       densityPeaking.value = 1.5
       Zeff.value = 2.5
       q0.value = 2.0
       fgw_ped.value = 0.5
       Tratio.value = 1.0
   if starting_tok.value == 'CAT-DEMO':
       Paux.value = 60.0
       aspect.value = 3.0
       Bt.value = 7.0
       plasmaCurrent.value = 10.5
       majorRadius.value = 4.0
       kappa.value = 2.0
       delta.value = 0.7
       densityPeaking.value = 1.6
       Zeff.value = 2.0
       q0.value = 2.0
       fgw_ped.value = 1.0
       Tratio.value = 1.0

   if starting_tok.value == 'ACT2':
       Paux.value = 100.0
       aspect.value = 4.0
       Bt.value = 8.75
       plasmaCurrent.value = 14.0
       majorRadius.value = 9.75
       kappa.value = 2.2
       delta.value = 0.6
       densityPeaking.value = 1.6
       Zeff.value = 2.1
       q0.value = 2.0
       fgw_ped.value = 1.0
       Tratio.value = 1.0
starting_tok.on_change('value', starting_tokamaks)

for w in [delta,kappa,aspect,majorRadius,plasmaCurrent,densityPeaking,tempShape,Paux,Bt,Zeff,q0,fgw_ped,error,Tratio]:
    w.on_change('value', update_data)

# Set up layouts and add to document
inputs = widgetbox(starting_tok,delta,kappa,aspect,majorRadius,plasmaCurrent,densityPeaking,Bt,Paux,Zeff,q0,fgw_ped,error,Tratio)

plots = column(plot,plot2,plot3)
plots2 = column(plot4,plot5,stats)
curdoc().add_root(row(inputs, plots,plots2, width=1200))
curdoc().title = "tenn"
