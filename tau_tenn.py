import sys
sys.path.append("/home/mcclenaghanj/neural")
import numpy as np
from brainfusetf import btf_connect
import numpy as np
import scipy
from scipy import integrate
import copy
import time

#defaultVars(inputs=OMFIT['inputs'],tglf_model='full') #Line for running in OMFIT

# Function builds input.tglf file and runs TGLF or TGLF-NN
def run_tglf(inputs,outputs,tglf_model,ne,Te,gridpoint,betap):

    Rmaj0 =inputs['R']*1.e2 # m to cm
    a = outputs['a']*1.e2 # m to cm
    Bt = inputs['Bt']*1.e4 # T to Gauss
    q0 = inputs['q0']
    q95 =outputs['q95']
    zeff = inputs['zeff']
    kappa1 = inputs['kappa']
    delta1 = inputs['delta']

    r = np.linspace(0,a,101)
    dr = r[1]-r[0]

    # Shafranov shift estimated as linear interpolation between R0 - (R0+betap*r*r/R), ~10% underestimated comparing to geqdsk
    Rmaj = np.ones(101)*(Rmaj0+(betap+0.5)*r[-1]*r[-1]/Rmaj0)-(betap+0.5)*r*r/Rmaj0
    dRmaj = np.gradient(Rmaj,dr)
    # New shifted elipse estimation of dRmaj

    kappa0 = 0.5*(1.+kappa1)
    kappa = kappa0*np.ones(101)+(kappa1-kappa0)*r**2/a**2
    delta = 0.5*delta1*(r/a+r**2/a**2)

    dlogTdr = np.gradient(Te,dr)/Te
    dlogndr = np.gradient(ne,dr)/ne

    skappa = r*np.gradient(kappa,dr)/kappa
    sdelta = r*np.gradient(delta,dr)

    q = q0 +(q95-q0)*r**2/(0.95*a)**2

    #rough estimate of psi for elongated plasma (found to be 1.11 times overestimated compared to geqdsk)
    #This error could cause ~10% difference in flux when converting to Gyro-Bohm units to real units

    #psi = kappa[1:]*Bt*scipy.integrate.cumtrapz(r/q,dx=dr)
    #psi = np.append(0.0,psi)
    #bunit = q/r*np.gradient(psi,dr)
    # New Miller psi estimate based on Stacey et al 2009
    bunit = np.zeros(101)
    for i in range(1,101):
        theta = np.linspace(0,2.0*np.pi,100)
        dtheta = theta[1]-theta[0]
        x = np.arcsin(delta[i])
        sdelta_stacey = sdelta[i]/np.sqrt(1.-delta[i]**2)
        bigtheta= np.sqrt(np.sin(theta +x*np.sin(theta))**2*(1.+x*np.cos(theta))**2+kappa[i]**2*np.cos(theta)**2)
        bigtheta*=1.0/(np.cos(x*np.sin(theta))+2.0*dRmaj[i]*np.cos(theta)+
                       (skappa[i]-sdelta_stacey*np.cos(theta)+(1.+skappa[i])*x*np.cos(theta))*np.sin(theta)*np.sin(theta+x*np.sin(theta)))
        line_integrand = 1.0/(1+r[i]/Rmaj[i]*np.cos(theta+x*np.sin(theta)))/bigtheta
        Ztmp = kappa[i]*r[i]*np.sin(theta)
        Rtmp = Rmaj[i]+r[i]*np.cos(theta+x*np.sin(theta))
        dRdt = np.gradient(Rtmp,dtheta)
        dZdt = np.gradient(Ztmp,dtheta)
        area = scipy.integrate.cumtrapz(np.sqrt(dRdt**2 + dZdt**2)*line_integrand, theta)
        bunit[i] = Bt*kappa[i]*(area[-1]/2.0/np.pi/r[i])
    bunit[0]=bunit[1]

    if tglf_model=='nn':
        input_tglf={}
    else:
        input_tglf=OMFIT['TGLF_GACODE']['FILES']['input.tglf']
    input_tglf['VEXB_SHEAR'] = 0.0

    mixing='dt'
    if mixing =='dt':
        input_tglf['NS']=5
        input_tglf['MASS_1'] = 0.000272313
        input_tglf['MASS_2'] = 1.
        input_tglf['MASS_3'] = 1.5
        input_tglf['MASS_4'] = 2.
        input_tglf['MASS_5'] = 10.

        for i in range(1,6):
            input_tglf['VPAR_'+str(i)] = 0.0
            input_tglf['VPAR_SHEAR_'+str(i)] = 0.0
            input_tglf['RLNS_'+str(i)] = -a*dlogndr[gridpoint]
            input_tglf['RLTS_'+str(i)] = -a*dlogTdr[gridpoint]

        input_tglf['TAUS_1'] = 1.
        input_tglf['TAUS_2'] = inputs['Tratio']
        input_tglf['TAUS_3'] = inputs['Tratio']
        input_tglf['TAUS_4'] = inputs['Tratio']
        input_tglf['TAUS_5'] = inputs['Tratio']

        input_tglf['ZS_2'] = 1
        input_tglf['ZS_3'] = 1
        input_tglf['ZS_4'] = zimp1 =  2
        input_tglf['ZS_5'] = zimp2 =  10
        impurityFraction2 = (zeff-1.)/(5*zimp1**2+zimp2**2-5*zimp1-zimp2)# setting He impurity density to 5 times Neon
        impurityFraction1 = 5.0*impurityFraction2
        input_tglf['AS_1'] = AS_1 = 1.
        input_tglf['AS_2'] = AS_2 = 0.5*(1.-zimp2*impurityFraction2-zimp2*impurityFraction2)
        input_tglf['AS_3'] = AS_3 = 0.5*(1.-zimp2*impurityFraction2-zimp2*impurityFraction2)
        input_tglf['AS_4']= AS_4  = impurityFraction1
        input_tglf['AS_5']= AS_5  = impurityFraction2
    else:
        input_tglf['NS'] = 3
        input_tglf['MASS_1'] = 0.000272313
        input_tglf['MASS_2'] = 1.
        input_tglf['MASS_3'] = 6.0

        for i in range(1,4):
            input_tglf['VPAR_'+str(i)] = 0.0
            input_tglf['VPAR_SHEAR_'+str(i)] = 0.0
            input_tglf['RLNS_'+str(i)] = -a*dlogndr[gridpoint]
            input_tglf['RLTS_'+str(i)] = -a*dlogTdr[gridpoint]

        input_tglf['TAUS_1'] = 1.
        input_tglf['TAUS_2'] = inputs['Tratio']
        input_tglf['TAUS_3'] = inputs['Tratio']


        input_tglf['TAUS_4'] = inputs['Tratio']
        input_tglf['TAUS_5'] = inputs['Tratio']

        zimp1 = 6.0
        impurityFraction1=(zeff-1.)/(zimp1*(zimp1-1.))
        input_tglf['AS_1'] = AS_1 = 1.
        input_tglf['AS_2'] = AS_2 = 1. - zimp1*impurityFraction1
        input_tglf['AS_3']= AS_3  = impurityFraction1
        AS_4=AS_5=0.0

    e = 4.8e-10
    k = 1.6e-12
    me = 2.0*1.6e-24*input_tglf['MASS_1']
    mi = 2.0*1.6e-24*input_tglf['MASS_2']
    c = 3.e10
    c_s = np.sqrt(k*Te/mi)
    input_tglf['BETAE'] =  8.*np.pi*ne[gridpoint]*k*Te[gridpoint]/bunit[gridpoint]**2
    loglam = 24.0-np.log(np.sqrt(ne[gridpoint]/Te[gridpoint]))
    input_tglf['XNUE'] = r[gridpoint]/c_s[gridpoint]*e**4*np.sqrt(2.)*np.pi*ne[gridpoint]*loglam/(np.sqrt(me)*(k*Te[gridpoint])**1.5)
    input_tglf['ZEFF'] = zeff
    rho_s  = c_s/(e*bunit/(mi*c))
    input_tglf['DEBYE'] = 7.43e2*np.sqrt(Te[gridpoint]/(ne[gridpoint]))/rho_s[gridpoint]

    input_tglf['RMIN_LOC'] = r[gridpoint]/a
    input_tglf['RMAJ_LOC'] = Rmaj[gridpoint]/a
    input_tglf['ZMAJ_LOC'] = 0.
    input_tglf['DRMINDX_LOC'] = 1.
    input_tglf['DRMAJDX_LOC'] = dRmaj[gridpoint]/a
    input_tglf['DZMAJDX_LOC'] = 0.

    input_tglf['Q_LOC']  = q[gridpoint]
    input_tglf['KAPPA_LOC'] = kappa[gridpoint]
    input_tglf['S_KAPPA_LOC'] = skappa[gridpoint]
    input_tglf['DELTA_LOC'] = delta[gridpoint]
    input_tglf['S_DELTA_LOC'] = sdelta[gridpoint]
    input_tglf['ZETA_LOC'] = 0.0
    input_tglf['S_ZETA_LOC'] = 0.0

    press = ne*Te*(AS_1+AS_2+AS_3+AS_4+AS_5)
    dpdr = np.gradient(press,dr)/press
    dqdr = np.gradient(q,dr)
    input_tglf['P_PRIME_LOC'] = q[gridpoint]/(r[gridpoint]/a)**2*r[gridpoint]*input_tglf['BETAE']/(8.*np.pi)*(
             (AS_1+AS_2*input_tglf['TAUS_2']+AS_3*input_tglf['TAUS_3']+AS_4*input_tglf['TAUS_4']+AS_5*input_tglf['TAUS_5'])*dpdr[gridpoint])
    s = r[gridpoint]/q[gridpoint]*dqdr[gridpoint]
    input_tglf['Q_PRIME_LOC'] = q[gridpoint]**2*a**2/r[gridpoint]**2*s

    q_gb = ne*k*Te*c_s*(rho_s/a)**2*1.0e-7

    nn_inputs = np.atleast_2d([[input_tglf['AS_2'], input_tglf['AS_3'], input_tglf['BETAE'],  input_tglf['DEBYE'],input_tglf['DELTA_LOC'], input_tglf['DRMAJDX_LOC'], input_tglf['KAPPA_LOC'],
                input_tglf['P_PRIME_LOC'], input_tglf['Q_LOC'], input_tglf['Q_PRIME_LOC'], input_tglf['RLNS_1'], input_tglf['RLNS_2'], input_tglf['RLNS_3'],
                input_tglf['RLTS_1'], input_tglf['RLTS_2'], input_tglf['RMAJ_LOC'], input_tglf['RMIN_LOC'], input_tglf['S_KAPPA_LOC'], input_tglf['TAUS_2'],
                 input_tglf['VEXB_SHEAR'], input_tglf['XNUE'], input_tglf['ZEFF']]]*2)

    if(tglf_model=='nn'):
        model='tglfnn/models/nn_SAT0_mb_1024_abs_reg_common_stair2x2x6.pb'
        with btf_connect(path=model) as tf:
            qturb = sum(tf.run(nn_inputs)[0,0:2])*q_gb[gridpoint]
    else:
        #OMFIT['TGLF_GACODE']['input.tglf']=OMFITgacode(input_tglf)
        OMFIT['TGLF_GACODE']['SCRIPTS']['runTGLF'].run()
        qturb = sum(OMFIT['TGLF_GACODE']['FILES']['gbflux']['data']['Q/Q_GB'])*q_gb[gridpoint]

    return qturb

# Function builds input.tglf file and runs EPED-NN
def run_epednn(inputs,outputs,betan):
    eped_inputs = {}
    eped_inputs['r']=inputs['R']
    eped_inputs['ip']=inputs['Ip']
    eped_inputs['a'] = outputs['a']
    eped_inputs['bt'] = inputs['Bt']
    eped_inputs['kappa'] = inputs['kappa']
    eped_inputs['delta'] = inputs['delta']
    eped_inputs['neped'] = 10.0*outputs['ne_ped']#convert 10^20 m^-3 to 10^19 m^-3
    eped_inputs['zeffped'] = inputs['zeff']
    eped_inputs['betan'] = betan
    mixing='dt'
    if mixing =='dt':
        eped_inputs['m']=2.5
    else:
        eped_inputs['m']=2.0

    model = 'eped1nn/models/EPED_mb_128_pow_norm_common_30x10.pb'
    nn_inputs=np.atleast_2d([[eped_inputs['a'], eped_inputs['betan'], eped_inputs['bt'], eped_inputs['delta'], eped_inputs['ip'],
               eped_inputs['kappa'], eped_inputs['m'], eped_inputs['neped'], eped_inputs['r'], eped_inputs['zeffped']]])
    with btf_connect(path=model) as tf:
        eped_run= tf.run(nn_inputs)
        pped = eped_run[0,4]
        wped = eped_run[0,-3]

    return wped,pped

# Function for H-mode profiles based on EPED  H-mode profile script
def Hmode_profiles(tedgekEV=.04,tpedkEV = .8,tcorekEV = 4.,rgrid= 101,nedge14 = 0.3,nped14 = .4,
    ncore14 = .55,nexpin=1.1,nexpout=1.1,texpin=2.2,texpout=2.2,widthp = 0.03):
    xphalf = 1.0-widthp
    xped=xphalf-widthp

    pconst=1.-np.tanh((1.-xphalf)/widthp)
    a_n=2.*(nped14-nedge14)/(1.+np.tanh(1.)-pconst)
    a_t=2.*(tpedkEV-tedgekEV)/(1.+np.tanh(1.)-pconst)

    ncoretanh=0.5*a_n*(1.-np.tanh(-xphalf/widthp)-pconst)+nedge14

    # core temperature
    tcoretanh=0.5*a_t*(1.-np.tanh(-xphalf/widthp)-pconst)+tedgekEV

    xpsi = np.linspace(0,1,rgrid)
    ones = np.ones(rgrid)

    nval = 0.5*a_n*(1.-np.tanh((xpsi-xphalf)/widthp)-pconst)+nedge14*ones
    tval = 0.5*a_t*(1.-np.tanh((xpsi-xphalf)/widthp)-pconst)+tedgekEV*ones

    xtoped=xpsi/xped
    for i in range(0,rgrid):
        if(xtoped[i]**nexpin<1.):
            nval[i] = nval[i]+(ncore14-ncoretanh)*(1.-xtoped[i]**nexpin)**nexpout
        if(xtoped[i]**texpin<1.):
            tval[i] = tval[i]+(tcorekEV-tcoretanh)*(1.-xtoped[i]**texpin)**texpout

    return nval,tval

# Solve for fusion and aux power and fluxes
def get_powerbalanceQ(inputs,outputs,ni,Ti,ne,Te,gridpoint):
    # Plasma estimated at ellipsoid torus for volume contribution (triangularity is small correction)
    kappa1=inputs['kappa']
    rtmp = np.linspace(0,outputs['a'],len(Ti))
    kappa0 = 0.5*(1.+kappa1)
    kappa = kappa0*np.ones(101)+(kappa1-kappa0)*rtmp**2/outputs['a']**2
    delta = np.linspace(0,inputs['delta'],101)
    length_elipse = np.pi*(3.*(1.+kappa)-np.sqrt((3.+kappa)*(1.+3.*kappa)))
    length_elipse*= np.ones(101)-0.14*delta

    # Fusion power model from Freidberg et al. PoP 2015
    k0 = -60.4593
    k1 = 6.1371
    k2 = -0.8609
    k3 = 0.0356
    k4 = -0.0045

    integrand = 0.2*17.6*1.6e-19*ni*ni*rtmp/4.0
    integrand *=  np.exp(k0+k1*np.log(Ti)+k2*np.log(Ti)**2+k3*np.log(Ti)**3+k4*np.log(Ti)**4)

    # Fluxes (Q) calculated in W/cm^2 
    # This assumes that alpha heating profile is same as alpha birth profile
    alphaQ =   scipy.integrate.cumtrapz(integrand,x=rtmp)/rtmp[1:]/1.0e-2
    alphaPower = alphaQ[-1]*rtmp[-1]*length_elipse[-1]*(2*np.pi*inputs['R'])*1.0e-2
    
    powerDensityBrem = 1.0e-6*ne*ne*inputs['zeff']*np.sqrt(Te*1000.)/(7.69e18)**2
    bremQ = -scipy.integrate.cumtrapz(powerDensityBrem*rtmp,x=rtmp)/rtmp[1:]/1.0e-2
    bremPower = bremQ[-1]*rtmp[-1]*length_elipse[-1]*(2*np.pi*inputs['R'])*1.0e-2

    
    # Auxiliary power estimated as exp function based on DIII-D beams
    auxHeatingProfile =scipy.integrate.cumtrapz(rtmp*np.exp(-4.5*rtmp),x=rtmp)
    auxHeatingProfile*=1.0/auxHeatingProfile[-1]
    auxQ = inputs['Paux']/(2.*np.pi*inputs['R'])/length_elipse[gridpoint]/rtmp[gridpoint]/1.0e-2*auxHeatingProfile[gridpoint-1]
    return auxQ+alphaQ[gridpoint-1]+bremQ[gridpoint-1],alphaPower,bremPower

def get_flux_error(inputs,outputs,tglf_model,T0,betaN):
    # Run EPED for pedestal
    wped,pped=run_epednn(inputs,outputs,betan=betaN)
    outputs['teped'] =pped/outputs['ne_ped']*1.0e4/1.602e2/(1.+inputs['Tratio'])
    # Get T,n profiles with new pedestal
    rgrid=101
    dT0=inputs['tempShape']
    nval,tval=Hmode_profiles(tcorekEV=T0,texpin=dT0,texpout=dT0,tpedkEV=outputs['teped'],nped14 = outputs['ne_ped'],nedge14 = 0.5*outputs['ne_ped'],ncore14 = outputs['ne_core'],widthp=wped*0.5,rgrid=rgrid)

    # Get fluxes and powers fusion fusion and Aux
    mixing='dt'
    if mixing == 'dt':
        zimp1 =  2
        zimp2 =  10
        impurityFraction2 = (inputs['zeff']-1.)/(5*zimp1**2+zimp2**2-5*zimp1-zimp2)# setting He impurity density to 5 times Neon
        impurityFraction1 = 5.0*impurityFraction2
        ni=nval*(1.-zimp2*impurityFraction2-zimp2*impurityFraction2)*1.0e20
    else:
        zimp = 6.0
        impurityFraction=(inputs['zeff']-1.)/zimp/(zimp-1.0)
        ni = nval*(1.-zimp*impurityFraction)*1.0e20

    exp_Q,Palpha,Pbrem = get_powerbalanceQ(inputs,outputs,ni=ni,Ti=inputs['Tratio']*tval,ne=nval*1.0e20,Te=tval,gridpoint=60)

    # Get betap for Shafranov shift estimate
    rtmp = np.linspace(0,1,rgrid)
    beta = np.pi*8.0e-7*(1.0+inputs['Tratio'])*np.trapz(rtmp*nval*1.6e1*tval*1.0e3)/np.trapz(rtmp)/inputs['Bt']/inputs['Bt']
    betaN = beta*outputs['a']*inputs['Bt']/inputs['Ip']*1.0e2
    betap = 25.*0.5*(1.+inputs['kappa']**2)*(betaN*0.01)**2/beta

    #Run TGLF
    tglf_Q = run_tglf(inputs,outputs,tglf_model,nval*1.0e14,tval*1.e3,betap=betap,gridpoint=60)
    error_rho60 = (exp_Q-tglf_Q)/exp_Q

    # Write output data
    Pheat = Palpha+inputs['Paux']
    m=2.0
    tau98 = (0.0562*inputs['Ip']**0.93*inputs['Bt']**0.15*Pheat**(-0.69)*(10.0*outputs['ne_ave'])**(0.41)*
             m**0.19*inputs['R']**1.97*(inputs['aspect'])**(-0.58)*inputs['kappa']**0.78)
    storedEnergy = 3.0*beta*(2*np.pi*inputs['R'])*(np.pi*inputs['kappa']*outputs['a']*outputs['a'])*inputs['Bt']*inputs['Bt']/(1.6*np.pi)

    outputs['Qplasma']=5.0*Palpha/inputs['Paux']
    outputs['Pbrem'] = Pbrem
    outputs['Pfus'] = 5.0*Palpha
    outputs['taue'] = storedEnergy/Pheat
    outputs['tau98'] = tau98
    outputs['H98'] = storedEnergy/Pheat/tau98
    outputs['fbs'] = 0.7/inputs['aspect']**0.5*betap
    outputs['T0'] = T0
    outputs['betaN'] = betaN
    outputs['tval'] = tval
    outputs['nval'] = nval

    return error_rho60,betaN


def tau_tenn(inputs,tglf_model='full'):

    outputs = {}
    outputs['a'] = inputs['R']/inputs['aspect']

    outputs['n_gw'] = inputs['Ip']/(np.pi*outputs['a']*outputs['a'])
    try:
        outputs['ne_ped'] = inputs['ne_ped']
    except:
        outputs['ne_ped'] = inputs['fgw_ped']*outputs['n_gw']

    outputs['ne_core'] = inputs['density_peaking']*outputs['ne_ped']

    rgrid=101
    rtmp=np.linspace(0,1,rgrid)
    nval,tval=Hmode_profiles(nped14 = outputs['ne_ped'],ncore14 = outputs['ne_core'],nedge14 = 0.5*outputs['ne_ped'],rgrid=rgrid)
    outputs['ne_ave'] = np.trapz(rtmp*nval)/np.trapz(rtmp)

    betaN=3.0
    #q95 estimate from L. Larent "Improvement of Tokamak concept" IAEA collection
    S = 0.5*(1.+inputs['kappa']**2*(1+0.3*inputs['delta']))*(1+1.5*(outputs['a']/inputs['R'])**2)
    outputs['q95'] = 5.*outputs['a']**2*inputs['Bt']/inputs['R']/inputs['Ip']*S
    Tguess = 0.01*betaN*inputs['Bt']*inputs['Ip']/outputs['a']/outputs['ne_core']/(2.0*1.6e1*4.0*np.pi*1.0e-4)
    Tak=1.0*Tguess
    Tbk=1.5*Tguess

    fak,betaN = get_flux_error(inputs=inputs,outputs=outputs,tglf_model=tglf_model,T0=Tak,betaN=betaN)
    fbk,betaN = get_flux_error(inputs=inputs,outputs=outputs,tglf_model=tglf_model,T0=Tbk,betaN=betaN)

    outputs['error'] =[]
    for i in range(0,20):
        Tck = max(outputs['teped']*1.05,Tbk-fbk*(Tbk-Tak)/(fbk-fak))
        fck,betaN = get_flux_error(inputs=inputs,outputs=outputs,tglf_model=tglf_model,T0=Tck,betaN=betaN)
        outputs['error'].append(fck)

        if abs(fck)<inputs['error'] and i>0:break
        if(np.sign(fak)==np.sign(fck)):
            Tak = copy.deepcopy(Tck)
            fak = copy.deepcopy(fck)
        else:
            Tbk = copy.deepcopy(Tck)
            fbk = copy.deepcopy(fck)
    return outputs

#Lines for running in omfit
#outputs=tau_tenn(inputs=inputs,tglf_model=tglf_model)
#OMFIT['outputs']=outputs
