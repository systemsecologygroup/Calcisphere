# energetics of annoplankton calcification
# this version for sage to compute steady state

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rc('text', usetex=True)

# initialize
fmax=100
dt=1.0

Ca_0=np.zeros(fmax)
EATP=np.zeros(fmax)
EATP_v=np.zeros(fmax)
EATP_m=np.zeros(fmax)

# cytoplasm
# r_cyt=10.0e-6 # diameter of calcispheres is around 20 micrometers
r_cyt=2.5e-6 # diameter of E. huxleyi is around 5 micrometers (Harvey et al 2015)
A_cyt=4.0*np.pi*r_cyt**2.0 # cell surface in m2

# set parameters
D_Ca=7.93e-6/10000.0 # cm2 s-1 -> m2 s-1 (after Li and Gregory 1974)

# calcification flux
QCa=6.11e-18 # 22 fmol h-1 -> 6.11e-18 mol s-1


Ca_out=10.0   # 10 mM=10e-3 mol L-1 -> 10e-6 mol m-3 could also be set to 20.0e-3 M
Ca_in=0.10e-3 # mM


# minimum permeability must be 22 fmol h-1 to precipitate one full coccolith with 22 fmol in one hour (Holtz et al 2013)

Ca_bd=Ca_out-QCa/(4.0*np.pi*r_cyt*D_Ca)
print 'Ca_bd=', Ca_bd
dCa=Ca_out-Ca_bd
fCa=dCa/Ca_out*100.0
print 'fCa=',fCa

Ca_hi=60.0
print Ca_hi/Ca_in

# if the calcification flux is 6.11e-18 mol s-1
# then the conductance must be at least this value

# ion flux per channel is 1pA, which corresponds to 3.0e6 divalent ions per second! (Tsien 1983)
i=3.0e6 # Ca2+ s-1 channel-1

# channel density up to 30-60 per mum-2 in snail axons (Tsien 1983)
N=1.0e12 # channels m-2

# this reults in:
I_Ca=N*i # 1.8e20 Ca ions m-2 s-1
N_A=6.0221367e23 # ions per mol
J_Ca=I_Ca/N_A # 2.98e-4 mol Ca m-2 s-1 (as maximum flux)
F_Ca=J_Ca*A_cyt # 2.3475e-14 mol s-1
print 'F_Ca=',F_Ca
print 'Q_Ca=',QCa
print 'percent of flux=',QCa/F_Ca*100.0
fV=QCa/F_Ca
print 'fV=',fV

PCa=F_Ca/(Ca_bd-Ca_in) # because PCa is the rate that leads to this flux given the ion gradient
print PCa

Ca_0 = (4.0*Ca_out*D_Ca*np.pi*r_cyt + Ca_in*PCa*fV)/(4.0*D_Ca*np.pi*r_cyt + PCa*fV)
print Ca_0, Ca_bd, Ca_0-Ca_bd # just to check if the equation is correct

ATP_per_Ca=0.5
EATP1 = PCa*(Ca_0-Ca_in)*ATP_per_Ca
EATP2 = PCa*(Ca_out-Ca_in)*ATP_per_Ca
print EATP1, EATP2
dEATP=EATP2-EATP1
print 'dEATP=',dEATP
fATP=dEATP/EATP2*100.0
print 'percent Energy savings=',fATP
print 'percent of flux saving=', (1.0-EATP1/EATP2)*100.0

for i in range(fmax):
    #print i
    f_V=i/100.0
    Ca_01 = (4.0*Ca_out*D_Ca*np.pi*r_cyt + Ca_in*PCa*f_V)/(4.0*D_Ca*np.pi*r_cyt + PCa*f_V)
    EATP[i]  =           PCa*(Ca_01-Ca_in)*ATP_per_Ca*60.0*60.0*24.0*1e12
    EATP_v[i]= f_V*      PCa*(Ca_01-Ca_in)*ATP_per_Ca*60.0*60.0*24.0*1e12
    EATP_m[i]= (1.0-f_V)*PCa*(Ca_01-Ca_in)*ATP_per_Ca*60.0*60.0*24.0*1e12

fV_ax=np.arange(0,1,0.01)

# plotting

plt.figure(3)
plt.plot(fV_ax,EATP_m,color='#4C832C',label='plasma membrane',lw=3)
plt.plot(fV_ax,EATP_v,color='#963932',label='vesicle',lw=3)
plt.plot(fV_ax,EATP,color='#362B67',label='combined',lw=3)
plt.plot(0.0156,PCa*((4.0*Ca_out*D_Ca*np.pi*r_cyt + Ca_in*PCa*0.0156)/(4.0*D_Ca*np.pi*r_cyt + PCa*0.0156)-Ca_in)*ATP_per_Ca*60.0*60.0*24.0*1e12,'o',color='#362B67',ms=10)
plt.plot(0.0156,0.0156*PCa*((4.0*Ca_out*D_Ca*np.pi*r_cyt + Ca_in*PCa*0.0156)/(4.0*D_Ca*np.pi*r_cyt + PCa*0.0156)-Ca_in)*ATP_per_Ca*60.0*60.0*24.0*1e12,'o',color='#963932',ms=10)
plt.plot(0.0156,(1.0-0.0156)*PCa*((4.0*Ca_out*D_Ca*np.pi*r_cyt + Ca_in*PCa*0.0156)/(4.0*D_Ca*np.pi*r_cyt + PCa*0.0156)-Ca_in)*ATP_per_Ca*60.0*60.0*24.0*1e12,'o',color='#4C832C',ms=10)
#plt.axis([0, 1, 0, 2.5e-16])
#plt.plot(tax,Ca_vs,'k',label='ves')
plt.ylabel(r'energy for calcium transport $\mathbf{(pmol ~ATP ~d^{-1})}$',fontsize=16)
plt.xlabel(r'fraction of $\mathbf{Ca^{2+}}$ pumped into vesicle ($\mathbf{f_v}$)',fontsize=16)
#plt.title('with vesicle')
plt.legend(loc='right')

plt.show()
