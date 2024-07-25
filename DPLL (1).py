import numpy as np
import control
import matplotlib.pyplot as plt
pi=np.pi
def Bode_10dB(H,block_Name,Hz,TF):
    mag, phase, w = control.bode_plot(H,plot=False,omega_limits=(100,10**9))
    mag=mag*TF
    magnitude = 10 * np.log10(mag)
    ww=w
    if Hz==1:
        ww=w/(2*pi)
    # Plot the Bode plot
    plt.semilogx(ww,magnitude,label=block_Name)
    plt.xlabel('Frequency (rad/sec)')
    if Hz==1:
        plt.xlabel('Frequency (Hz)')
    plt.ylabel('Magnitude (dBc/Hz)')
    plt.legend(fontsize="7", loc='lower center')
    plt.grid(True)
    return mag, phase, w
def Bode2_10dB(mag,w,block_Name,Hz):
    # Calculate magnitude in dB using 10log(x)
    magnitude = 10 * np.log10(mag)
    ww=w
    if Hz==1:
        ww=w/(2*pi)
    # Plot the Bode plot
    plt.semilogx(ww,magnitude,label=block_Name)
    plt.xlabel('Frequency (rad/sec)')
    plt.legend(fontsize="7", loc='lower center')
    if Hz==1:
        plt.xlabel('Frequency (Hz)')
    plt.ylabel('Magnitude (dBc/Hz)')
    plt.grid(True)


s=control.tf('s')
Fref=(150*10**6)
Tref=1/Fref
TDC_res=0.1*10**-12
Fc=1.5*10**6 #assume for now
PM=pi/2.6 #phase margin = 69 degrees
P=4
N=64.5
N1=18
N2=5
Kvco=55*10**6 #in Hz/V
Kdco=(1.2/2**N1)*Kvco*2*pi #in rad/s/LSB
Fout=Fref*N/P
Fvco=Fref*N
DTC_res=1/((2**10)*Fvco)

## CALCULATING ALPHA AND BETA
B_A=(2*pi*Fc*Tref)/np.tan(PM)
Beta=(2*pi*TDC_res*((2*pi*Fc)**2)*N)/(np.sqrt((np.tan(PM)**2)+1)*Kdco)
alpha=Beta*(1/B_A)
print('Beta = ',Beta)
print('Alpha = ',alpha)

## LOOP GAIN EQUATION
Wz=(Beta*Fref)/alpha
b=(Kdco*Beta)/(2*pi*TDC_res*N)
LG=b*(1+s/Wz)*(1/(s*s))
plt.figure()
control.bode_plot(LG, dB=True,Hz=1,omega_limits=(100,10**9),plot=1)
plt.title('Open loop transfer function')
x=control.margin(LG)
print('Phase Margin = ',x[1])
print('Bandwidth = ',x[3]/(2*pi*10**6),' Mhz')


## CLOSED LOOP GAIN EGUATION
CL=LG*N/(1+LG)
plt.figure()
control.bode_plot(CL, dB=True,Hz=1,omega_limits=(100,10**9),plot=1)
plt.title('Closed loop transfer function')


###########################.noise transfer function of each block.####################################

#Refrence
Ref=CL/P
plt.figure()
mag_Ref, phase, w =control.bode_plot(Ref,dB=1, plot=1,omega_limits=(100,10**9))
plt.title('Transfer function of Refrence')

#N-Divider
N_div=CL/P
plt.figure()
mag_N_div, phase, w = control.bode_plot(N_div,dB=1, plot=1,omega_limits=(100,10**9))
plt.title('Transfer function of N-Divider')

#VCO
VCO=(1/(1+LG))*(1/P)
plt.figure()
mag_VCO, phase, w = control.bode_plot(VCO,dB=1, plot=1,omega_limits=(100,10**9))
plt.title('Transfer function of VCO')

#TDC
TDC=(CL*(2*pi*TDC_res/Tref))/P
plt.figure()
mag_TDC, phase, w = control.bode_plot(TDC,dB=1, plot=1,omega_limits=(100,10**9))
plt.title('Transfer function of TDC')

#DAC
DAC=(1/(1+LG))*(Kvco/s)*(1/P)
plt.figure()
mag_DAC, phase, w = control.bode_plot(DAC,dB=1, plot=1,omega_limits=(100,10**9))
plt.title('Transfer function of DCO QN')

#DTC
DTC=CL/P
plt.figure()
mag_DTC, phase, w = control.bode_plot(DTC,dB=1, plot=1,omega_limits=(100,10**9))
plt.title('Transfer function of DTC')


###############################.modeled phase noise of each block.#####################################
plt.figure()

#DSM_Noise
Sg=1/(12*Fref)
num=(4*(pi**2))*(s**2)
den=(N**2)*(Fref**2)
DSM=Sg*num/den
Bode_10dB(DSM,'DSM noise spectrum',1,1)

#DAC_QN
delta_f=(1.2*Kvco)/(2**N2)
DAC_QN=((4*pi**2)*(s**2)*delta_f**2)/((Fref**5)*12)
Bode_10dB(DAC_QN,'DCO QN noise spectrum',1,1)


#TDC_QN
TDC_QN=(((2*pi*TDC_res)/Tref)**2)*(s**0)/(12*Fref)
Bode_10dB(TDC_QN,'TDC QN noise spectrum',1,1)

#DTC_QN
DTC_QN=((2*pi*DTC_res*Fvco)**2)*(Tref/12)*(s**2)*(Tref**2)
Bode_10dB(DTC_QN,'DTC QN noise spectrum',1,1)

#Buffer & Dividers
num1 =(1+(s/(2*pi*2*(10**6))))
den1 =s
H =num1/den1
H=H*1.22*10**-9
#print ('H(s) =', H)
Bode_10dB(H,'Buffer & Dividers noise spectrum',1,1)

#VCO
w1=2*pi*10**6
w2=pi*2*3*10**7
num2=(1+s/w1)*(1+s/w2)*(1+s/w2)
den2=(s*s*s)
g=(num2/den2)*10**9.05
Bode_10dB(g,'VCO noise spectrum',1,1)

#Crystal_Refrence
w11=2*pi*10**3
w22=2*pi*30*10**3
num3=(1+s/w11)*(1+s/w22)*(1+s/w22)
den3=(s*s*s)
G=(num3/den3)*10**-1.85
Bode_10dB(G,'Crystal_Refrence noise spectrum',1,1)

#DTC and TDC noise
Conv_noise = 2.50732E-9 * ((s/(2*np.pi*2E6)) + 1)/s
Bode_10dB(Conv_noise,'Data Converters Noise spectrum',1,1)



#####################################.close noise VS VCO noise.#########################################
plt.figure()

#P-divider phase noise curve
mag1_Pdiv,phase1,w=Bode_10dB(H,'P-divider phase noise curve',Hz=1,TF=N**2)

#N-divider phase noise curve
mag1_Ndiv,phase1,w=Bode_10dB(H,'N-divider phase noise curve',Hz=1,TF=N**2)

#buffer phase noise curve
mag1_buff,phase1,w=Bode_10dB(H,'Buffer phase noise curve',Hz=1,TF=N**2)

#Refrence phase noise curve
mag1_Ref,phase1,w=Bode_10dB(G,'Crystal-Refrence phase noise curve',Hz=1,TF=N**2)

#DCO QN curve
mag1_DAC_QN,phase1,w=Bode_10dB(DAC_QN,'DCO QN curve',Hz=1,TF=N**2)

#TDC QN curve
mag1_TDC_QN,phase1,w=Bode_10dB(TDC_QN,'TDC QN noise curve',Hz=1,TF=N**2)

#DTC QN curve
mag1_DTC_QN,phase1,w=Bode_10dB(DTC_QN,'DTC QN noise curve',Hz=1,TF=N**2)

#VCO phase noise curve
mag1_VCO,phase,W=Bode_10dB(g,'VCO phase noise curve',Hz=1,TF=1)

#TDC and DTC phase noise curve
mag1_Conv_N,phase1,w=Bode_10dB(Conv_noise,'TDC and DTC phase noise curve',Hz=1,TF=N**2)

#DSM phase noise curve
mag1_DSM,phase1,w=Bode_10dB(DSM,'DSM phase noise curve',Hz=1,TF=N**2)


Mag1_close_noise=2*mag1_Conv_N+mag1_Ref+mag1_Ndiv+mag1_buff+mag1_Pdiv+mag1_TDC_QN+mag1_DAC_QN+mag1_DTC_QN#+mag1_DSM

plt.figure()
plt.title('Wu Optimized at intersection ')
Bode2_10dB(Mag1_close_noise,W,'Close noise',Hz=1)
Bode2_10dB(mag1_VCO,W,'VCO',Hz=1)

##################################.Total Noise at output.############################################
plt.figure()
plt.title('Phase Noise of each block at output')

#P-divider phase noise curve
mag2_Pdiv,phase1,w=Bode_10dB(H,'P-divider phase noise curve',Hz=1,TF=1)

#N-divider phase noise curve
mag2_Ndiv,phase1,w=Bode_10dB(H,'N-divider phase noise curve',Hz=1,TF=mag_N_div**2)

#buffer phase noise curve
mag2_buff,phase1,w=Bode_10dB(H,'Buffer phase noise curve',Hz=1,TF=1)

#Refrence phase noise curve
mag2_Ref,phase1,w=Bode_10dB(G,'Crystal-Refrence phase noise curve',Hz=1,TF=mag_Ref**2)

#DCO QN curve
mag2_DAC_QN,phase1,w=Bode_10dB(DAC_QN,'DCO QN curve',Hz=1,TF=mag_DAC**2)

#TDC QN curve
mag2_TDC_QN,phase1,w=Bode_10dB(TDC_QN,'TDC QN noise curve',Hz=1,TF=mag_DTC**2)

#DTC QN curve
mag2_DTC_QN,phase1,w=Bode_10dB(DTC_QN,'DTC QN noise curve',Hz=1,TF=mag_DTC**2)

#VCO phase noise curve
mag2_VCO,phase,W=Bode_10dB(g,'VCO phase noise curve',Hz=1,TF=mag_VCO**2)

#DTC phase noise curve
mag2_Conv_N,phase1,w=Bode_10dB(Conv_noise,'DTC phase noise curve',Hz=1,TF=mag_DTC**2)

#TDC phase noise curve
mag22_Conv_N,phase1,w=Bode_10dB(Conv_noise,'TDC phase noise curve',Hz=1,TF=mag_TDC**2)

#DSM phase noise curve
#mag2_DSM,phase1,w=Bode_10dB(DSM,'DSM phase noise curve',Hz=1,TF=mag_N_div**2)


Mag2_close_noise=mag22_Conv_N+mag2_Conv_N+mag2_Ref+mag2_Ndiv+mag2_buff+mag2_TDC_QN+mag2_Pdiv+mag2_DAC_QN+mag2_DTC_QN#+mag2_DSM
Mag_total=Mag2_close_noise+mag2_VCO
plt.figure()
plt.title('Total noise at output')
Bode2_10dB(Mag_total,W,'output reffered noise',Hz=1)

#plt.show()

filterr=np.zeros(1000)
filterr[411:872]=1
W_Int=W*filterr/(2*pi)
Mag_Int=Mag_total*filterr

RMS_sq=np.trapz(Mag_Int,W_Int)
RMS_Jitter=np.sqrt(2*(RMS_sq))/(2*pi*2.42*10**9)


#############################################.percentage of each block.########################################
#p-divider
P_divider=mag2_Pdiv*filterr
pdiv=np.trapz(P_divider,W_Int)

#N-divider
N_divider=mag2_Ndiv*filterr
ndiv=np.trapz(N_divider,W_Int)

#Crystal-Refrence
refrence=mag2_Ref*filterr
ref=np.trapz(refrence,W_Int)

#VCO
VCO_out=mag2_VCO*filterr
vco=np.trapz(VCO_out,W_Int)

#DAC
DACC=mag2_DAC_QN*filterr
dac=np.trapz(DACC,W_Int)

#TDC
TDCC=(mag2_TDC_QN+mag22_Conv_N)*filterr
tdc=np.trapz(TDCC,W_Int)

#DTC
DTCC=(mag2_DTC_QN+ mag2_Conv_N)*filterr
dtc=np.trapz(DTCC,W_Int)

#Buffer
Buffer=mag2_buff*filterr
buff=np.trapz(Buffer,W_Int)

print('\n')
print('Total output integrated RMS jitter is ', RMS_Jitter)
print('\n')
print('P-divider percentage in the total output integrated RMS jitter is ',(pdiv/RMS_sq)*100,'%')
print('N-divider percentage in the total output integrated RMS jitter is ',(ndiv/RMS_sq)*100,'%')
print('Buffer percentage in the total output integrated RMS jitter is ',(buff/RMS_sq)*100,'%')
print('DAC percentage in the total output integrated RMS jitter is ',(dac/RMS_sq)*100,'%')
print('TDC percentage in the total output integrated RMS jitter is ',(tdc/RMS_sq)*100,'%')
print('DTC percentage in the total output integrated RMS jitter is ',(dtc/RMS_sq)*100,'%')
print('VCO percentage in the total output integrated RMS jitter is ',(vco/RMS_sq)*100,'%')
print('Refrence percentage in the total output integrated RMS jitter is ',(ref/RMS_sq)*100,'%')