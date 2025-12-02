import numpy as np

dt = 1.0 # set time in ps. 
chi = 100 # parameter for setting the cutoff for reduing noise from FFT, "t_thres" > the response function decay timescale!!!

rdata = np.loadtxt('ANALYSIS_Separations-timeseries.dat')
d_sensor = rdata[:,1]
d_effector = rdata[:,2]

### Calculates cross (auto)-correlation of x, y (x) using NumPy "correlate" function. 
def corr(x, y):
    C = np.correlate(x, y, mode='same')
    N = len(C)
    C = C[N//2:]
    lengths = np.arange(N, N//2, -1)
    C_norm = C / lengths
    #C_norm /= C_norm[0]
    return C_norm

### Padding with zeros once a function crosses zeros for the first time after a time threshold Tmax (and returs the index)
def padZeros(arrayFn, Tmax):
    arrayFn[Tmax:] = 0 # Set values to zeros for indices >= Tmax  
    return  

############################################################################################################################################
### Calculations ... 
############################################################################################################################################
### Self-correlation and functions
AC_s = corr(d_sensor, d_sensor)
AC_e = corr(d_effector, d_effector)

### Time lags 
t_lags = dt * (np.arange(0, AC_s.size, 1))

### Cross-correlation function
CC = (corr(d_sensor, d_effector) + corr(d_effector, d_sensor))/2 

### Self-Response function: JS(t) = -{1/kT}{dAC(t)/dt} 
kTJS_s = -np.gradient(AC_s, t_lags, edge_order=2)
kTJS_e = -np.gradient(AC_e, t_lags, edge_order=2)

### Cross-Response function: JC(t) = -{1/kT}{dCC(t)/dt} 
kTJC = -np.gradient(CC, t_lags, edge_order=2)

### Direct Fourier transform of Slef-Response functions (JSd).
# We need to set a cut-off for the derivative (e.g. 5e4 ps) as the error grows with decreasing statistics which affects the FFT!!!
# Padding the array with zeros once JS crosses zero for the first time after "t_thres" and return the index !!!
t_thres = t_lags.size//chi
cutoff_JS_s = padZeros(kTJS_s, t_thres)  
cutoff_JS_e = padZeros(kTJS_e, t_thres)  

fs = np.fft.rfftfreq(2*t_lags.size, d=dt)
JSd_s = np.fft.rfft(kTJS_s, n=2*t_lags.size)
JSd_e = np.fft.rfft(kTJS_e, n=2*t_lags.size)

### Direct Fourier transform of Cross-Response functions (JCd).
cutoff_JC = padZeros(kTJC, t_thres)  

JCd = np.fft.rfft(kTJC, n=2*t_lags.size)

### Force-transmit functions (TF): TF_e = -F_e/F_s = JC/JS_e 
# First, we have to change the zeros in the denominator to ones. 
JSd_s[JSd_s == 0] = 1
JSd_e[JSd_e == 0] = 1

TF_s2e = JCd / JSd_e # sensor-to-effector force transfer
TF_e2s = JCd / JSd_s # effector-to-sensor force transfer
np.savetxt('ForceTransmitFunction_real-and-imaginary_parts.txt', list(zip(fs, TF_s2e.real, TF_s2e.imag, TF_e2s.real, TF_e2s.imag)), fmt='%.12f', header='<Freq.(1/ps)>  <Re_TF_s2e>  <Im_TF_s2e>  <Re_TF_e2s>  <Im_TF_e2s>')

### Mod of the complex transmit function
# Note that "mod" misses though the linear property, and hence the non-additivity of different transmission modes!!! 
np.savetxt('modForceTransmitFunctions.txt', list(zip(fs, np.abs(TF_s2e), np.abs(TF_e2s))), fmt='%.12f', header='<Freq.(1/ps)>  <|TF_s2e|>  <|TF_e2s|>')

