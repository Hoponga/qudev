# notes 
# 1. convert out of log mag to lin mag & do izt with complex data 
# - Kind of done, but doesn't work as intended I think 
# 2. look at window functions for smoothing before -- hamming window
# convert time domain data back to log mag  

# 3. increase SNR -- try demodulation on the real data -- multiply by frequency of signal 
# or moving average
# 
# scikit-rf 

# read data 

# 
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import get_window 

import pandas as pd
#datafile = 'Trace02freq.csv'
datafile = 'Trace03.csv' # this one has real/imag values 
trace_data = pd.read_csv(datafile)
print(trace_data)
X_COL = "Stimulus(Hz)"
SIGNAL_COL = "S21(dB)"

REAL_COL = "S21(Real)"
IMAG_COL = "S21(Imag)"
DISPLAY_LOG_MAG = True 
trace_frequencies = trace_data[X_COL].to_numpy()
delta_F = (trace_frequencies[-1] - trace_frequencies[0])/len(trace_frequencies)
trace_signal_real = trace_data[REAL_COL].to_numpy() 
trace_signal_imag = trace_data[IMAG_COL].to_numpy()
trace_signal = trace_signal_real + 1j*trace_signal_imag
trace_signal /= delta_F


print(trace_data.columns)
#trace_signal = trace_data[SIGNAL_COL].to_numpy()

#plt.plot(trace_signal)

# apply windowing function to smooth signal in frequency domain 
def windowing_function(signal, function = "hamming", timesteps = 1000): 
    m = len(signal)
    t = np.arange(m)
    w = get_window(function, m)
    return signal * w 



def np_ifft(signal): 
    itx = np.fft.ifft(signal) 
    return itx 

# to time domain, x(n) = 1/N * sum_0^{n-1} x_f[k] * exp((2j*pi/N)*kn)
# algorithm is exact same as FFT just with different factors
def ifft(signal): 
    
    n = len(signal) 
    print(len(signal))
    
    if n == 1:
        return signal
    else:
        
        signal_even = ifft(signal[::2])
        signal_odd = ifft(signal[1::2])
        factor = \
          np.exp(2j*np.pi*np.arange(n)/ n)
        
        
        signal_ret = np.concatenate(\
            [signal_even+factor[:int(n/2)]*signal_odd,
             signal_even+factor[int(n/2):]*signal_odd])
       
        return 1/n * signal_ret
    

def demodulate(signal): 
    freq_response = np.fft.fft(signal) 
    freqs = np.fft.fftfreq(signal)
    for c, freq in zip(freq_response, freqs): 
        if c: 
            print(f"{c} * exp(2 pi t {freq})") 




# # for reference, look at gnu octave implementation 
# def czt(x, m=None, w=None, a=None):

#     n = len(x)
#     if m is None: m = n
#     if w is None: w = np.exp(-2j * np.pi / m)
#     if a is None: a = 1

#     # In our chirp z-transform, we have values k - n range from k = 1 (1-n) to k = 
#     chirp = w ** (np.arange(1 - n, max(m, n)) ** 2 / 2.0)
#     N2 = int(2 ** ceil(log2(m + n - 1)))  # next power of 2

#     # assume dft -> convolution of sequences an & bn 
#     #an = append(x * a ** -arange(n) * chirp[n-1:(n+n-1)])

#     xp = append(x * a ** -arange(n) * chirp[n - 1 : n + n - 1], zeros(N2 - n))
#     ichirpp = append(1 / chirp[: m + n - 1], zeros(N2 - (m + n - 1)))


#     r = ifft(fft(xp) * fft(ichirpp))
#     return r[n - 1 : m + n - 1] * chirp[n - 1 : m + n - 1]


def iczt(x, m = None, w = None, a = None): 
    n = len(x)
    if m is None: m = n
    if w is None: w = np.exp(-2j * np.pi / m)
    if a is None: a = 1

    # In our chirp z-transform, we have values k - n range from k = 1 (1-n) to k = 
    chirp = w ** (-np.arange(1 - n, max(m, n)) ** 2 / 2.0)
    N2 = int(2 ** np.ceil(np.log2(m + n - 1)))  # next power of 2

    # assume dft -> convolution of sequences an & bn 
    #an = append(x * a ** -arange(n) * chirp[n-1:(n+n-1)])
    #xp = x * a ** arange(n) * chirp[n-1 : n + n - 1]
    #ichirpp = 1 / chirp[: m + n - 1]

    xp = np.append(x * a ** np.arange(n) * chirp[n - 1 : n + n - 1], np.zeros(N2 - n))
    ichirpp = np.append(1 / chirp[: m + n - 1], np.zeros(N2 - (m + n - 1)))

    xp_transformed = np.fft.fft(xp)
    ichirpp_transformed = np.fft.fft(ichirpp)
    result = np.fft.ifft(xp_transformed * ichirpp_transformed) 
    return result[n-1 : m + n-1] * chirp[n-1 : n + n-1]
    
    r = np.fft.ifft(np.fft.fft(xp) * np.fft.fft(ichirpp))
    return r[n - 1 : m + n - 1] * chirp[n - 1 : m + n - 1]

# to freq domain, x(w) = sum_0^{n-1} x(k) * exp[(-2j*pi/N)*kn]


def fft(signal):
    
    n = len(signal) 
    
    if n == 1:
        return signal
    else:
        signal_even = fft(signal[::2])
        signal_odd = fft(signal[1::2])
        assert len(signal_even) == len(signal_odd)
        factor = \
          np.exp(-2j*np.pi*np.arange(n)/ n)
        
        
        signal_ret = np.concatenate(\
            [signal_even+factor[:int(n/2)]*signal_odd,
             signal_even+factor[int(n/2):]*signal_odd])
        return signal_ret


#plt.plot(np_ifft(trace_signal))
print(len(trace_signal))
windowed_signal = windowing_function(trace_signal)
time_domain_signal = iczt(windowed_signal)
time_domain_signal = np.abs(time_domain_signal)


# demodulation stuff 



DISPLAY_LOG_MAG = True
if DISPLAY_LOG_MAG: 
    plt.ylabel("log mag amplitude")
    time_domain_signal = 20*np.log10(time_domain_signal)

else: 
    plt.ylabel("lin mag amplitude")


print(time_domain_signal)
plt.plot(time_domain_signal)
plt.title("time domain signal")
plt.xlabel("time (ns)")
plt.ylabel("Log mag amplitude")

plt.savefig("time_domain.png")
plt.show()


# plt.title("Sine wave plotted using inverse Fourier transform");
# plt.xlabel('Time')
# plt.ylabel('Amplitude')
# plt.grid(True)
# plt.show();