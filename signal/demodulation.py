import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
import scipy

# pozar microwave engineering 

# do weighted average/smoothing function over demodulated data 

def demodulate(signal, times): 
    rectified = np.abs(signal) 

    freq_response = np.fft.fft(signal)
    timestep = times[1] - times[0]
    freqs = np.fft.fftfreq(len(signal), d = timestep) # calculate the freq values corresponding to freq_response magnitudes
    print(freqs)
    #plt.plot(signal)
    #plt.title("unfiltered signal")
    #plt.show()
    
    print(freqs)
    max_coeff = np.argmax(freq_response)
    max_freqs = freqs[max_coeff]
    print(max_freqs)


    actual_freq = np.abs(max_freqs)
    print(actual_freq)
    actual_freq = 4.37e9
    print(f"Actual determined peak frequency: {actual_freq}")
    #b, a = scipy.signal.butter(10, actual_freq, 'lp', fs = 1/timestep)
    #w, h = scipy.signal.freqs(b, a)

    # Let us try applying a gaussian filter to the signal in the frequency domain
    signal_ft = np.fft.fft(signal)
    #gaussian_window = scipy.signal.gaussian(len(signal_ft)*2, std = 2000)[10000:]
    #signal_filtered_ft = signal_ft * gaussian_window 
    signal_filtered_ft = np.zeros(len(signal_ft))
    signal_filtered = np.fft.ifft(signal_filtered_ft)

    other_filtered = signal_filtered*np.exp(-2j*np.pi*actual_freq*times)
    other = signal*np.exp(-2j*np.pi*actual_freq*times)
    plt.plot(signal, times)
    plt.show()
    plt.plot(other, times)
    plt.show()
    return other, other_filtered
    #plt.plot(w, 20*np.log10(np.abs(h)))
    #plt.show()
    sos = scipy.signal.butter(10, actual_freq, 'lp', output = 'sos', fs = 1/timestep)

    return scipy.signal.sosfilt(sos, rectified)

    freq_response = np.fft.fft(rectified)
    low_passed_response = freq_response*h 
    return np.fft.ifft(low_passed_response)
    freqs = np.fft.fftfreq(len(signal))
    max_coeff = np.argmax(freq_response)
    max_freqs = freqs[max_coeff]
    print(max_freqs)
    max_freqs = 4.35e9

    return signal*np.exp(2j*np.pi*np.arange(len(signal))*max_freqs/len(signal))
    # for c, freq in zip(freq_response, freqs): 
    #     if c: 
    #         print(f"{c} * exp(2 pi t {freq})") 


X_COL = "Stimulus(s)"
SIGNAL_COL = "S21(dB)"

REAL_COL = "S21(Real)"
IMAG_COL = "S21(Imag)"

real_tddata_file = "real_timedomaindata_Trace2.csv"
real_td_data = pd.read_csv(real_tddata_file)
signal = real_td_data[REAL_COL].to_numpy()
times = real_td_data[X_COL].to_numpy()
#plt.plot(demodulate(signal, times))

complex_tddata_file = "ComplexTrace2.csv"
complex_td_data = pd.read_csv(complex_tddata_file)
complex_times = complex_td_data[X_COL].to_numpy()
print(complex_times)
complex_signal_real = complex_td_data[REAL_COL].to_numpy()
complex_signal_imag = complex_td_data[IMAG_COL].to_numpy()

complex_signal = complex_signal_real + 1j*complex_signal_imag
#out_signal = demodulate(complex_signal, complex_times)

#out_signal = 20*np.log(np.abs(out_signal))
#plt.plot(demodulate(out_signal, complex_times)) # complex
#out_1 = demodulate(complex_signal_real, complex_times)

#out_2 = demodulate(complex_signal_imag, complex_times)

out_3, out_3_filtered = demodulate(complex_signal, complex_times)
WINDOW_SIZE = 150
def weighted_avg_smoothing(signal, window = WINDOW_SIZE): 
    queue = []
    smoothed = [] 
    for i in range(window, len(signal)): 
        if len(queue) < window: 
            queue.append(signal[i - window])
        elif len(queue) == window: 
            queue.pop(0)
            queue.append(signal[i - window])
        #print(len(queue))
        
        smoothed.append( sum(queue)/window)
    print(f"{len(signal)}-step smoothing done")

    return smoothed 

def low_pass(signal): 
    pass 

# plt.plot(complex_times[1000:], complex_signal[1000:])
# plt.plot(complex_times[1000:], out_3[1000:])
# plt.title("real part of signal before & after demodulation")

smoothed_out_3 = weighted_avg_smoothing(out_3)
smoothed_out_3 = 20*np.log10(np.abs(smoothed_out_3))
out_3 = 20*np.log10(np.abs(out_3))
out_3_filtered = 20*np.log10(np.abs(out_3_filtered))


#out_complex = 20*np.log10(np.abs(out_1 + 1j*out_2))


plt.plot(complex_times[:-WINDOW_SIZE], out_3[:-WINDOW_SIZE])

plt.plot(complex_times[:-WINDOW_SIZE], smoothed_out_3)
plt.title("filtered real signal")


print("test")
plt.show()
plt.savefig("filtered.png")

