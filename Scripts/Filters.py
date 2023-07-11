from scipy import signal


def low_pass(s,cutoff,N,fs):
    sos = signal.butter(N, cutoff, fs = fs,btype ='low', output = 'sos')

    return signal.sosfilt(sos, s)
