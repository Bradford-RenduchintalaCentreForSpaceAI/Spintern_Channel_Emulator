from scipy import signal


class Filters():
    from scipy import signal
    
    
    def Low_pass(self,s,cutoff,N,fs):
        sos = signal.butter(N, cutoff, fs = fs,btype ='low', output = 'sos')
        return signal.sosfilt(sos, s)

        
        
    
    
    
