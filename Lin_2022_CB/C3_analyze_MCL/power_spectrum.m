
function [freq, dataPS] = power_spectrum(data, frame_interval)
    
    N = length(data);
    dataFFT = fft(data);
    dataFFT = dataFFT(1:floor(N/2)+1);
        
    sample_freq = 1/frame_interval; 
    freq = (0: (sample_freq/N): (sample_freq/2) )'; 
    
    % Power spectrum
    dataPS = (1/(sample_freq*N))* abs(dataFFT).^2;    
    
    % Adjust the power spectrum by 2, except for the 0 and highest frequency
    dataPS(2:end-1) = 2*dataPS(2:end-1);    

end