function [f, power] = odeoutputfft(T, Y)
% calculates an fft from a set of data sampled at nonuniformally spaced
% sample points
%
% Syntax
%
% [f, power] = odeoutputfft(T, Y)
%
% Input
%
% T - the sample points of the data in Y (usually a time series)
%
% Y - a set of data points sampled at the positions in T 
%
% 

% Copyright Richard Crozier 2010

% fs = 100;                          % Sample frequency (Hz)
% t = 0:1/fs:10-1/fs;                % 10 sec sample
% x = (1.3)*sin(2*pi*15*t) ...       % 15 Hz component
%   + (1.7)*sin(2*pi*40*(t-2)) ...   % 40 Hz component
%   + (2.5)*randn(size(t));          % Gaussian noise;
% % 
% m = length(x);          % Window length
% n = pow2(nextpow2(m));  % Transform length
% y = fft(x,n);           % DFT
% f = (0:n-1)*(fs/n);     % Frequency range
% power = y.*conj(y)/n;   % Power of the DFT

    tdiff = mean(T(2:end) - T(1:end-1))/2;
    
    fs = 1/tdiff;
    
    x = interp1(T, Y, 0:tdiff:max(T));
    
    m = length(x);
    
    n = pow2(nextpow2(m/4)); % Next power of 2 from length of y
    
    %n = 2048;
    
    y = fft(x,n);
    
    y0 = fftshift(y);          % Rearrange y values
    f = (-n/2:n/2-1)*(fs/n);  % 0-centered frequency range
    power = y0.*conj(y0)/n;   % 0-centered power
% 
%     f = (0:n-1)*(fs/n);
%     
%     power = y.*conj(y)/n;
%     
end