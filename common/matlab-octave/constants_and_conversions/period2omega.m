function omega = period2omega (T)
% convert period in seconds to frequency in radians/s
    
    omega = 2 .* pi .* (1./T);

end