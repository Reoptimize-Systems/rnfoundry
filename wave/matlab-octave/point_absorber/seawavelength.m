function L = seawavelength(f, h, g)
% seawavelength: finds the wavelength of an ocean wave based on linear wave
% theory
%
% Input
%
% f - matrix of wavelengths for which the wave length will be calculated
%
% h - the mean water level (i.e water depth). This can be either a scalor
%     value, or a matrix of values of the same size as f
%
% g - (optional) either a scalor or matrix of the same size as f containing the
%      acceleration due to gravity, if omitted 9.81 m / s^2 is used
%
% Output
% 
% L - matrix of wave lengths of the same size as f
%

    % in 'Lecture notes for the course in water wave mechanics' by Thomas
    % Lykke Anderson & Peter Frigaard (Aalborg University) page 28 the
    % equation for the wavelength of a wave in linear wave theory is given
    % as
    %
    % L = ( g T^2 / 2 pi ) tanh( 2 pi h / L )                 (eq 3.42)
    % 
    % which must be solved iteratively for L, this is achieved here by
    % finding the zero of the function
    %
    % ( g T^2 / 2 pi ) tanh( 2 pi h / L ) - L
 
    % If not supplied, set g to 9.81 m / s^2
    if nargin < 3
        g = 9.81;
    end
    
    % first we will create an initial guess at the wavelength using the
    % formula from Guo 2002 which is based on logarithmic matching
    beta = 2.4908;
    
    x = h .* (2 .* pi .* f) ./ (sqrt(g .* h));
    
    L = 2 .* pi .* h ./ (x.^2 .* (1 - exp(-x.^beta)).^(1./beta));
    
    % next precalculate some parts of the equation to be solved to save time
    term1 = g .* (1./f).^2 ./ (2 .* pi);
    
    % We preallocate a matrix of ones of the appropriate size so the user
    % can pass in a single value of h or one for each frequency if desired
    term2 = ones(size(term1)) .* 2 .* pi .* h;
    
    % then, for each combination of variables find the zero of the function
    % wavemin given below, using the initial guess calculated previously
    for i = 1:numel(term1)
        
        L(i) = fzero(@(x) wavemin(term1(i), term2(i), x), L(i));
        
    end
    
end


function result = wavemin(term1, term2, L)
    % evaluates the formula
    %
    % ( g T^2 / 2 pi ) tanh( 2 pi h / L ) - L
    %
    % given some precalculated parts of the equation
    result = term1 .* tanh(term2 ./ L) - L;
    
end