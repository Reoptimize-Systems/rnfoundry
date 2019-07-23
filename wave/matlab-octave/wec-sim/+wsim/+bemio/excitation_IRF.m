function hydro = excitation_IRF (hydro,t_end,n_t,n_w,w_min,w_max)
% Calculates the normalized excitation impulse response function.
% 
% hydro = wsim.bemio.excitation_IRF (hydro, t_end, n_t, n_w, w_min, w_max)
%
% Input
%
%  hydro - data structure containing the field 'B'
%
%  t_end - calculation range for the IRF, where the IRF is calculated from
%   t = 0 to t_end. If t_end is empty ([]), a value of 100 s is used
%
%  n_t - number of time steps in the IRF. If n_t is empty ([]), a value of 
%   1001 is used
%
%  n_w - number of frequency steps used in the IRF calculation 
%   (hydrodynamic coefficients are interpolated to correspond). If n_w is
%   empty ([]), a value of 1001 is used
%
%  w_min - minimum frequency to use in the IRF calculation. If w_min is
%   empty ([]), the minimum frequency from the BEM data is used
%
%  w_max - maximum frequency to use in the IRF calculation. If w_max is
%   empty ([]), the maximum frequency from the BEM data is used
% 
% Output
%
%  hydro - modified data structure with the fields ra_K, ra_t and ra_w
%   added (or replaced)
%
% See also: 
%


    p = waitbar (0,'Calculating excitation IRFs...');  % Progress bar

    CC = onCleanup (@() close (p));

    % Set defaults if empty
    if isempty(t_end)==1;  t_end = 100;           end
    if isempty(n_t)==1;    n_t = 1001;            end
    if isempty(n_w)==1;    n_w = 1001;            end
    if isempty(w_min)==1;  w_min = min(hydro.w);  end
    if isempty(w_max)==1;  w_max = max(hydro.w);  end

    % Interpolate to the given t and w
    t = linspace (-t_end, t_end, n_t);
    w = linspace (w_min, w_max, n_w);  
    N = sum(hydro.dof) * hydro.Nh;

    % Calculate the impulse response function for excitation
    n = 0;

    for i = 1:sum(hydro.dof)

        for j = 1:hydro.Nh

            ex_re = interp1 (hydro.w, squeeze (hydro.ex_re(i,j,:)), w);

            ex_im = interp1 (hydro.w, squeeze (hydro.ex_im(i,j,:)), w);

            hydro.ex_K(i,j,:) = (1/pi) * trapz (w, ex_re .* cos ( w .* t(:)) - ex_im .* sin (w .* t(:)), 2);
            
            n = n + 1;
            
        end
        
        waitbar(n/N)
        
    end
    
    hydro.ex_t = t;
    hydro.ex_w = w;

end