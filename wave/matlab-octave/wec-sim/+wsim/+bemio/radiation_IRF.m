function hydro = radiation_IRF (hydro, t_end, n_t, n_w, w_min, w_max)
% Calculates the normalized radiation impulse response function.
% 
% Syntax
%
% hydro = wsim.bemio.radiation_IRF (hydro, t_end, n_t, n_w, w_min, w_max)
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


    p = waitbar(0,'Calculating radiation IRFs...');  % Progress bar

    CC = onCleanup (@() close (p));

    % Set defaults if empty
    if isempty(t_end)==1;  t_end = 100;           end
    if isempty(n_t)==1;    n_t = 1001;            end
    if isempty(n_w)==1;    n_w = 1001;            end
    if isempty(w_min)==1;  w_min = min(hydro.w);  end
    if isempty(w_max)==1;  w_max = max(hydro.w);  end

    % Interpolate to the given t and w
    t = linspace(0,t_end,n_t);
    w = linspace(w_min,w_max,n_w);
    % N = length(t)*sum(hydro.dof)*sum(hydro.dof);
    N = sum(hydro.dof) * sum(hydro.dof);

    % Calculate the impulse response function for radiation
    n = 0;

    hydro.ra_K = nan (sum(hydro.dof), sum(hydro.dof), length(t));

    for i = 1:sum(hydro.dof)

        for j = 1:sum(hydro.dof)

            ra_B = interp1 (hydro.w, squeeze (hydro.B(i,j,:)), w);

            % note that the following depends on 'broadcasting', w.*t(:)
            % results in a (n_w x n_t) matrix, each row corresponding to
            % one time value, we then integrate along the columns to get
            % the ra_K for every time point
            hydro.ra_K(i,j,:) = (2/pi) * trapz (w, ra_B .* (cos (w.*t(:)) .* w), 2);
            
            % hydro.ra_L(i,j,:) = (2/pi) * trapz (w, ra_B .* (sin (w.*t(:))), 2);  %Not used

    %         for k = 1:length(t)
    %             hydro.ra_K(i,j,k) = (2/pi) * trapz (w, ra_B .* (cos (w*t(k)) .* w));
    %             % hydro.ra_L(i,j,k) = (2/pi)*trapz(w,ra_B.*(sin(w*t(k))));  %Not used
                n = n + 1;
    %         end

        end
        waitbar(n/N)
    end

    hydro.ra_t = t;
    hydro.ra_w = w;

end