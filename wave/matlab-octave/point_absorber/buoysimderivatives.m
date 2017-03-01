function dx = buoysimderivatives (x, vBh, vBs, buoysimoptions)
% calculates the derivatives of the radiation force components for a
% heaving buoy
%
% Syntax
%
% dx = buoysimderivatives (x, buoysimoptions)
%
% Inputs
%
%  x - Vector containing the previous values of the radiation force
%    components in heave and surge, calculated through an integration of
%    continuous states (based on prony's method).
%
%  vBh - buoy velocity in heave
%
%  vBs - buoy velocity in surge
%
%  buoysimoptions - structure containing the point absorber simulation
%    parameters and sea state.
%
% Ouputs
%
%  dx - derivatives of the radiation force components in heave and surge.
%    The length of the vector will be determined by the number of radiation
%    coefficients specified in the buoysimoptions structure in the field
%    'NRadiationCoefs'. dx will be of length 2*NRadiationCoefs, the first
%    half being the heave radiation components , the rest being the surge
%    radiation components.
%

    % copy over some values (they are not modified, and so there will be no
    % memory/speed penalty from this (in theory). The purpose of this is
    % purely to provide cleaner, more readable code below.
    Hbeta = buoysimoptions.BuoyParameters.Hbeta;
    Halpha = buoysimoptions.BuoyParameters.Halpha;
    Sbeta = buoysimoptions.BuoyParameters.Sbeta;
    Salpha = buoysimoptions.BuoyParameters.Salpha;

    onetoncoeffs = 1:buoysimoptions.NRadiationCoefs;
    heaveradcoeffinds = onetoncoeffs;
    surgeradcoeffinds = buoysimoptions.NRadiationCoefs+1:2*buoysimoptions.NRadiationCoefs;
    
    % Calculate the derivative of the radiation force components in heave
    % and surge
    dx = [ real( ...
              bsxfun (@times, Hbeta(heaveradcoeffinds,1), x(heaveradcoeffinds,:)) ...
                + bsxfun (@times, Halpha(heaveradcoeffinds,:), vBh) ...
               );
           real( ...
             bsxfun (@times, Sbeta(onetoncoeffs,1), x(surgeradcoeffinds,:)) ...
               + bsxfun (@times, Salpha(onetoncoeffs,1), vBs) ...
               ) ...
         ];

end