function ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, RgVRc)
% generates tables of performance data a multiple speed and load points for
% a slotted radial flux machine design
%
% Syntax
%
% [ptables] = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, RgVRc)
%
% Input
%
%   design, simoptions - the machine design and simoptions structures
%
%   rpm - vector of rpm values
%
%   RgVRc - vector of load resistance to phase resistance ratios
%
% Output
%
%   ptables - structure containing tables of output data at each
%     corresponding combination of speed and loading
%
%

    % run the simulations and return the results using the generic axial
    % flux evaluation function
    
    outfields = {'PowerLossIronMean'};
    
    % set up the simulation functions
    simoptions.simfun = 'simfun_RADIAL_SLOTTED';
    simoptions.finfun = 'prescribedmotfinfun_RADIAL_SLOTTED';
    simoptions.odeevfun = 'prescribedmotodetorquefcn_ROTARY';
    simoptions.torquefcn = 'torquefcn_ROTARY';
    simoptions.resfun = 'prescribedmotresfun_ROTARY';
    
    % call the common radial performance tables function
    ptables = performancetables_RADIAL(design, simoptions, rpm, RgVRc, outfields);
    
end