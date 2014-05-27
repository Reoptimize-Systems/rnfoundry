function ptables = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, RlVRp)
% generates tables of performance data a multiple speed and load points for
% a slotted radial flux machine design
%
% Syntax
%
% [ptables] = performancetables_RADIAL_SLOTTED(design, simoptions, rpm, RlVRp)
%
% Input
%
%   design, simoptions - the machine design and simoptions structures
%
%   rpm - vector of rpm values
%
%   RlVRp - vector of load resistance to phase resistance ratios
%
% Output
%
%   ptables - structure containing tables of output data at each
%     corresponding combination of speed and loading. The fields of ptables
%     will contain at least the following:
%
%       PowerLoadMean
%       Efficiency
%       IPhasePeak
%       IPhaseRms
%       ICoilPeak
%       ICoilRms
%       EMFPhasePeak
%       EMFPhaseRms
%       JCoilPeak
%       JCoilRms
%       EnergyLoadTotal
%       PowerPhaseRMean
%       TorquePtoMean
%       TorquePtoPeak
%       PowerLossEddyMean
%       PowerLossMean
%       FrequencyPeak
%
%     each table will be a (n x m) matrix where n is the number of speed
%     points and m is the number of load ratio points supplied.
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
    ptables = performancetables_RADIAL(design, simoptions, rpm, RlVRp, outfields);
    
end