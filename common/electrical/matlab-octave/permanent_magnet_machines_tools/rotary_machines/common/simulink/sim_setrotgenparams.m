function sim_setrotgenparams(design, simoptions, modelname, blockname, makepersistent)
% Adds the design and simoptions structures to the UserData field of a
% rotary generator simulink block.
%
% Syntax
%
%   sim_setrotgenparams(design, simoptions, modelname, blockname)
%
% Description 
%
% Adds the machine and simulation information necessary for a simulink
% simulation of a rotary generator to a rotary generator simulink block.
%
% Input
% 
%   design and simoptions are structures containing information on the
%     machine design to be added. 
%
%   modelname - The path name of the simulink model which contains the
%     linear generator block to be modified.
%
%   blockname - The block name of the generator block
%
%   makepersistent - (optional) treu/false  flag determining if the
%     generator data should be made persistent. If true, the data is saved
%     in the model file of the model containing the block. Defaults to
%     false if not supplied.
%
% See also: sim_addrotgen

% Created Richard Crozier 2013

    if nargin < 5
        makepersistent = false;
    end
    
    % first break the library link if there is any, as this would prevent
    % modifying the block by adding the new UserData
    set_param([modelname, '/', blockname], 'LinkStatus', 'inactive')
    
    % Setup the data structure for the model inputs
    data.design     = design;
    data.simoptions  = simoptions;
    
    % Check if initial currents have been supplied
    if ~isfield(simoptions, 'InitCurrents')
        data.InitCurrents = zeros(design.Phases, 1);
    else
        data.InitCurrents = simoptions.SimulinkInitCurrents;
    end
    
    % construct the name of the S-Function sub-block
    model = [modelname, '/', blockname, '/Three Phase Rotary PM Machine S-Function'];
    
    % Set the block parameters through the UserData field
    set_param(model, 'UserData', data);
    
    if makepersistent
        set_param(model, 'UserDataPersistent', 'on');
    else
        set_param(model, 'UserDataPersistent', 'off');
    end
    
    % add the windings info
    model = [modelname, '/', blockname, '/Windings'];
    
    % Set the winding block parameters
    set_param(model, 'SelfImpedance1', sprintf('[%g %g]', design.PhaseResistance(1), design.PhaseInductance(1)));
    set_param(model, 'SelfImpedance2', sprintf('[%g %g]', design.PhaseResistance(2), design.PhaseInductance(1)));
    set_param(model, 'SelfImpedance3', sprintf('[%g %g]', design.PhaseResistance(3), design.PhaseInductance(1)));
     
    if numel(data.design.PhaseInductance) == 2
        set_param(model, 'MutualImpedance', sprintf('[0 %g]', data.design.PhaseInductance(2)));
    else
        set_param(model, 'MutualImpedance', '[0 0]');
    end

end