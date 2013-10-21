function name = sim_addgen(design, simoptions, modelname)
% Adds a rotary generator block to a simulink model
%
% Syntax
%
%   name = sim_addgen(design, simoptions, modelname)
%
% Description
%
% sim_addlg adds a rotary generator model to a simulink system populated
% with the appropriate data.
% 
% Input
%
%   design and simparams are structures containing information on the
%     machine design to be added. 
%
%   modelname - (optional) The path name of the simulink model to which the
%     block is to be added. If not supplied the current system will be
%     used. The current system is the current system is the system or
%     subsystem most recently clicked on.
%
% See also: sim_setlgparams, completedesign_te, checkgendesign_te

% Checked: AWM 29-11-12

    if nargin < 3
        % if a model name is not supplied, use the current system
        modelname = gcs;
    end

    % Load the simulink models
    load_system('simulink');

    % Add a level-2 S-Function block with an appropriate name, and the
    % level 2 S-Function name ('lgsfcn_lvl2')
    block = add_block( 'simulink/User-Defined Functions/Level-2 M-file S-Function', ...
                       [modelname, '/TE Lin Gen'], ...
                       'MakeNameUnique', 'on', ...
                       'mfile', 'lgsfcn_lvl2');
    
    % Get the bock name so LG parameters can be set
    name = get_param(block, 'Name');
    
    % Add the machine and simparams data to the block
    sim_setrotgenparams(design, simparams, modelname, name)
    
end