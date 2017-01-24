% Test_rotgen_three_phase_PM_machine_mdl.m
%
%

clear data

data.design = test_design_RADIAL_SLOTTED ();

data.design.RlVRp = 10;

data.simoptions = struct();
data.simoptions.GetVariableGapForce = false;

data.design = completedesign_RADIAL_SLOTTED (data.design, data.simoptions);

[data.design, data.simoptions] = simfun_RADIAL_SLOTTED(data.design, data.simoptions);

[data.design, data.simoptions] = finfun_RADIAL_SLOTTED(data.design, data.simoptions);

%% Now set up the simulink simulation

open_system('Test_rotgen_three_phase_pm_machine');

% set_param('Test_rotgen_three_phase_pm_machine/Three Phase Rotary PM Machine', 'UserData', data);

% use the ode15s solver which is more appropriate for electrical machines
set_param('Test_rotgen_three_phase_pm_machine','Solver','ode15s','StopTime','5')
