% Test_lineargen_three_phase_pm_machine_s_function_mdl
%
%


clear data

data.design = test_design_TM_SLOTLESS ();

data.design.RlVRp = 10;

data.simoptions = struct();
data.simoptions.GetVariableGapForce = false;

data.design = completedesign_TM_SLOTLESS (data.design, data.simoptions);

[data.design, data.simoptions] = simfun_TM_SLOTTED (data.design, data.simoptions);

[data.design, data.simoptions] = finfun_TM_SLOTTED (data.design, data.simoptions);

%% Now set up the simulink simulation

open_system('Test_lineargen_three_phase_pm_machine_s_function');

% use the ode15s solver which is more appropriate for electrical machines
set_param('Test_lineargen_three_phase_pm_machine_s_function','Solver','ode15s','StopTime','5')
