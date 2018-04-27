

mpgaoptions.OBJ_F = 'obj_foc_pid';
mpgaoptions.SEL_F = 'sus';
mpgaoptions.XOV_F = 'recint';
mpgaoptions.MUT_F = 'mutbga';

FieldBounds = [ 0.01, 10;   % d_Kp
                1, 1000;    % d_Ki
                0.01, 10;   % q_Kp
                1, 1000; ]; % q_Ki

objfunargs = {FieldBounds, design, simoptions};

[mpgastate, mpgaoptions] = mpga (mpgaoptions, 'ObjectiveArgs', objfunargs);

