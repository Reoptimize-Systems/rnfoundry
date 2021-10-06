function [FT, FF] = forcefcn_linear_mvgfield_system(design, simoptions, xT, vT, xF, vF, xBh, xBs, vBh, vBs)

    % Currently prescribed motion forces are identical to system simulation forces
    [FT, FF] = forcefcn_linear_mvgfield_pscbmot(design, simoptions, xT, vT, xF, vF, xBh, xBs, vBh, vBs);
 
end