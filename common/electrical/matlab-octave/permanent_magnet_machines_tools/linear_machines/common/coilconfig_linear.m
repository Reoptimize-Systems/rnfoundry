function [branches, coilsperbranch] = coilconfig_linear(activecoilsperphase, branchfac)
% finds a suitible coil configuration given the desired number of machine
% active coils per phase and a desired fraction of the coils per phase into
% which coil blocks are to be split
%
% Syntax
%
% [branches, coilsperbranch] = coilconfig_linear(activecoilsperphase, branchfac)
%
% Input
%
%   activecoilsperphase - number active coils per phase
%
%   branchfrac - This should be a number between 0 and 1. Higher
%     numbers will result in more parallel branches, and less series coils
%     per branch. Lower numbers yield less branches and more series coils.
%
% Output
%
%   branches - closest possible number of parallel branches to desired
%     fraction of coils per phase
%
%   coilsperbranch - number of series coils in each parallel branch
%

    % call common function to do the conversion
    [branches, coilsperbranch] = branchfac2coilconfig_AM(activecoilsperphase, branchfac);
    
end