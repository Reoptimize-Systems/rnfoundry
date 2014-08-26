function [branches, coilsperbranch] = branchfac2coilconfig_AM(activecoilsperphase, branchfac)
% finds a suitible coil configuration given the desired number of machine
% active coils per phase and a desired fraction of the coils per phase into
% which coil blocks are to be split
%
% Syntax
%
% [branches, coilsperbranch] = branchfac2coilconfig_AM(activecoilsperphase, branchfac)
%
% Input
%
%   activecoilsperphase - number active coils per phase, in most cases this
%     will simply be the total number of coils in one phase.
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

    % now determine the number of parallel coil branches to use
    branchcomp = branchfac * activecoilsperphase;
    
    factors = factor2(activecoilsperphase)';
    
    NearestFacStruct = ipdm(branchcomp, factors, ...
                            'Subset', 'NearestNeighbor', ...
                            'Result', 'Structure');
                        
    branches = factors(NearestFacStruct.columnindex, NearestFacStruct.rowindex);
    
    coilsperbranch = activecoilsperphase / branches;
    
end