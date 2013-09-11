function [branches, coilsperbranch] = coilconfig_linear(coilsperphase, branchfac)
% finds a suitible coil configuration given the desired number of machine
% active coils per phase and a desired fraction of the coils per phase into
% which coil blocks are to be split
%
% Syntax
%
% [branches, coilsperbranch] = coilconfig_linear(coilsperphase, branchfac)
%
% Input
%
%   coilsperphase - number active coils per phase
%
%   branchfrac - desired fraction of coils per phase to make up a desired
%   block of coils
%
% Output
%
%   branches - closest possible number of parallel branches to desired
%   fraction of coils per phase
%
%   coilsperbranch - number of series coils in each parallel branch
%

    % now determine the number of parallel coil branches to use
    branchcomp = branchfac * coilsperphase;
    
    factors = factor2(coilsperphase)';
    
    NearestFacStruct = ipdm(branchcomp, factors, ...
                            'Subset', 'NearestNeighbor', ...
                            'Result', 'Structure');
                        
    branches = factors(NearestFacStruct.columnindex, NearestFacStruct.rowindex);
    
    coilsperbranch = coilsperphase / branches;
    
end