function [design, simoptions] = chrom2design_TORUS_SLOTLESS(simoptions, Chrom, varargin)
% converts a chromosomal representation of a slotless torus machine to a
% full machine design in preparation for simulation
%
% Syntax
%
% [design, simoptions] = chrom2design_TORUS_SLOTLESS(simoptions, Chrom, 'Parameter', Value)
%
%

% Copyright 2012-2013 Richard Crozier

% 
%        FieldBounds = [ 0.1,   0.9;%  1. taucoVtaucsm 
                                    %  2. tyVtm
%                    0.1,   6.0;    %  3. tcVtm
%                    0.0,   0.5;    %  4. gVtm
                                    %  5. NSupports
                                    %  6. tausuppVsuppspace
%                                   %  7  tsuppbVtbio
%                    0.1,   1.0;    %  8. tbiiVtbio
%                    0.1,   2.0;    %  9. tbioVtm
%                    1.005, 4.0;    % 10. tmVtaumm
%                    0.5,   0.95;   % 11. taummVtaupm
%                    0.5,   0.99;   % 12. RsVRbi 
%                    0.8,   0.99;   % 13. RbiVRmi
%                    0.1,   0.9;    % 14. RmiVRmo
%                    0.5,   3.0;    % 15. Rmo
%                    0.5,   15.0;   % 16. RlVRp
%                    0.2,   0.85;   % 17. kfill
%                    0,     1;      % 18. DcAreaFac
%                    5,     50;     % 19. pole pairs
%                    0,     1;      % 20. branchfac                     
%                    0.5,   1;      % 21. modulefac
%                    1,     10 ];   % 22. NStages
%                                   % 23. 

    % number of Phases
    options.Phases = 3;
    % number of coils per pole and phase
    options.qc = fr(3,3);
    options.RlVRp = 10;
    options.CoilFillFactor = 0.86;
    options.BranchFac = 0;
    options.ModuleFac = 0;
    options.NStages = 1;
    options.Maxtc = 0.2;
    options.Maxtm = 0.2;
    options.MinAirGap = 0.5/1000;
    
    options = parseoptions(options, varargin);
    
    design.CoilFillFactor = options.CoilFillFactor;
    design.Phases = max(1, round(options.Phases));
    design.qc = options.qc;
    design.RlVRp = options.RlVRp;
    design.ModuleFac = options.ModuleFac;
    
    % convert machine ratios to actual dimensions
    design.taucoVtaucsm = Chrom(1);
    design.tyVtm = Chrom(2);
    design.tcVtm = Chrom(3);
    design.gVtm = Chrom(4);
    design.NModuleSupports = round(Chrom(5));
    design.tausuppVsuppspace = Chrom(6);
    design.tsuppbVtbio = Chrom(7);
    design.tbiiVtbio = Chrom(8);
    design.tbioVtm = Chrom(9);
    design.tmVtaumm = Chrom(10);
    design.taummVtaupm = Chrom(11);
    design.RsVRbi = Chrom(12);
    design.RbiVRmi = Chrom(13);
    design.RmiVRmo = Chrom(14);
    design.Rmo = Chrom(15);
    design.NBasicWindings = round(Chrom(16));

    factors = factor2(design.NBasicWindings)';
    
    % now determine the number of modules to use
    modulecomp = design.ModuleFac * design.NBasicWindings;
    
    NearestFacStruct = ipdm(modulecomp, factors, ...
                            'Subset', 'NearestNeighbor', ...
                            'Result', 'Structure');
                        
    design.NModules = factors(NearestFacStruct.columnindex, NearestFacStruct.rowindex);
    
    design = completedesign_TORUS_SLOTLESS(design, simoptions);
    
    if design.tc > options.Maxtc
        design.tc = options.Maxtc;
    end
    
    if design.tm > options.Maxtm
        design.tm = options.Maxtm;
        design.tbio = design.tm * design.tbioVtm;
        design.tbii = design.tbio * design.tbiiVtbio;
        design.tbi = [design.tbio, design.tbii];
        design.g = design.tm * design.gVtm;
        design.tsuppb = design.tbio * design.tsuppbVtbio;
        design.ty = design.tm * design.tyVtm;
    end
    
    if design.g < options.MinAirGap
        design.g = options.MinAirGap;
    end
    
    design.DcAreaFac = Chrom(17);
    design.BranchFac = Chrom(18);
    design.ModuleFac = Chrom(19);

    design.Hc = design.tc;
    design.Wc = design.tauco;
    
    design = preprocsystemdesign_AF(design, simoptions);

end



