function [design, simoptions] = chrom2design_TORUS_CORELESS(simoptions, Chrom, varargin)
% Converts a chromosome representing a coreless torus machine design to the
% full design and simoptions structure
%
% Syntax
%
% [design, simoptions] = chrom2design_TORUS_CORELESS_nova(simoptions, Chrom)
% [design, simoptions] = chrom2design_TORUS_CORELESS_nova(..., 'Parameter', Value)
%
%
% Input
% 
%  simoptions - simulation options structure
%
%  Chrom - vector of values representing machine design parameters
%
%  A number of parameter-value pairs in the form of a string, and the
%  associated input variable can be supplied to set additional options. If
%  not supplied, these will be given default values. The possible
%  parameter-value pairs and their default values are as follows:
%
%
%    'phases' - The number of phases in the winding, defaults to  3 
%
%    'qc' - fr (fractions) object representing the number of coils per pole
%       and phase. Defaults to fr(1,4).
%
%    'fillfactor' - Wire fill-factor in the coils, defaults to 0.86
%
%    'taucmoVtaucsm' - Outer coil pitch to max possible outer coil pitch
%       ratio (percentage coils fills available space). Defaults to 0.99
%
%    'ModuleFac' - 
%
%    'RgVRc' - Machine phase resistance to load resistance ratio. Defaults
%       to 10.
% 
%    'Maxtc' - Maximum allowed coil height (dimension tc). Designs which 
%       exceed this value will be modified to have this value. Defaults to
%       0.2.
%
%    'Maxtm' - Maximum allowed magnet height (dimension tm). Designs which 
%       exceed this value will be modified to have this value. Defaults to
%       0.04.
%
%    'MinAirGap' - The minimum possible air gap (dimension g). Designs with
%       a smaller air gap will be modified to have this value. Defaults 10
%       0.5e-3 (half a mm).
%
% Output
%
% Fully constructed design and simoptions structures for a coreless torus
% machine suitible for completion by the associated simfun and finfun
% functions.
%

% 
%        FieldBounds = [ 0.1,   0.9;    %  1. taucmiVtaucmo 
%                    0.1,   6.0;    %  2. tcVtm
%                    0.0,   0.5;    %  3. gVtm
                                    %   . NSupports
                                    %   . tausuppVsuppspace
%                                      5. tsuppbVtbio
%                    0.1,   1.0;    %  4. tbiiVtbio
%                    0.1,   2.0;    %  5. tbioVtm
%                    1.005, 4.0;    %  6. tmVtaumm
%                    0.5,   0.95;   %  7. taummVtaupm
%                    0.5,   0.99;   %  8. RsVRbi 
%                    0.8,   0.99;   %  9. RbiVRmi
%                    0.1,   0.9;    % 10. RmiVRmo
%                    0.5,   3.0;    % 11. Rmo
%                    0.5,   15.0;   % 12. RgVRc
%                    0.2,   0.85;   % 13. kfill
%                    0,     1;      % 14. DcAreaFac
%                    5,     50;     % 15. pole pairs


    % number of phases
    options.phases = 3;
    % number of coils per pole and phase
    options.qc = fr(1,4);
    % grid resistance to phase resistance ratio
    options.RgVRc = 10;
    options.fillfactor = 0.86;
    options.ModuleFac = 0;
    options.NStages = 1;
    options.taucmoVtaucsm = 0.99;
    options.Maxtc = 0.2;
    options.Maxtm = 0.04;
    options.MinAirGap = 0.5/1000;
    
    options = parseoptions(options, varargin);
    
    simoptions = setfieldifabsent(simoptions, 'WindingType', 'nonoverlapping');
    
    if strcmpi(simoptions.WindingType, 'overlapping')
        design.CoilLayers = 2;
    elseif strcmpi(simoptions.WindingType, 'nonoverlapping')
        design.CoilLayers = 1;
    end
    
    design.WindingType = simoptions.WindingType;
    
    design.fillfactor = options.fillfactor;
    design.phases = max(1, round(options.phases));
    design.qc = options.qc;
    design.taucoVtaucsm = options.taucmoVtaucsm;
    design.RgVRc = options.RgVRc;
    design.ModuleFac = options.ModuleFac;
    design.NStages = options.NStages;
    
    % convert machine ratios to actual dimensions
    design.tauciVtauco = Chrom(1);
    design.tcVtm = Chrom(2);
    design.gVtm = Chrom(3);
    design.NModuleSupports = round(Chrom(4));
    design.tausuppVsuppspace = Chrom(5);
    design.tsuppbVtbio = Chrom(6);
    design.tbiiVtbio = Chrom(7);
    design.tbioVtm = Chrom(8);
    design.tmVtaumm = Chrom(9);
    design.taummVtaupm = Chrom(10);
    design.RsVRbi = Chrom(11);
    design.RbiVRmi = Chrom(12);
    design.RmiVRmo = Chrom(13);
    design.Rmo = Chrom(14);
    design.NBasicWindings = round(Chrom(15));

    factors = factor2(design.NBasicWindings)';
    
    % now determine the number of modules to use
    modulecomp = design.ModuleFac * design.NBasicWindings;
    
    NearestFacStruct = ipdm(modulecomp, factors, ...
                            'Subset', 'NearestNeighbor', ...
                            'Result', 'Structure');
                        
    design.NModules = factors(NearestFacStruct.columnindex, NearestFacStruct.rowindex);
    
    design = completedesign_TORUS_CORELESS(design, simoptions);
    
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
    end
    
    if design.g < options.MinAirGap
        design.g = options.MinAirGap;
    end
    
    design.DcAreaFac = Chrom(16);
    design.BranchFac = Chrom(17);

    design.Hc = design.tc;
    design.Wc = (design.tauco - design.tauci) / 2;
    
    design = preprocsystemdesign_AF(design, simoptions);

end




