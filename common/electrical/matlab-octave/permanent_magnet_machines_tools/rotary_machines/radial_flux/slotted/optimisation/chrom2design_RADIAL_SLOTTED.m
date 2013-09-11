function [design, simoptions] = chrom2design_RADIAL_SLOTTED(simoptions, Chrom, varargin)
% converts a chromosomal representation of a slotted radial flux machine to
% a full machine design in preparation for simulation
%
% Syntax
%
% [design, simoptions] = chrom2design_RADIAL_SLOTTED(simoptions, Chrom, 'Parameter', Value)
%
%

% Copyright 2012-2013 Richard Crozier

    % number of phases
    options.phases = 3;
    % number of coils per pole and phase
    options.qc = fr(3,3);
    % grid resistance to phase resistance ratio
    options.RgVRc = 10;
    % coil fill factor
    options.fillfactor = 0.65;
    % branch factor, determines number of parallel and series coils
    options.BranchFac = 0;
%     options.ModuleFac = 0;
    % stator type, determines if we have an outer facing or inner facing
    % stator
    options.StatorType = 'so';
    % number of coil layers
    options.CoilLayers = 2;
    % maximum permitted radial coil height
    options.Max_tc = 0.2;
    % maximum permitted radial magnet height
    options.Max_tm = 0.05;
    % maximum stack length
    options.Max_ls = 4;
    % maximum ratio of ty to tm
    options.Max_tyVtm = 3;
    % minimum permitted air gap
    options.Min_g = 0.5/1000;
    options.Max_tsbVtc = 0.1;
    
    % parse the input options, replacing defaults
    options = parseoptions(options, varargin);
    
    % copy the resulting options into the appropriate places
    design.StatorType = options.StatorType;
    design.fillfactor = options.fillfactor;
    design.phases = max(1, round(options.phases));
    design.qc = options.qc;
    design.RgVRc = options.RgVRc;
    design.CoilLayers = options.CoilLayers;
%     design.ModuleFac = options.ModuleFac;
    
    if strcmp(design.StatorType, 'si')

        error('Stator Type ''si'' not yet supported.')

    elseif strcmp(design.StatorType, 'so')

        % convert machine ratios to actual dimensions
        design.Rbo = Chrom(1);
        design.RmoVRbo = Chrom(2);
        design.RmiVRmo = Chrom(3);
        design.RsoVRmi = Chrom(4);
        design.RtsbVRso = Chrom(5);
        design.tsgVtsb = Chrom(6);
        design.RyoVRtsb = Chrom(7);
        design.RyiVRyo = Chrom(8);
        design.thetamVthetap = Chrom(9);
        design.thetacVthetas = Chrom(10);
        design.thetasgVthetac = Chrom(11);
        design.lsVtm = Chrom(12);
        design.NBasicWindings = round(Chrom(13));
        design.DcAreaFac = Chrom(14);
        design.BranchFac = Chrom(15);
        if numel(Chrom) > 15
            design.MagnetSkew = Chrom(16);
        end

%         factors = factor2(design.NBasicWindings)';
% 
%         % now determine the number of modules to use
%         modulecomp = design.ModuleFac * design.NBasicWindings;
% 
%         NearestFacStruct = ipdm(modulecomp, factors, ...
%                                 'Subset', 'NearestNeighbor', ...
%                                 'Result', 'Structure');
% 
%         design.NModules = factors(NearestFacStruct.columnindex, NearestFacStruct.rowindex);

        design = completedesign_RADIAL_SLOTTED(design, simoptions);
        
        if (design.tsb > 0) && (design.tsb / design.tc) > options.Max_tsbVtc
            % shift the shoe base outward
            rshift = (design.tsb - (design.tc*options.Max_tsbVtc));
            design.Rtsb = design.Rtsb + rshift;
            design.tsb = design.tc*options.Max_tsbVtc;
            % recalculate the shoe gap size
            design.tsg = design.tsb * design.tsgVtsb;
            design.Rstg = design.Rso - design.tsg;
            design = updatedimsso(design);
        end

        if design.tc > options.Max_tc
            % move the stator yoke outwards to reduce the size of the slot
            rshift = (design.tc - options.Max_tc);
            design.Ryi = design.Ryi + rshift;
            design.Ryo = design.Ryo + rshift;
            design = updatedimsso(design);
        end
        
        if (design.ty / design.tm) > options.Max_tyVtm
            % move the stator yoke internal radius outwards to reduce the
            % thickness of the yoke
            rshift = design.ty - (design.tm * options.Max_tyVtm);
            design.Ryi = design.Ryi + rshift;
            design = updatedimsso(design);
        end

        if design.tm > options.Max_tm
            rshift = design.tm - options.Max_tm;
            design.tm = options.Max_tm;
            design.Ryi = design.Ryi + rshift;
            design.Ryo = design.Ryo + rshift;
            design.Rtsb = design.Rtsb + rshift;
            design.Rstg = design.Rstg + rshift;
            design.Rso = design.Rso + rshift;
            design.Rmi = design.Rmi + rshift;
            design = updatedimsso(design);
        end
        
        if design.tsb > 0 && (design.tsg < design.tsb)

            x = ((design.thetac - design.thetasg)/2) * (design.Rso - design.tsb);
            y = design.tsb - design.tsg;

            tsbangle = rad2deg(atan( y / x ));
            
            if design.tsg < 1e-5
                x = ((design.thetac - design.thetasg)/2) * (design.Rso - design.tsb);
                y = design.tsb;
                tsgangle = rad2deg(atan( y / x ));
            else
                tsgangle = inf;
            end

            if tsbangle < 15 || tsgangle < 15

                rshift = design.tsb;
                design.tsb = 0;
                design.tsg = 0;
                design.Ryi = design.Ryi + rshift;
                design.Ryo = design.Ryo + rshift;
                design.Rci = design.Ryo;
                design.Rco = design.Rco + rshift;
                design.Rso = design.Rco;
                design.Rtsb = design.Rco;
                design.Rstg = design.Rco;

                design = updatedimsso(design);
            end

        end

        if design.g < options.Min_g
            % increase the outer diameter
            rshift = (options.Min_g - design.g);
            design.Rmi = design.Rmi + rshift;
            design.Rmo = design.Rmo + rshift;
            design.Rbi = design.Rbi + rshift;
            design.Rbo = design.Rbo + rshift;
            design = updatedimsso(design);
        end
        
    end
    
    % prevent overly deep stack lengths
    if design.ls > options.Max_ls
       design.ls = options.Max_ls;
       design.lsVtm = design.ls / design.tm;
    end
    
    design.Hc = design.tc;
    design.Wc = design.thetac * design.Rcm;
    
    design = preprocsystemdesign_RADIAL(design, simoptions);

    % recalcualte the design ratios in case the have been modified
    design = recalculateratios(design);
    
end


function design = updatedimsso(design)

    % some additional radial variables
    design.Rco = design.Rtsb;
    design.Rci = design.Ryo;
    design.Rbi = design.Rmo;

    % lengths in radial direction
    design.ty = design.Ryo - design.Ryi;
    design.tc = design.Rco - design.Rci;
    design.tsb = design.Rso - design.Rtsb;
    design.g = design.Rmi - design.Rso;
    design.tm = design.Rmo - design.Rmi;
    design.tbi = design.Rbo - design.Rbi;

    % the shoe tip length
    design.Rstg = design.Rso - design.tsg;

    % mean radial position of magnets
    design.Rmm = mean([design.Rmo, design.Rmi]);
    design.Rcm = mean([design.Rci, design.Rco]);
    design.Rbm = mean([design.Rbo, design.Rbi]);
    design.Rym = mean([design.Ryi, design.Ryo]);

end


function design = recalculateratios(design)
% recalculates the design ratios of the slotted torus machine from the
% dimensions

    design.RmoVRbo = design.Rmo / design.Rbo;
    design.RmiVRmo = design.Rmi / design.Rmo;
    design.RsoVRmi = design.Rso / design.Rmi;
    design.RtsbVRso = design.Rtsb / design.Rso;
    design.RstgVRso = design.Rstg / design.Rso;
    design.RyoVRtsb = design.Ryo / design.Rtsb;
    design.RyiVRyo = design.Ryi / design.Ryo;
    design.thetamVthetap = design.thetam / design.thetap;
    design.thetacVthetas = design.thetac / design.thetas;
    design.thetasgVthetac = design.thetasg / design.thetac;
    design.lsVtm = design.ls / design.tm;

end



