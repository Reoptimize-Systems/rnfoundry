function [FemmProblem, design] = pmsmfemmprob(design, pos, varargin)
% generates a two-pole model of the linear pmsm in an mfemm FemmProblem
% structure
%
% PMSM machine diagram with dimensions shown below:
%
%           hba        ht           g      hm  hbf
%         <-----><------------><---------><--><---->
%         |      |_____________            ___|    |  ^
%         |                    | ^     ^  |   |    |  :
%         |                    | : Wt  :  |   |    |  :
%         |       _____________| v     :  |   |    |  :
%         |      |                     :  |   |    |  :
%         |    ^ |_____________        :  |   |    |  :
%         |    :               |       :  |   |    |  :
%         | Wc :               |    Wm :  |   |    |  : Wp
%         |	   :  _____________|       :  |   |    |  :
%         |    v |               ^Ws   :  |   |    |  :
%         |      |_____________  v     :  |   |    |  :
%         |                    |       :  |   |    |  :
%         |                    |       :  |   |    |  :
%         |	      _____________|       v  |___|    |  :
%         |      |                            |    |  v   
%
%

    Inputs.NWindingLayers = 1;
    Inputs.CoilCurrents = [0, 0, 0];
    
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    % Add some labels
    simwidth = sqrt((design.ht+design.hba+design.hm+design.g+design.hbf)^2 + (2*design.Wp)^2);
    mshsz = min(simwidth / 75, sqrt(simwidth^2 + (2*design.Wp)^2) / 150); % changed mesh density
    
    design = setfieldifabsent(design, 'Wsg', 0);
    design = setfieldifabsent(design, 'hsb', 0);
    design = setfieldifabsent(design, 'hsg', 0);

    [FemmProblem, boundnames, magnames, ...
        wirename, maginds, wireind, steelind] = femmprob_FM(design);
    
    % draw wrapped mag region
    wrapperthickness = [ 0, design.hbf; 
                         0, 2*(design.hbf+design.hm);
                         0, 3*(design.hbf+design.hm); ];
                         
                         
    [FemmProblem, wrapperthickness, leftcentres, rightcentres, wmagouternodeids] = ...
        wrappedrectmagaperiodic(FemmProblem, ...
                                design.Wp, ...
                                design.Wm, ...
                                design.hm, ...
                                design.hm/2 + design.hba + design.ht + design.g, ...
                                pos, ...
                                wrapperthickness, ...
                                'MagnetMaterial', maginds, ...
                                'MagnetGroup', 3, ...
                                'MeshSize', mshsz);
                            
	FemmProblem = addblocklabel_mfemm(FemmProblem, ...
                                      rightcentres(1,1), ...
                                      rightcentres(1,2), ...
                                      'BlockType', '1117 Steel', ...
                                      'InGroup', 3, ...
                                      'MaxArea', mshsz);
    
    FemmProblem = addblocklabel_mfemm(FemmProblem, ...
                                      rightcentres(2,1), ...
                                      rightcentres(2,2), ...
                                      'BlockType', 'Air', ...
                                      'MaxArea', 15 * mshsz);
                                  
    FemmProblem = addblocklabel_mfemm(FemmProblem, ...
                                      rightcentres(3,1), ...
                                      rightcentres(3,2), ...
                                      'BlockType', 'Air', ...
                                      'MaxArea', 100 * mshsz);
                                  
    % draw rhs axial flux stator with appropriate values replicate the pmsm
    % values
    [FemmProblem, outernodes, coillabellocs] = ...
        axialfluxstatorhalf2dfemmprob(6, 2, ...
                                      design.Wp, ...
                                      design.Wc, ...
                                      design.Wsg, ...
                                      design.hba * 2, ...
                                      design.ht, ...
                                      design.hsb, ...
                                      design.hsg, ...
                                      0, ...
                                      'r', ...
                                      'ToothMaterial', steelind, ...
                                      'NWindingLayers', Inputs.NWindingLayers, ...
                                      'FemmProblem', FemmProblem, ...
                                      'ToothRegionMeshSize', mshsz);
                                  
	% 
    circuits = {'Coil A', 'Coil C', 'Coil B'};
    
    coilblockprops = struct( 'BlockType', wirename, ...
                             'InCircuit', [ circuits, circuits ], ...
                             'Turns', [repmat({design.CoilTurns}, 1,3), repmat({-design.CoilTurns}, 1,3)], ...
                             'InGroup', {11, 31, 21, 12, 32, 22} ...
                            );

    for i = 1:3

        FemmProblem = addcircuit_mfemm(FemmProblem, circuits{i}, 'TotalAmps_re', Inputs.CoilCurrents(i));

        if Inputs.NWindingLayers == 1
            
            x = coillabellocs(i,1);
            y = coillabellocs(i,2);
            FemmProblem = addblocklabel_mfemm(FemmProblem, x, y, coilblockprops(i));
            
            x = coillabellocs(i+3,1);
            y = coillabellocs(i+3,2);
            FemmProblem = addblocklabel_mfemm(FemmProblem, x, y, coilblockprops(i+3));
            
        elseif Inputs.NWindingLayers == 2
            
            x = coillabellocs(2*(i-1) + 2,1);
            y = coillabellocs(2*(i-1) + 2,2);
            FemmProblem = addblocklabel_mfemm(FemmProblem, x, y, coilblockprops(i));
            
            x = coillabellocs(2*(i-1) + 1,1);
            y = coillabellocs(2*(i-1) + 1,2);
            FemmProblem = addblocklabel_mfemm(FemmProblem, x, y, 'BlockType', 'Air');
            
            x = coillabellocs(2*(i-1) + 7,1);
            y = coillabellocs(2*(i-1) + 7,2);
            FemmProblem = addblocklabel_mfemm(FemmProblem, x, y, coilblockprops(i+3));
            
            x = coillabellocs(2*(i-1) + 8,1);
            y = coillabellocs(2*(i-1) + 8,2);
            FemmProblem = addblocklabel_mfemm(FemmProblem, x, y, 'BlockType', 'Air');
            
        end

    end
    
	% complete the armature
	[FemmProblem, nodeinds,  nodeids] = addnodes_mfemm( FemmProblem, ...
                                                        [0, 0, design.hba + design.ht, design.hba + design.ht], ...
                                                        [0, 2*design.Wp, 0, 2*design.Wp] );
    
	% 
	[FemmProblem, boundind, boundnames{end+1}] = addboundaryprop_mfemm(FemmProblem, 'anti-periodic', 5);
    
    segprops.BoundaryMarker = boundnames{1};
    segprops(2).BoundaryMarker = boundnames{end};
    segprops(3).BoundaryMarker = boundnames{end};
    segprops(4).BoundaryMarker = '';
    segprops(5).BoundaryMarker = '';
    
    [FemmProblem ] = addsegments_mfemm(FemmProblem, ...
                                       [nodeids(1), nodeids(1), nodeids(2), nodeids(3), nodeids(4)], ...
                                       [nodeids(2), nodeids(3), nodeids(4), outernodes(2), outernodes(3)], ...
                                       segprops);
    
	% 
    [FemmProblem, boundind, boundnames{end+1}] = addboundaryprop_mfemm(FemmProblem, 'anti-periodic', 5);
    
    [FemmProblem ] = addsegments_mfemm(FemmProblem, ...
                                       [nodeids(3), nodeids(4)], ...
                                       [wmagouternodeids(1), wmagouternodeids(4)], ...
                                       'BoundaryMarker', boundnames{end});
    
    % Add steel label to armature
    [ FemmProblem ] = addblocklabel_mfemm(FemmProblem, design.hba/2, design.Wp, ...
                                          'BlockType', '1117 Steel', ...
                                          'maxarea', mshsz, ...
                                          'InGroup', 1);
                                      
	% Add air gap label
    [ FemmProblem ] = addblocklabel_mfemm(FemmProblem, ...
                                          design.hba + design.ht + design.g/2, ...
                                          design.Wp, ...
                                          'BlockType', 'Air');
    
    
end

function [FemmProblem, boundnames, magnames, wirename, maginds, wireind, steelind]  = femmprob_FM(design)

    FemmProblem = newproblem_mfemm('p', 'LengthUnits', 'meters', 'Depth', design.ls);
    
    % Add some boundary types
    [FemmProblem, boundind, boundnames{1}] = addboundaryprop_mfemm(FemmProblem, 'Pros A', 0, 'Phi', 90);
    
    % Add some materials
    % First Steel
    [FemmProblem, steelind] = addmaterials_mfemm(FemmProblem, '1117 Steel');
    
    % Next magnets
    for i = 1:numel(design.HcMag)
        [FemmProblem, magnames{i}, maginds(i)] = addmagnet_mfemm(FemmProblem, design.HcMag(i));
    end

    % Then Copper wire
    wirename = 'wire';
    [FemmProblem, wireind] = addmagnetwire_mfemm(FemmProblem, wirename, design.Dc);

end

function FemmProblem = femmprob_AM(design)

    

end