function [fens, gcells, BeamInfo, OuterSupportCells] = createmeshfebeamdef2_FM(x, y, z, n, znodes, zguidecon, zsupp, zframe, zextreme, guideyshift, tolerance, BeamInfo)


    if ~isfield(BeamInfo.SupportBeams, 'label')

        BeamInfo.SupportBeams.label = 1;
        BeamInfo.SupportBeams.x1x2_vector = [1 0 0];
        
        BeamInfo.OuterPoleSupports.label = 2;
        BeamInfo.OuterPoleSupports.x1x2_vector = [0 0 1];
        
        BeamInfo.OuterWebs.label = 3;
        BeamInfo.OuterWebs.x1x2_vector = [0 0 1];
        
        BeamInfo.GuideRails.label = 4;
        BeamInfo.GuideRails.x1x2_vector = [1 0 0];
        
        BeamInfo.GuideBearings.label = 5;
        BeamInfo.GuideBearings.x1x2_vector = [0 0 1];
       
    end
    
    % Preallocate a cell array to hold the beam finite element blocks
    OuterSupportCells = cell(BeamInfo.OuterPoleSupports.NoPerSide, BeamInfo.OuterPoleSupports.Sections, 3);
    
    % consolidate outer frame nodes, to a tolerance and categorise as pole
    % support only, guide bearing only, or outer web only, or a combination
    % of these
    znodemat = unifyznodes_FM(zguidecon, zsupp, zframe, tolerance);
    
    webnodes = -ones(numel(zsupp), 5);
    webcount = 1;
    
    guidebearingnodes = -ones(numel(zguidecon), 5);
    guidebearingcount = 1;
    
    supportnodes = zeros(size(znodemat,1), 5);
    
    
    startid = 1;
    endid = startid + BeamInfo.OuterPoleSupports.Sections;
    
    idsperside = numel(zframe) * (BeamInfo.OuterPoleSupports.Sections + 1);
    
    maxouterid = (2 * idsperside) + (1:4);
    
    Parameters = BeamInfo.OuterPoleSupports;
    
    BeamInfo.distn1 = [];
    BeamInfo.distn2 = [];
    
    sbeamcount = 1;
    
    for i = 1:size(znodemat, 1)
        
        if znodemat(i,4)
            
            % calculate the locations of the nodes making up the outer pole
            % support beam starting at the point
            xyz = zeros(BeamInfo.OuterPoleSupports.Sections+1, 3);

            xyz(:,1) = linspace(-x, x, BeamInfo.OuterPoleSupports.Sections+1);
            xyz(:,2) = y;
            xyz(:,3) = znodemat(i,1);

            % get the beam connection info to create the gcells with
            zframecon = [(startid:endid-1)', (startid+1:endid)'];

            Parameters.conn = zframecon;

            newnodes = (startid:endid)';
            
            if i == 1
                % create the fe nodes for the beam
                fens = fenodeset(struct('id', newnodes, 'xyz', xyz));
                newgcells = gcellset_beam3(Parameters);
                gcells = newgcells;
            else
                % create the fe nodes for the beam
                fens = cat(fens, fenodeset(struct('id', newnodes, 'xyz', xyz)));
                newgcells = gcellset_beam3(Parameters);
                gcells = cat(gcells, newgcells);
            end
            
            j = 1;
            for bs = 1:BeamInfo.OuterPoleSupports.Sections

                selcell = subset(gcells, ((sbeamcount-1) * 2 * BeamInfo.OuterPoleSupports.Sections) + bs);

                febsel = feblock_defor_ss_beam3 (struct ('gcells', selcell));

                OuterSupportCells{sbeamcount,j,3} = febsel;
                
                j = j + 1;

            end
            
            BeamInfo.distn2 = [BeamInfo.distn2; newnodes];

            % Now repeat for the other side of the frame
            xyz(:,2) = -y;

            newnodes = (startid:endid)' + idsperside;
            
            % create the fe nodes for the beam
            fens = cat(fens, fenodeset(struct('id', newnodes, 'xyz', xyz)));

            BeamInfo.distn1 = [BeamInfo.distn1; newnodes];
            
            % get the beam connection info to create the gcells with
            zframecon = [((startid:endid-1) + idsperside)', ((startid+1:endid) + idsperside)'];

            Parameters.conn = zframecon;

            newgcells = gcellset_beam3(Parameters);
            gcells = cat(gcells, newgcells);
            
            j = 1;
            for bs = (BeamInfo.OuterPoleSupports.Sections+1):(2*BeamInfo.OuterPoleSupports.Sections)

                selcell = subset(gcells, ((sbeamcount-1) * 2 * BeamInfo.OuterPoleSupports.Sections) + bs);

                febsel = feblock_defor_ss_beam3 (struct ('gcells', selcell));

                OuterSupportCells{sbeamcount,j,1} = febsel;
                
                j = j + 1;

            end
            
            supportnodes(i,:) = [znodemat(i,1), startid, endid, startid + idsperside, endid + idsperside];
            
            if znodemat(i,3)
                % node is also a guide bearing node
                guidebearingnodes(guidebearingcount, :) = [znodemat(i,1), startid, endid, startid + idsperside, endid + idsperside]; 
                guidebearingcount = guidebearingcount + 1;
            end
            
            if znodemat(i,2)
                % node is also an outer web node
                webnodes(webcount, :) = [znodemat(i,1), startid, endid, startid + idsperside, endid + idsperside];

                webcount = webcount + 1;
            end

            % Get appropriate start and end ids for the next outer pole support
            % beam
            startid = endid + 1;

            endid = startid + BeamInfo.OuterPoleSupports.Sections;
            
            sbeamcount = sbeamcount + 1;

        elseif znodemat(i,2) && znodemat(i,3)

            xyz = zeros(4,3);
            xyz(:,1) = [-x; x; -x; x];
            xyz(:,2) = [y; y; -y; -y];
            xyz(:,3) = znodemat(i,1);
            
            % create the fe nodes for the beam
            fens = cat(fens, fenodeset(struct('id', maxouterid', 'xyz', xyz)));

            % node is a guide bearing node
            guidebearingnodes(guidebearingcount,:) = [znodemat(i,1), maxouterid];
            guidebearingcount = guidebearingcount + 1;

            % node is also an outer web node
            webnodes(webcount, :) = [znodemat(i,1), maxouterid];
            webcount = webcount + 1;
            
            maxouterid = maxouterid + 4;
            supportnodes(i,:) = [znodemat(i,1), maxouterid];
            
            
        elseif znodemat(i,3)
            
            xyz = zeros(4,3);
            xyz(:,1) = [-x; x; -x; x];
            xyz(:,2) = [y; y; -y; -y];
            xyz(:,3) = znodemat(i,1);
            
            % create the fe nodes for the beam
            fens = cat(fens, fenodeset(struct('id', maxouterid', 'xyz', xyz)));
            
            % node is a guide bearing node only
            guidebearingnodes(guidebearingcount,:) = [znodemat(i,1), maxouterid];
            guidebearingcount = guidebearingcount + 1;
            
            supportnodes(i,:) = [znodemat(i,1), maxouterid];
            maxouterid = maxouterid + 4;

        elseif znodemat(i,2)
            
            xyz = zeros(4,3);
            xyz(:,1) = [-x; x; -x; x];
            xyz(:,2) = [y; y; -y; -y];
            xyz(:,3) = znodemat(i,1);
            
            % create the fe nodes for the beam
            fens = cat(fens, fenodeset(struct('id', maxouterid', 'xyz', xyz)));

            % node is an outer web node only
            webnodes(webcount, :) = [znodemat(i,1), maxouterid];
            webcount = webcount + 1;
            
            supportnodes(i,:) = [znodemat(i,1), maxouterid];
            maxouterid = maxouterid + 4;

        end
        
    end
    
    % Get the next useable node ids
    maxouterid = maxouterid(1):(maxouterid(1)+BeamInfo.OuterWebs.Sections-2);

    Parameters = BeamInfo.OuterWebs;
    
    for i = 1:size(webnodes, 1)
        
        xyz = zeros(BeamInfo.OuterWebs.Sections-1, 3);
        
        ywebs = linspace(y, -y, BeamInfo.OuterWebs.Sections+1);
        
        xyz(:,1) = -x;
        xyz(:,2) = ywebs(2:end-1);
        xyz(:,3) = webnodes(i,1);

        fens = cat(fens, fenodeset(struct('id', maxouterid', 'xyz', xyz)));

        Parameters.conn = [[webnodes(i, 2), maxouterid(1)]; ...
                           [maxouterid(1:end-1)', maxouterid(2:end)']; ...
                           [maxouterid(end), webnodes(i, 4)]];
                       
        gcells = cat(gcells, gcellset_beam3(Parameters));               
        
        % create new ids for the nex set of beams
%         maxouterid = maxouterid + BeamInfo.OuterWebs.Sections - 1;
        maxouterid = maxouterid(end)+1:(maxouterid(end)+BeamInfo.OuterWebs.Sections-1);
        
        xyz(:,1) = x;

        fens = cat(fens, fenodeset(struct('id', maxouterid', 'xyz', xyz)));

        Parameters.conn = [[webnodes(i, 3), maxouterid(1)]; ...
                           [maxouterid(1:end-1)', maxouterid(2:end)']; ...
                           [maxouterid(end), webnodes(i, 5)]];
                       
        gcells = cat(gcells, gcellset_beam3(Parameters));               
        
%         maxouterid = maxouterid + BeamInfo.OuterWebs.Sections - 1;
        maxouterid = maxouterid(end)+1:(maxouterid(end)+BeamInfo.OuterWebs.Sections-1);
    
    end
    
    
    maxouterid = maxouterid(1)-1;
    
    % create guide rail nodes for connection to bearings and therefore
    % frame
    
    Parameters = BeamInfo.GuideBearings;
    
    for i = 1:size(guidebearingnodes, 1)
        
        maxouterid = maxouterid(end)+1:(maxouterid(end) + 4);
        
        xyz = [-x,  y + guideyshift, guidebearingnodes(i,1);
                x,  y + guideyshift, guidebearingnodes(i,1);
               -x, -y - guideyshift, guidebearingnodes(i,1);
                x, -y - guideyshift, guidebearingnodes(i,1)];

        fens = cat(fens, fenodeset(struct('id', maxouterid', 'xyz', xyz)));
        
        Parameters.conn = [guidebearingnodes(i,2:5)',  maxouterid'];
        
        gcells = cat(gcells, gcellset_beam3(Parameters));
        
        guidebearingnodes(i,2:end) = maxouterid;

    end

    Parameters = BeamInfo.GuideRails;
    
    BeamInfo.clampedn = [];
    
    % Get the next useable node ids
    maxouterid = maxouterid(end)+1:(maxouterid(end)+BeamInfo.GuideRails.Sections);

    % first create a set of nodes going from the extreme end of the guide
    % rails to the top guide bearing
    xyz = zeros(BeamInfo.GuideRails.Sections, 3);
    
    zguiderail = linspace(guidebearingnodes(1,1), -zextreme, BeamInfo.GuideRails.Sections+1);
    
    xyz(:,1) = -x;
    xyz(:,2) = y + guideyshift;
    xyz(:,3) = zguiderail(2:end);

    fens = cat(fens, fenodeset(struct('id', maxouterid', 'xyz', xyz)));

    Parameters.conn = [[guidebearingnodes(1, 2), maxouterid(1)]; ...
                       [maxouterid(1:end-1)', maxouterid(2:end)']];

    BeamInfo.clampedn = [BeamInfo.clampedn, maxouterid(end)];
    
    gcells = cat(gcells, gcellset_beam3(Parameters));

    % create new ids for the nex set of beams
%     maxouterid = maxouterid + BeamInfo.GuideRails.Sections;
    maxouterid = maxouterid(end)+1:(maxouterid(end)+BeamInfo.GuideRails.Sections);

    xyz(:,1) = x;

    fens = cat(fens, fenodeset(struct('id', maxouterid', 'xyz', xyz)));

    Parameters.conn = [[guidebearingnodes(1, 3), maxouterid(1)]; ...
                       [maxouterid(1:end-1)', maxouterid(2:end)']];

    BeamInfo.clampedn = [BeamInfo.clampedn, maxouterid(end)];
    
    gcells = cat(gcells, gcellset_beam3(Parameters));

%     maxouterid = maxouterid + BeamInfo.GuideRails.Sections;
    maxouterid = maxouterid(end)+1:(maxouterid(end)+BeamInfo.GuideRails.Sections);

    xyz(:,1) = -x;
    xyz(:,2) = -y - guideyshift;

    fens = cat(fens, fenodeset(struct('id', maxouterid', 'xyz', xyz)));

    Parameters.conn = [[guidebearingnodes(1, 4), maxouterid(1)]; ...
                       [maxouterid(1:end-1)', maxouterid(2:end)']];

    BeamInfo.clampedn = [BeamInfo.clampedn, maxouterid(end)];
    
    gcells = cat(gcells, gcellset_beam3(Parameters));

%     maxouterid = maxouterid + BeamInfo.GuideRails.Sections;
    maxouterid = maxouterid(end)+1:(maxouterid(end)+BeamInfo.GuideRails.Sections);

    xyz(:,1) = x;

    fens = cat(fens, fenodeset(struct('id', maxouterid', 'xyz', xyz)));

    Parameters.conn = [[guidebearingnodes(1, 5), maxouterid(1)]; ...
                       [maxouterid(1:end-1)', maxouterid(2:end)']];

    BeamInfo.clampedn = [BeamInfo.clampedn, maxouterid(end)];
    
    gcells = cat(gcells, gcellset_beam3(Parameters));

%     maxouterid =
%     maxouterid(1):(maxouterid(1)+BeamInfo.GuideRails.Sections-2);
    maxouterid = maxouterid(end)+1:(maxouterid(end)+BeamInfo.GuideRails.Sections-1);

    % if there's more than one guide bearing, we must link in between them
    if size(guidebearingnodes, 1) > 1

        % go between guide bearings creating nodes between each one
        for i = 2:size(guidebearingnodes, 1)

            xyz = zeros(BeamInfo.GuideRails.Sections-1, 3);

            zguiderail = linspace(guidebearingnodes(i-1,1), guidebearingnodes(i,1), BeamInfo.GuideRails.Sections+1);

            xyz(:,1) = -x;
            xyz(:,2) = y + guideyshift;
            xyz(:,3) = zguiderail(2:end-1);

            fens = cat(fens, fenodeset(struct('id', maxouterid', 'xyz', xyz)));

            Parameters.conn = [[guidebearingnodes(i-1, 2), maxouterid(1)]; ...
                               [maxouterid(1:end-1)', maxouterid(2:end)']; ...
                               [maxouterid(end), guidebearingnodes(i, 2)]];

            gcells = cat(gcells, gcellset_beam3(Parameters));

            % create new ids for the nex set of beams
%             maxouterid = maxouterid + BeamInfo.GuideRails.Sections;
            maxouterid = maxouterid(end)+1:(maxouterid(end)+BeamInfo.GuideRails.Sections-1);


            xyz(:,1) = x;

            fens = cat(fens, fenodeset(struct('id', maxouterid', 'xyz', xyz)));

            Parameters.conn = [[guidebearingnodes(i-1, 3), maxouterid(1)]; ...
                               [maxouterid(1:end-1)', maxouterid(2:end)']; ...
                               [maxouterid(end), guidebearingnodes(i, 3)]];

            gcells = cat(gcells, gcellset_beam3(Parameters));

%             maxouterid = maxouterid + BeamInfo.GuideRails.Sections - 1;
            maxouterid = maxouterid(end)+1:(maxouterid(end)+BeamInfo.GuideRails.Sections-1);

            xyz(:,1) = -x;
            xyz(:,2) = -y - guideyshift;

            fens = cat(fens, fenodeset(struct('id', maxouterid', 'xyz', xyz)));

            Parameters.conn = [[guidebearingnodes(i-1, 4), maxouterid(1)]; ...
                [maxouterid(1:end-1)', maxouterid(2:end)']; ...
                [maxouterid(end), guidebearingnodes(i, 4)]];

            gcells = cat(gcells, gcellset_beam3(Parameters));

%             maxouterid = maxouterid + BeamInfo.GuideRails.Sections - 1;
            maxouterid = maxouterid(end)+1:(maxouterid(end)+BeamInfo.GuideRails.Sections-1);

            xyz(:,1) = x;

            fens = cat(fens, fenodeset(struct('id', maxouterid', 'xyz', xyz)));

            Parameters.conn = [[guidebearingnodes(i-1, 5), maxouterid(1)]; ...
                [maxouterid(1:end-1)', maxouterid(2:end)']; ...
                [maxouterid(end), guidebearingnodes(i, 5)]];

            gcells = cat(gcells, gcellset_beam3(Parameters));

%             maxouterid = maxouterid + BeamInfo.GuideRails.Sections;
            maxouterid = maxouterid(end)+1:(maxouterid(end)+BeamInfo.GuideRails.Sections-1);

        end

    end
    
    % Next create a set of nodes going from the bottom guide rail
    % bearing to the bottom of the guide rail

    maxouterid = maxouterid(1):(maxouterid(1)+BeamInfo.GuideRails.Sections-1);
    
    xyz = zeros(BeamInfo.GuideRails.Sections, 3);
    
    zguiderail = linspace(guidebearingnodes(end,1), zextreme, BeamInfo.GuideRails.Sections+1);
    
    xyz(:,1) = -x;
    xyz(:,2) = y + guideyshift;
    xyz(:,3) = zguiderail(2:end);

    fens = cat(fens, fenodeset(struct('id', maxouterid', 'xyz', xyz)));

    Parameters.conn = [[guidebearingnodes(end, 2), maxouterid(1)]; ...
                       [maxouterid(1:end-1)', maxouterid(2:end)']];

    BeamInfo.clampedn = [BeamInfo.clampedn, maxouterid(end)];
    
    gcells = cat(gcells, gcellset_beam3(Parameters));

    % create new ids for the nex set of beams
%     maxouterid = maxouterid + BeamInfo.GuideRails.Sections - 1;
    maxouterid = maxouterid(end)+1:(maxouterid(end)+BeamInfo.GuideRails.Sections);

    xyz(:,1) = x;

    fens = cat(fens, fenodeset(struct('id', maxouterid', 'xyz', xyz)));

    Parameters.conn = [[guidebearingnodes(end, 3), maxouterid(1)]; ...
                       [maxouterid(1:end-1)', maxouterid(2:end)']];

    BeamInfo.clampedn = [BeamInfo.clampedn, maxouterid(end)];
    
    gcells = cat(gcells, gcellset_beam3(Parameters));

%     maxouterid = maxouterid + BeamInfo.GuideRails.Sections - 1;
    maxouterid = maxouterid(end)+1:(maxouterid(end)+BeamInfo.GuideRails.Sections);

    xyz(:,1) = -x;
    xyz(:,2) = -y - guideyshift;

    fens = cat(fens, fenodeset(struct('id', maxouterid', 'xyz', xyz)));

    Parameters.conn = [[guidebearingnodes(end, 4), maxouterid(1)]; ...
                       [maxouterid(1:end-1)', maxouterid(2:end)']];

    BeamInfo.clampedn = [BeamInfo.clampedn, maxouterid(end)];

    gcells = cat(gcells, gcellset_beam3(Parameters));

%     maxouterid = maxouterid + BeamInfo.GuideRails.Sections - 1;
    maxouterid = maxouterid(end)+1:(maxouterid(end)+BeamInfo.GuideRails.Sections);

    xyz(:,1) = x;

    fens = cat(fens, fenodeset(struct('id', maxouterid', 'xyz', xyz)));

    Parameters.conn = [[guidebearingnodes(end, 5), maxouterid(1)]; ...
                       [maxouterid(1:end-1)', maxouterid(2:end)']];

    BeamInfo.clampedn = [BeamInfo.clampedn, maxouterid(end)];
    
    gcells = cat(gcells, gcellset_beam3(Parameters));
    
    % Connect up the outer frame beams, guide bearings and outer web nodes
    firstfreeid = maxouterid(end) + 1;
    
    BeamInfo.SupportBeams.Sections = 10;
    Parameters = BeamInfo.SupportBeams;
    
    for i = 2:size(supportnodes, 1)
        
        nL = max(1, round((supportnodes(i,1) - supportnodes(i-1,1)) * BeamInfo.SupportBeams.Sections));
        
        [fens, gcells, newids] = add_beam(fens, gcells, supportnodes(i-1,2), supportnodes(i,2), nL, Parameters, firstfreeid);

        if ~isempty(newids)
            firstfreeid = newids(end) + 1;
        end

        [fens, gcells, newids] = add_beam(fens, gcells, supportnodes(i-1,3), supportnodes(i,3), nL, Parameters, firstfreeid);

        if ~isempty(newids)
            firstfreeid = newids(end) + 1;
        end

        [fens, gcells, newids] = add_beam(fens, gcells, supportnodes(i-1,4), supportnodes(i,4), nL, Parameters, firstfreeid);

        if ~isempty(newids)
            firstfreeid = newids(end) + 1;
        end

        [fens, gcells, newids] = add_beam(fens, gcells, supportnodes(i-1,5), supportnodes(i,5), nL, Parameters, firstfreeid);

        if ~isempty(newids)
            firstfreeid = newids(end) + 1;
        end

    end
    
    %% Draw result

%     % Finite element block containing all the elements (useful for drawing)
%     feb = feblock_defor_ss_beam3 (struct ('gcells', gcells));
% 
%     geom = field(struct ('name', 'geom', ...
%                                  'dim', 3, ...
%                                  'fens', fens));
% 
%     ur  = field(struct ('name', 'ur', ...
%                                 'dim', 6, ...
%                                 'data', zeros(get(geom,'nfens'),6)));
%     % Graphics display
%     gv = graphic_viewer;
% 
%     gv = reset (gv,[]);
% 
%     % get the displacements
%     u = slice(ur, (1:3), 'u');
% 
%     % Get the rotations
%     rot = slice(ur, (4:6), 'rot');
% 
%     % Now draw the results
%     draw(feb, gv, struct ('x', geom, 'u', 0*u, 'rot', 0*rot, 'facecolor', 'none', 'drawscale', 0.2));

    
end
    
    
    
    
    