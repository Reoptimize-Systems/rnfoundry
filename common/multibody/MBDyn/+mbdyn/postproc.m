classdef postproc < handle
    % postproc   class for post-processing and display of mbdyn output
    %
    % Description
    %
    % postproc is a class for the post-processing and visualisation of
    % mbdyn output files.
    %
    % postproc Methods:
    %  postproc - constructor
    %  loadResultsFromFiles - loads mbdyn output data
    %  clear - clears all data (loaded by loadResultsFromFiles)
    %  plotNodeTrajectories - create plot of the trajectories of all or a
    %     subset of nodes
    %  drawStep - draw the system at the given time step index
    %  animate - animate the system
    %
    % 
    
    properties (GetAccess = public, SetAccess = private)
        
        preProcSystem;
        resultsLoaded;
        nNodes;
        nodes;
        nodeNames;
        simInfo;
        mBDynOutFileName;
        
    end
    
    methods
        
        function self = postproc (mbdoutfilename, info)
            % constructor for the mbdyn postproc (post-processing) class
            %
            % Syntax
            %
            % postproc(mbdoutfilename, info)
            %
            % Description
            %
            % mbdyn.postproc(mbdoutfilename, info) creates a new postpoc
            % class
            %
            % Input
            %
            %  mbdoutfilename - root path of the mbdyn output files, e.g.
            %   if the files are mysim_results.out, mysim_results.move etc.
            %   this should be 'mysim_results' ( or the full path without
            %   extension e.g. /home/jblogs/mysim_results )
            %
            %  info - information about the system which has been solved.
            %    This can be either a structure or a class of type
            %    mbdyn.pre.system. If a mbdyn.pre.system this must
            %    represent the system which has been simulated, it is
            %    intended that it is the system which was actually used to
            %    generate the mbdyn input file for the simulation to be
            %    post-processed. Alternatively, if it is a structure, it
            %    must have the following fields:
            %
            %    InitialTime - scalar value of the intial simualtion time
            %    TimeStep - scalar value of the time step used in the
            %      simulation
            %
            
            
            if nargin > 1
                if isa (info, 'mbdyn.pre.system')
                    self.preProcSystem = info;
                    
                    self.simInfo = struct ( 'InitialTime', self.preProcSystem.problems{1}.initialTime, ...
                                            ...'FinalTime', self.preProcSystem.problems{1}.finalTime, ...
                                            'TimeStep', self.preProcSystem.problems{1}.timeStep );
                    
                elseif isstruct (info)
                    if all (isfield (info, {'InitialTime', 'TimeStep'}))
                        self.simInfo = info;
                    else
                        error ('Simulation info structure provided does not have all required fields');
                    end
                else
                    error ('preprocsys must be a mbdyn.pre.system object or a structure')
                end
            else
                self.simInfo = struct ('InitialTime', 0, 'TimeStep', 1);
                warning ('No simulation info supplied, using initial time of 0 and time step of 1.')
            end
            
            % load the results
            loadResultsFromFiles (self, mbdoutfilename)
            
        end
        
        function loadResultsFromFiles (self, mbdoutfilename)
            % load results from mbdyn output files
            %
            % Syntax
            %
            % loadResultsFromFiles (mbp, mbdoutfilename)
            %
            % Description
            %
            % loadResultsFromFiles loads the data output by MBDyn during a
            % simulation in preparation for post-processing. Any previous
            % data loded from output files is cleared.
            %
            % Input
            %
            %  mbp - mbdyn.postproc object
            %
            %  mbdoutfilename - root path of the mbdyn output files, e.g.
            %   if the files are mysim_results.out, mysim_results.move etc.
            %   this should be 'mysim_results' ( or the full path without
            %   extension e.g. /home/jblogs/mysim_results )
            %
            %
            
            [pathstr, name, ext] = fileparts (mbdoutfilename);
            
            self.clearData ();
            
            self.mBDynOutFileName = mbdoutfilename;
            
%             self.parseLogFile (logfile);
            
            % load .mov file
            movfile = fullfile (pathstr, [name, '.mov']);
            
            if ~exist (movfile, 'file')
                % in case file root name has dots in it
                movfile = fullfile (pathstr, [name, ext, '.mov']);
                if ~exist  (movfile, 'file')
                    error ('MBDyn output file %s not found', movfile);
                end
            else
                % make ext empty sp we can add it to the other files below
                % whether it's needed or not
                ext = '';
            end

            movdata = dlmread(movfile);
            
            self.nNodes = 0;
            
            for ind = 1:size(movdata,1)
                
                label = movdata(ind,1);
                
                nodename = ['node_', int2str(movdata(ind,1))];
                
                if ~isfield (self.nodes, nodename)
                    
                    self.nodes.(nodename) = struct ( ...
                                               'Label', label, ...
                                               'Name', nodename, ...
                                               'Position', movdata(ind,2:4), ...
                                               'Orientation', movdata(ind,5:13), ...
                                               'Velocity', movdata(ind,14:16), ...
                                               'AngularVelocity', movdata(ind,17:19) ...
                                              );
                                          
                    self.nNodes = self.nNodes + 1;
                    
                else
                    
                    self.nodes.(nodename).Position = [self.nodes.(nodename).Position; movdata(ind,2:4)];
                    
                    self.nodes.(nodename).Orientation = [self.nodes.(nodename).Orientation; movdata(ind,5:13)];
                    
                    self.nodes.(nodename).Velocity = [self.nodes.(nodename).Velocity; movdata(ind,14:16)];
                    
                    self.nodes.(nodename).AngularVelocity = [self.nodes.(nodename).AngularVelocity; movdata(ind,17:19)];
                                
                end
                
            end
            
            % determine the final time
            self.simInfo.FinalTime = self.simInfo.InitialTime + (size(movdata,1)-1)*self.simInfo.TimeStep;
            
            % load .log file
            %
            % The log file contains a description of the problem in it's
            % intial state including bodies and their labels
            logfile = fullfile (pathstr, ext, [name, '.log']);
            
            % load .ine file
            inefile = fullfile (pathstr, ext, [name, '.ine']);
            
            % load .frc file
            frcfile = fullfile (pathstr, ext, [name, '.frc']);
            
            self.nodeNames = fieldnames (self.nodes);
            
            self.resultsLoaded = true;
            
        end
        
        function clearData (self)
            % clears all the data in the object, resetting it
            
            self.nodes = [];
            self.nNodes = [];
            self.resultsLoaded = false;
            self.nodeNames = {};
            self.mBDynOutFileName = '';
            
        end
        
        function clearAll (self)
            % clears all the data in the object, resetting it
            
            self.clearData ();
            
            self.preProcSystem = [];
            simInfo = [];
        end
        
        function [hfig, hax] = plotNodeTrajectories (self, varargin)
            % plot the trajectories of the nodes
            %
            % Syntax
            %
            % [hfig, hax] = plotNodeTrajectories (mbp, 'Parameter', value, ...)
            %
            % Input
            %
            %  mbp - mbdyn.postproc object
            %
            %  Addtional optional arguments are supplied as paramter-value
            %  pairs. The avaialable options are:
            %
            %  'AxLims' - (3 x 2) matrix containing limits the plot axes
            %    should be set to, each row is the x, y and z axes limits
            %    respectively. Uselful to focus on a particular region of
            %    the system. If not supplied, suitable axes limits will be
            %    determined based on the motion data.
            %
            %  'Legend' - flag determining whether to add a legend to the
            %    trajectory plot. Default is true.
            %
            %  'Title' - flag determining whether to add a title to the
            %    trajectory plot. Default is true.
            %
            %  'OnlyNodes' - vector of indices of nodes to plot, only these
            %    nodes will have their trajectoris plotted. By default all
            %    node trajectories are plotted.
            %
            % Output
            %
            %  hfig - handle to figure created
            %
            %  hax - handle to plot axes created
            %
            %
            
            options.AxLims = [];
            options.Legend = true;
            options.OnlyNodes = 1:self.nNodes;
            options.Title = true;
            
            options = parse_pv_pairs (options, varargin);
            
            % don't repeat nodes
            options.OnlyNodes = unique (options.OnlyNodes(:));
           
            if ~self.resultsLoaded
                error ('No results have been loaded yet for plotting')
            end
            
            hfig = figure;
            hax = axes;
            
            hold on;
            legstrings = cell (1, numel(options.OnlyNodes));
            
            for ind = 1:numel(options.OnlyNodes)
                
                plot3 ( self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(:,1), ...
                        self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(:,2), ...
                        self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(:,3), ...
                        '-' );
                
                if isempty (options.AxLims)
                    if ind == 1
                        axXlim = [ min(self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(:,1)), ...
                                  max(self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(:,1)) ];
                        axYlim = [ min(self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(:,2)), ...
                                  max(self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(:,2)) ];
                        axZlim = [ min(self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(:,3)), ...
                                  max(self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(:,3)) ];

                    else

                        axXlim = [ min( axXlim(1), min(self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(:,1)) ), ...
                                  max( axXlim(2), max(self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(:,1)) ) ];
                        axYlim = [ min( axYlim(1), min(self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(:,2)) ), ...
                                  max( axYlim(2), max(self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(:,2)) ) ];
                        axZlim = [ min( axZlim(1), min(self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(:,3)) ), ...
                                  max( axZlim(2), max(self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(:,3)) ) ];

                    end
                end
                
                if options.Legend
                    legstrings{ind} = self.nodeNames{options.OnlyNodes(ind)};
                end
            
            end
            
            hold off;
            
            % set data aspect ratios so that there is no distortion
            daspect (hax, [1,1,1]);
            
            if isempty (options.AxLims)
                stretchfactor = 0.25;
                
                minsize = mean ([diff(axXlim), diff(axYlim), diff(axZlim)]) / 10;
                
                axXlim(1) = axXlim(1) - stretchfactor*abs(axXlim(1));
                axXlim(2) = axXlim(2) + stretchfactor*abs(axXlim(2));
                
                if diff (axXlim) < minsize
                    axXlim(1) = axXlim(1) - (minsize/2);
                    axXlim(2) = axXlim(2) + (minsize/2);
                end
                
                axYlim(1) = axYlim(1) - stretchfactor*abs(axYlim(1));
                axYlim(2) = axYlim(2) + stretchfactor*abs(axYlim(2));
                
                if diff (axYlim) < minsize
                    axYlim(1) = axYlim(1) - (minsize/2);
                    axYlim(2) = axYlim(2) + (minsize/2);
                end
                
                axZlim(1) = axZlim(1) - stretchfactor*abs(axZlim(1));
                axZlim(2) = axZlim(2) + stretchfactor*abs(axZlim(2));
                
                if diff (axZlim) < minsize
                    axZlim(1) = axZlim(1) - (minsize/2);
                    axZlim(2) = axZlim(2) + (minsize/2);
                end
                
            else
                axXlim = options.AxLims (1,1:2);
                axYlim = options.AxLims (2,1:2);
                axZlim = options.AxLims (3,1:2);
            end
            
            set (hax, 'XLim', axXlim);
            set (hax, 'YLim', axYlim);
            set (hax, 'ZLim', axZlim);
            
            xlabel (hax, 'x');
            ylabel (hax, 'y');
            zlabel (hax, 'z');
            
            if options.Legend
                legend (hax, legstrings, 'Interpreter', 'none');
            end
            
            if options.Title
                title (sprintf ('Node trajectories plot for results file:\n%s', strrep (self.mBDynOutFileName, '_', '\_')));
            end
            
            view(3);
            drawnow;
            
        end
        
        
        function animate (self, varargin)
            % animate the nodes and bodies of the system
            %
            % Syntax
            %
            % animate (mbp, 'Parameter', value) 
            %
            % Input
            %
            %  mbp - mbdyn.postproc object
            %
            %  Addtional optional arguments are supplied as parameter-value
            %  pairs. The available options are:
            %
            %  'PlotAxes' - handle to plot axes in which to draw the
            %    trajectories. If not supplied a new figure and axes are
            %    created for the plot.
            %
            %  'AxLims' - (3 x 2) matrix containing limits the plot axes
            %    should be set to, each row is the x, y and z axes limits
            %    respectively. Uselful to focus on a particular region of
            %    the system. If not supplied, suitable axes limits will be
            %    determined based on the motion data.
            %
            %  'Title' - flag determining whether to add a title to the
            %    system plot. Default is true.
            %
            %  'DrawMode' - string indicating the drawing mode. Can be one
            %    of: 'solid', 'ghost', 'wire', 'wireghost'.
            %
            %  'DrawLabels' - flag indicating whether to draw node labels,
            %    default is false.
            %
            %  'DrawNodes' - flag determining whether to draw nodes,
            %    default is true
            %
            %  'DrawBodies' - flag determining whether to draw bodies,
            %    default is true.
            %
            %  'OnlyNodes' - vector of indices of nodes to plot, only these
            %    nodes will be plotted. By default all nodes are plotted if
            %    'DrawNodes' is true. This setting is ignored if
            %    'DrawNodes' is false.
            %
            %  'Light' - flag determining whether to light the scene to
            %    give a more realistic look. Best used with 'solid' option.
            %
            %  'NodeSize' - scalar value indicating how large the nodes
            %    should be drawn (in the same units as the plot). If not
            %    supplied, the nodes are drawn with a size based on the
            %    axes limits.
            %
            %
            
            if ~self.resultsLoaded
                error ('No results have been loaded yet for plotting')
            end
            
            options.DrawLabels = false;
            options.AxLims = [];
            options.PlotTrajectories = isempty (self.preProcSystem);
            options.NodeSize = [];
            options.DrawMode = 'wireghost';
            options.DrawNodes = true;
            options.DrawBodies = true;
            options.Skip = 1;
            options.Light = false;
            options.VideoFile = [];
            options.VideoSpeed = 1;
            options.OnlyNodes = 1:self.nNodes;
            
            options = parse_pv_pairs (options, varargin);
            
            if self.resultsLoaded == false
                error ('No results have been loaded yet');
            end
            
            plotdata.HAx = [];
            
            for tind = 1:options.Skip:size(self.nodes.(self.nodeNames{1}).Position,1)
                
                if tind > 1
                    % clear old node label drawings
                    if options.DrawLabels
                        delete (plotdata.HNodeLabels);
                    end
                    % clear the axes
                    cla (plotdata.HAx);
                end
                
                plotdata = self.drawStep(tind, ...
                              'NodeSize', options.NodeSize, ...
                              'DrawLabels', options.DrawLabels, ...
                              'AxLims', options.AxLims, ...
                              'DrawMode', options.DrawMode, ...
                              'DrawNodes', options.DrawNodes, ...
                              'DrawBodies', options.DrawBodies, ...
                              'Light', options.Light, ...
                              'PlotAxes', plotdata.HAx, ...
                              'Title', false, ...
                              'OnlyNodes', options.OnlyNodes);
                          
                if tind == 1
                    
                    if isempty (options.AxLims)
                        % determine suitible axes limits for the entire sim
                        % duration
                        
                        for ind = 1:numel(options.OnlyNodes)

                            if ind == 1
                                xExcursion = [ min( self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(:,1) ...
                                                    - self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(1,1) ), ...
                                               max( self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(:,1) ...
                                                    - self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(1,1)) ...
                                             ];
                                          
                                yExcursion = [ min( self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(:,2) ...
                                                   - self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(1,2) ), ...
                                              max( self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(:,2) ...
                                                   - self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(1,2)) ...
                                             ];
                                               
                                zExcursion = [ min(self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(:,3) ...
                                                   - self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(1,3)), ...
                                              max( self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(:,3) ...
                                                   - self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(1,3) ) ...
                                             ];

                            else

                                xExcursion = [ min( xExcursion(1), ...
                                                    min( self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(:,1) ...
                                                         - self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(1,1) ) ), ...
                                               max( xExcursion(2), ...
                                                    max( self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(:,1) ...
                                                         - self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(1,1) ) ) ];
                                                     
                                yExcursion = [ min( yExcursion(1), ...
                                                    min(self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(:,2) ...
                                                        - self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(1,2)) ), ...
                                               max( yExcursion(2), ...
                                                    max(self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(:,2) ...
                                                        - self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(1,2)) ) ];
                                          
                                zExcursion = [ min( zExcursion(1), ...
                                                    min(self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(:,3) ...
                                                        - self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(1,3)) ), ...
                                               max( zExcursion(2), ...
                                                    max(self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(:,3) ...
                                                        - self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(1,3)) ) ];

                            end

                        end
                        
                        axXlim = get (plotdata.HAx, 'XLim');

                        axYlim = get (plotdata.HAx, 'YLim');

                        axZlim = get (plotdata.HAx, 'ZLim');
                        
                        % stretch the axes a bit and add node excursions to
                        % axes limits
                        stretchfactor = 0.1;

                        axXlim(1) = axXlim(1) - stretchfactor*abs(axXlim(1)) - min (xExcursion(1), 0);
                        axXlim(2) = axXlim(2) + stretchfactor*abs(axXlim(2)) + max (xExcursion(2), 0);

                        axYlim(1) = axYlim(1) - stretchfactor*abs(axYlim(1)) - min (yExcursion(1), 0);
                        axYlim(2) = axYlim(2) + stretchfactor*abs(axYlim(2)) + max (yExcursion(2), 0);

                        axZlim(1) = axZlim(1) - stretchfactor*abs(axZlim(1)) - min (zExcursion(1), 0);
                        axZlim(2) = axZlim(2) + stretchfactor*abs(axZlim(2)) + max (zExcursion(2), 0);
                        
                        options.AxLims = [axXlim; axYlim; axZlim];
                    end
                    
                end
                
                for indii = 1:self.nNodes
                    
                    if options.PlotTrajectories
                        if tind == 1
                            htraj(indii) = animatedline ( plotdata.HAx, ...
                                self.nodes.(self.nodeNames{indii}).Position(tind,1), ...
                                self.nodes.(self.nodeNames{indii}).Position(tind,2), ...
                                self.nodes.(self.nodeNames{indii}).Position(tind,3) ...
                                );
                            
                        else
                            
                            addpoints ( htraj(indii), ...
                                self.nodes.(self.nodeNames{indii}).Position(tind,1), ...
                                self.nodes.(self.nodeNames{indii}).Position(tind,2), ...
                                self.nodes.(self.nodeNames{indii}).Position(tind,3) );
                        end
                    end
                    
                end
                
                set (plotdata.HAx, 'XLim', options.AxLims(1,1:2));
                set (plotdata.HAx, 'YLim', options.AxLims(2,1:2));
                set (plotdata.HAx, 'ZLim', options.AxLims(3,1:2));
                
                if ~isempty(options.VideoFile)

                    if tind == 1
                        FPS = options.VideoSpeed ...
                               * numel (1:options.Skip:size(self.nodes.(self.nodeNames{1}).Position,1)) ...
                                        / (self.simInfo.FinalTime - self.simInfo.InitialTime);
                        
                        mov = VideoWriter(options.VideoFile); %#ok<TNMLP>
                        mov.FrameRate = FPS;
                        mov.Quality = 100;

                        F = getframe(plotdata.HFig);
                        open(mov);
                        writeVideo(mov,F);
                    else
                        F = getframe(plotdata.HFig);
                        writeVideo(mov,F);
                    end
                    
                end
            
            end
            
            if ~isempty (options.VideoFile)
                close(mov);
            end
            
        end
    
        function plotdata = drawStep (self, tind, varargin)
            % draw the system at the given time step index
            %
            % Syntax
            %
            % plotdata = drawStep (mbp, tind, 'Parameter', value)
            %
            % Input
            %
            %  mbp - mbdyn.postproc object
            %
            %  Addtional optional arguments are supplied as parameter-value
            %  pairs. The available options are:
            %
            %  'PlotAxes' - handle to plot axes in which to draw the
            %    trajectories. If not supplied a new figure and axes are
            %    created for the plot.
            %
            %  'AxLims' - (3 x 2) matrix containing limits the plot axes
            %    should be set to, each row is the x, y and z axes limits
            %    respectively. Uselful to focus on a particular region of
            %    the system. If not supplied, suitable axes limits will be
            %    determined based on the motion data.
            %
            %  'Title' - flag determining whether to add a title to the
            %    system plot. Default is true.
            %
            %  'DrawMode' - string indicating the drawing mode. Can be one
            %    of: 'solid', 'ghost', 'wire', 'wireghost'.
            %
            %  'DrawLabels' - flag indicating whether to draw node labels,
            %    default is false.
            %
            %  'DrawNodes' - flag determining whether to draw nodes,
            %    default is true
            %
            %  'DrawBodies' - flag determining whether to draw bodies,
            %    default is true.
            %
            %  'OnlyNodes' - vector of indices of nodes to plot, only these
            %    nodes will be plotted. By default all nodes are plotted if
            %    'DrawNodes' is true. This setting is ignored if
            %    'DrawNodes' is false.
            %
            %  'Light' - flag determining whether to light the scene to
            %    give a more realistic look. Best used with 'solid' option.
            %
            %  'NodeSize' - scalar value indicating how large the nodes
            %    should be drawn (in the same units as the plot). If not
            %    supplied, the nodes are drawn with a size based on the
            %    axes limits.
            %
            % Output
            %
            %  hfig - handle to figure created
            %
            %  hax - handle to plot axes created
            %
            %
            
            if ~self.resultsLoaded
                error ('No results have been loaded yet for plotting')
            end
            
            options.PlotAxes = [];
            options.DrawLabels = false;
            options.AxLims = [];
            options.NodeSize = [];
            options.DrawMode = 'wireghost';
            options.DrawNodes = true;
            options.DrawBodies = true;
            options.Light = false;
            options.Title = true;
            options.OnlyNodes = 1:self.nNodes;
            
            options = parse_pv_pairs (options, varargin);
            
            if options.DrawNodes == false
                options.OnlyNodes = false;
            end
            
            if isempty (options.PlotAxes)
                plotdata.HFig = figure;
                plotdata.HAx = axes;
            else
                plotdata.HAx = options.PlotAxes;
                plotdata.HFig = get (plotdata.HAx, 'Parent');
            end
            
            if ~isempty (self.preProcSystem)
                
                for indii = 1:self.nNodes
                    % set the node positions and orientations
                    self.preProcSystem.setNodePosition (self.nodes.(self.nodeNames{indii}).Label, self.nodes.(self.nodeNames{indii}).Position(tind,:)');
                    % get the transpose of the matrix as mbdyn writes it
                    % out row-wise, not columnwise
                    om = mbdyn.pre.orientmat ('orientation', reshape (self.nodes.(self.nodeNames{indii}).Orientation(tind,:), 3, 3)');
                    self.preProcSystem.setNodeOrientation (self.nodes.(self.nodeNames{indii}).Label, om);
                end
                
                % draw the system
                self.preProcSystem.draw ( 'AxesHandle', plotdata.HAx, ...
                    'Mode', options.DrawMode, ...
                    'Bodies', options.DrawBodies, ...
                    'StructuralNodes', options.OnlyNodes, ...
                    'Joints', false, ...
                    'Light', options.Light );
                
            end
            
            if isempty (options.AxLims)
                axXlim = get (plotdata.HAx, 'XLim');
                axYlim = get (plotdata.HAx, 'YLim');
                axZlim = get (plotdata.HAx, 'ZLim');
            else
                axXlim = options.AxLims (1,1:2);
                axYlim = options.AxLims (2,1:2);
                axZlim = options.AxLims (3,1:2);
            end
            
            if isempty (options.NodeSize)
                % choose suitable node size based on the axis size if not
                % supplied
                plotdata.nodeSize = repmat (0.03 * mean ( diff ([axXlim; axYlim; axZlim ], 1, 2) ), ...
                             1, self.nNodes);
            else
                if isscalar (options.NodeSize)
                    plotdata.nodeSize = repmat (options.NodeSize, 1, self.nNodes);
                elseif numel (options.NodeSize) ~= self.nNodes
                    error ('You supplied a vector of node sizes which is not the same size as the number of nodes.')
                end
            end
            
            for indii = 1:self.nNodes

                if options.DrawLabels
                    plotdata.HNodeLabels(indii) ...
                        = text ( plotdata.HAx, ...
                                 self.nodes.(self.nodeNames{indii}).Position(tind,1) + plotdata.nodeSize(indii), ...
                                 self.nodes.(self.nodeNames{indii}).Position(tind,2), ...
                                 self.nodes.(self.nodeNames{indii}).Position(tind,3) + plotdata.nodeSize(indii), ...
                                 self.nodeNames{indii}, ...
                                 'Interpreter', 'none');
                end
                
            end
            
            set (plotdata.HAx, 'XLim', axXlim);
            set (plotdata.HAx, 'YLim', axYlim);
            set (plotdata.HAx, 'ZLim', axZlim);
            
            if options.Title
                filestr = strrep (self.mBDynOutFileName, '_', '\_');
                if ~isempty (self.preProcSystem)
                    title ( plotdata.HAx, ...
                            sprintf ( 'System plot at time t = %.2fs for MBDyn results file:\n%s', ...
                                      (tind-1)*self.preProcSystem.problems{1}.timeStep, ...
                                      filestr ) ...
                          );
                else
                    title ( plotdata.HAx, ...
                            sprintf ( 'System plot at time index %d for MBDyn results file:\n%s', ...
                                      tind, filestr ) ...
                          );
                end
            end
            
            daspect (plotdata.HAx, [1,1,1]);
            
            view(3);
            
            drawnow;
        end

    end
    
    
    % internal methods
    methods (Access = protected)

        function h = drawnode (self, pos, rot, sx, sy, sz)
            
            [X,Y,Z] = cylinder([0 ones(1,30) 0.98 ones(1,10) 0]/sqrt(2), 4);
            
            [X,Y,Z] = rotsurf(X, Y, Z, [0 0 pi/4]);
            
            X = sx*X; 
            Y = sy*Y; 
            Z = sz*(Z-0.5);
            
            [X,Y,Z] = rotsurf(X, Y, Z, rot);
            
            X = X + pos(1); 
            Y = Y + pos(2); 
            Z = Z + pos(3);
            
            color = 'b';
            alph = 1;
            
            h = surf(X, Y, Z, 'FaceColor', color, 'FaceAlpha', alph, ...
                'EdgeAlpha', alph, 'EdgeColor', 'none');
            
        end
%         
%         
%         function h = moveobj (self, h, pos, rot)
%             
%             
%         end

    end
    
end