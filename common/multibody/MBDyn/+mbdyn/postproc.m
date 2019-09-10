classdef postproc < handle
% postproc class for post-processing and display of mbdyn output
%
% Description
%
% postproc is a class for the post-processing and visualisation of
% mbdyn output files.
%
% postproc Methods:
%
%  postproc - constructor
%  loadResultsFromFiles - loads mbdyn output data
%  clear - clears all data (loaded by loadResultsFromFiles)
%  plotNodeTrajectories - create plot of the trajectories of all or a
%    subset of nodes
%  plotNodeAngularVelocities - plot the angular velocities of the nodes
%  plotNodePositions - plot the positions of the nodes
%  plotNodeTrajectories - plot the trajectories of the nodes
%  plotNodeVelocities - plot the velocities of the nodes
%  availableNetCDFOutput - get the output avaialble in the netcdf file (if
%    it exists)
%  displayNetCDFVarNames - display the names of variables in the netcdf
%    file on the command line
%  getNetCDFVariable - get the contents of a variable stored in the netcdf file
%  drawStep - draw the system at the given time step index
%  animate - animate the system
%
% 
    
    properties (GetAccess = public, SetAccess = private)
        
        preProcSystem;
        resultsLoaded;
        haveNetCDF;
        nNodes;
        nodes;
        time;
        nodeNames;
        simInfo;
        mBDynOutFileName;
        movFile;
        outFile;
        logFile;
        ineFile;
        jntFile;
        txtFile;
        ncFile;
        
    end
    
    methods
        
        function self = postproc (mbdoutfilename, info)
            % constructor for the mbdyn postproc (post-processing) class
            %
            % Syntax
            %
            % mbdyn.postproc (mbdoutfilename, info)
            %
            % Description
            %
            % mbdyn.postproc constructor, creates a new mbdyn.postproc
            % object
            %
            % Input
            %
            %  mbdoutfilename - root path of the mbdyn output files, e.g.
            %   if the files are mysim_results.out, mysim_results.mov etc.
            %   this can be 'mysim_results' ( or the full path without
            %   extension e.g. /home/jblogs/mysim_results ). If a file
            %   extension is present in mbdoutfilename, the postproc class
            %   will first check for files with a base name without the
            %   extension, e.g. if mbdoutfilename was 'mysim_results.mbd',
            %   it would still check first for mysim_results.out,
            %   mysim_results.mov etc. if these were not found, it will
            %   then check for mysim_results.mbd.out, mysim_results.mbd.mov
            %   etc.
            %
            %   MBDyn can produce netcdf format ouput files, which are
            %   preferred for post-processing using this class. Using the
            %   netcd output format will make much more simulation
            %   information avaialable as otherwise only the .mov file can
            %   be imported. To make MBDyn export a netcdf format file use
            %   the 'output results' keyword in the MBDyn input file, or,
            %   if using the mbdyn.pre.system class, use the OutputResults
            %   option when creating the mbdyn.pre.system object.
            %
            %   In general, if the prefix is
            %
            %   /home/jbloggs/mysim_results
            %
            %   MBDyn will create the files:
            %
            %   /home/jbloggs/my_mbdyn_sim.frc
            %   /home/jbloggs/my_mbdyn_sim.ine
            %   /home/jbloggs/my_mbdyn_sim.out
            %   /home/jbloggs/my_mbdyn_sim.mov
            %   /home/jbloggs/my_mbdyn_sim.jnt
            %   /home/jbloggs/my_mbdyn_sim.log
            %
            %   and/or a netcdf format file:
            %
            %   /home/jbloggs/my_mbdyn_sim.nc
            %
            %   A windows example might look like
            %   C:\Users\JBloggs\Documents\my_mbdyn_sim 
            %   producing the files:
            %
            %   C:\Users\JBloggs\Documents\my_mbdyn_sim.frc
            %   C:\Users\JBloggs\Documents\my_mbdyn_sim.ine
            %   C:\Users\JBloggs\Documents\my_mbdyn_sim.out
            %   C:\Users\JBloggs\Documents\my_mbdyn_sim.mov
            %   C:\Users\JBloggs\Documents\my_mbdyn_sim.jnt
            %   C:\Users\JBloggs\Documents\my_mbdyn_sim.log
            %
            %   and/or:
            %
            %   C:\Users\JBloggs\Documents\my_mbdyn_sim.nc
            %
            %   The netcdf format is preferred for the mbdyn.postproc class
            %   as otherwise only .mov file can be used/parsed. All data is
            %   available if the netcdf format is used.
            %
            %  info - information about the system which has been solved.
            %    This can be either a structure or a class of type
            %    mbdyn.pre.system. If a mbdyn.pre.system this must
            %    represent the system which has been simulated, it is
            %    intended that it is the system which was actually used to
            %    generate the mbdyn input file for the simulation to be
            %    post-processed. 
            %
            %    Alternatively, if it is a structure, it must have the
            %    following fields:
            %
            %    InitialTime - scalar value of the intial simulation time,
            %    this 
            %
            %    TimeStep - scalar value of the time step used in the
            %      simulation
            %
            %    DefaultOrientation - string describing the default orientation
            %      format used in '.mov' output file (set using the
            %      'default orientation' keyword in the MBDyn input file,
            %      MBDyn's default is 'euler123' if this is not set). The
            %      possible values are: 'euler123', 'euler313', 'euler321',
            %      'orientation vector' or 'orientation matrix'. This
            %      field is ignored if a netcdf format output file is
            %      avaialable as the default orientation is determined from
            %      the contents of the file.
            %
            
            self.mBDynOutFileName = mbdoutfilename;
            
            self.ncFile = self.checkFile ('.nc');
            
            if nargin > 1
                if isa (info, 'mbdyn.pre.system')
                    self.preProcSystem = info;
                    
                    self.simInfo = struct ( 'InitialTime', self.preProcSystem.problems{1}.initialTime, ...
                                            ...'FinalTime', self.preProcSystem.problems{1}.finalTime, ...
                                            'TimeStep', self.preProcSystem.problems{1}.timeStep, ...
                                            'DefaultOrientation', self.preProcSystem.controlData.DefaultOrientation);
                                        
                elseif ~isempty (self.ncFile)
                    
                    warning ('A simulation info structure was provided, but this will be ignored as the simulation info will be determined from a netcdf output file.');
                    
                elseif isstruct (info)
                    if all (isfield (info, {'InitialTime', 'TimeStep', 'DefaultOrientation'}))
                        self.simInfo = info;
                    else
                        error ('Simulation info structure provided does not have all required fields');
                    end
                else
                    error ('preprocsys must be a mbdyn.pre.system object or a structure')
                end
            else
                
                if isempty (self.ncFile)
                    % no simInfo and no netcdf file to get it from
                    self.simInfo = struct ( 'InitialTime', 0, ...
                                            'TimeStep', 1, ...
                                            'DefaultOrientation', 'euler123' );
                    warning ( 'No simulation info supplied, using initial time of 0, time step of 1, and default orientation %s.', ...
                              self.simInfo.DefaultOrientation );
                else
                    % simInfo will be determined from the netcdf file
                    self.simInfo = struct ();
                end
                
            end
            
            % load the results
            loadResultsFromFiles (self, mbdoutfilename);
            
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
            %   if the files are mysim_results.out, mysim_results.mov etc.
            %   this should be 'mysim_results' ( or the full path without
            %   extension e.g. /home/jblogs/mysim_results ). If a file
            %   extension is present mbdoutfilename, the postproc class
            %   will first check for files with a base name without the
            %   extension, e.g. if mbdoutfilename was 'mysim_results.mbd',
            %   it would still check first for mysim_results.out,
            %   mysim_results.mov etc. if these were not found, it will
            %   then check for mysim_results.mbd.out, mysim_results.mbd.mov
            %   etc.
            %
            %
            
            self.clearData ();
            
            self.mBDynOutFileName = mbdoutfilename;
            
            % try to get the MBDyn redirected output file
            self.txtFile = self.checkFile ('.txt');
                
            self.ncFile = self.checkFile ('.nc');
            
            if isempty (self.ncFile)
                % no netcdf format file was found
                
                self.movFile = self.checkFile ('.mov');

                if isempty (self.movFile)
                    error ('Neither the netcdf format file or the .mov output file cound be found');
                end
                
                movdata = dlmread(self.movFile);

                self.nNodes = 0;

                posinds = 2:4;

                switch self.simInfo.DefaultOrientation

                    case {'euler123', 'euler313', 'euler321', 'orientation vector'}

                        orientinds = 5:7;
                        velinds = 8:10;
                        omegainds = 11:13;

                    case 'orientation matrix'

                        orientinds = 5:13;
                        velinds = 14:16;
                        omegainds = 17:19;

                    otherwise

                        error ('Invalid DefaultOrientation: %s', ...
                            self.simInfo.DefaultOrientation);

                end

                % parse the data
                for ind = 1:size(movdata,1)

                    label = movdata(ind,1);

                    nodename = ['node_', int2str(movdata(ind,1))];

                    if ~isfield (self.nodes, nodename)

                        self.nodes.(nodename) = struct ( ...
                                                   'Label', label, ...
                                                   'Name', nodename, ...
                                                   'Position', movdata(ind,posinds), ...
                                                   'Orientation', movdata(ind,orientinds), ...
                                                   'Velocity', movdata(ind,velinds), ...
                                                   'AngularVelocity', movdata(ind,omegainds) ...
                                                  );
                                              
                        switch self.simInfo.DefaultOrientation

                            case {'euler123', 'euler313', 'euler321', 'orientation vector'}

                                self.nodes.(nodename).Orientation = deg2rad (self.nodes.(nodename).Orientation);

                            case 'orientation matrix'

                                self.nodes.(nodename).Orientation = reshape(self.nodes.(nodename).Orientation, 3, 3)';

                        end


                        self.nNodes = self.nNodes + 1;

                    else

                        self.nodes.(nodename).Position = [self.nodes.(nodename).Position; movdata(ind,posinds)];

                        switch self.simInfo.DefaultOrientation
                            
                            case {'euler123', 'euler313', 'euler321', 'orientation vector'}
                                
                                self.nodes.(nodename).Orientation = [self.nodes.(nodename).Orientation; deg2rad (movdata(ind,orientinds))];
                                
                            case 'orientation matrix'
                                
                                self.nodes.(nodename).Orientation = cat (3, self.nodes.(nodename).Orientation, reshape(movdata(ind,orientinds), 3, 3)');
                                
                        end

                        self.nodes.(nodename).Velocity = [self.nodes.(nodename).Velocity; movdata(ind,velinds)];

                        self.nodes.(nodename).AngularVelocity = [self.nodes.(nodename).AngularVelocity; movdata(ind,omegainds)];

                    end

                end
                 

                % determine the final time
                nsteps = size(movdata,1)/self.nNodes;
                
                self.time = (1:nsteps)*self.simInfo.TimeStep;
                
                self.simInfo.FinalTime = self.time(end);

            else
                % there was a netcdf format file to load
                
                % get a description of what's in the file
                nci = ncinfo (self.ncFile);
                
                structlabels = ncread ( nci.Filename, 'node.struct');
                
                self.nNodes = numel (structlabels);
                
                % open netcdf file in read only mode (the default). Will
                % throw an error if the file cannot be opened
                ncfileid = self.netcdf_open (self.ncFile, 'NOWRITE');
                
                % make sure the file's closed when we're done
                CC = onCleanup (@() self.netcdf_close (ncfileid));
                
                % get the simulation time steps
                varid = self.netcdf_inqVarID (ncfileid, 'time');
                
                self.time = self.netcdf_getVar (ncfileid,varid);
                
                if isempty (self.time)
                    error ('Simulation time in netcdf file is empty, aborting loading of simulation data as none is available');
                end
                
                self.simInfo.InitialTime = self.time(1);
                
                self.simInfo.FinalTime = self.time(end);
                
                self.simInfo.DefaultOrientation = [];
                
                for ind = 1:numel (structlabels)
                    
                    label = structlabels(ind);
                    
                    labelstr = int2str(label); 
                    
                    nodename = [ 'node_', labelstr ];
                    
                    self.nodes.(nodename) = struct ( ...
                           'Label', double(label), ...
                           'Name', nodename, ...
                           'Position', [], ...
                           'Orientation', [], ...
                           'Velocity', [], ...
                           'AngularVelocity', [], ...
                           'OrientationType', '' ...
                          );
                      
                    
                    nodestr = sprintf('node.struct.%s', labelstr );

                    varid = self.netcdf_inqVarID (ncfileid, sprintf('%s.X', nodestr));
                    
                    self.nodes.(nodename).Position = self.netcdf_getVar (ncfileid,varid).';
                    
                    varid = self.netcdf_inqVarID (ncfileid, sprintf('%s.XP', nodestr));
                    
                    self.nodes.(nodename).Velocity = self.netcdf_getVar (ncfileid,varid).';
                    
                    varid = self.netcdf_inqVarID (ncfileid, sprintf('%s.Omega', nodestr));
                   
                    self.nodes.(nodename).AngularVelocity = self.netcdf_getVar (ncfileid,varid).';
                    
                    % now get orientation type by examining what variables
                    % are in the file
                    if isempty (self.simInfo.DefaultOrientation)
                        
                        try

                            varid = self.netcdf_inqVarID (ncfileid, sprintf('%s.R', nodestr));

                            self.simInfo.DefaultOrientation = 'orientation matrix';

                        catch

                            try
                                
                                varid = self.netcdf_inqVarID (ncfileid, sprintf('%s.Phi', nodestr));

                                self.simInfo.DefaultOrientation = 'orientation vector';

                            catch

                                try
                                    
                                    varid = self.netcdf_inqVarID (ncfileid, sprintf('%s.E', nodestr));
                                    
                                    attrvalue = netcdf.getAtt(ncfileid,varid,'description');

                                    if ~isempty(strfind (attrvalue, '123'))
                                        
                                        self.simInfo.DefaultOrientation = 'euler123';
                                        
                                    elseif ~isempty(strfind (attrvalue, '313'))
                                        
                                        self.simInfo.DefaultOrientation = 'euler313';
                                        
                                    elseif ~isempty(strfind (attrvalue, '321'))
                                        
                                        self.simInfo.DefaultOrientation = 'euler321';
                                        
                                    else
                                        
                                        error ('Could not determine euler angle type (123, 313 or 321)');
                                        
                                    end

                                catch

                                    error ('Could not get structural node orientation');

                                end

                            end

                        end

                    end
                    
                    % extract the orientation    
                    switch self.simInfo.DefaultOrientation

                        case {'euler123', 'euler313', 'euler321'}

                            varid = self.netcdf_inqVarID (ncfileid, sprintf('%s.E', nodestr));
                            
                            self.nodes.(nodename).Orientation = self.netcdf_getVar (ncfileid,varid).' ;

                        case 'orientation vector'
                            
                            varid = self.netcdf_inqVarID (ncfileid, sprintf('%s.Phi', nodestr));
                            
                            self.nodes.(nodename).Orientation = self.netcdf_getVar (ncfileid,varid).' ;

                        case 'orientation matrix'
                            
                            varid = self.netcdf_inqVarID (ncfileid, sprintf('%s.R', nodestr));
                            
                            self.nodes.(nodename).Orientation = self.netcdf_getVar (ncfileid,varid) ;

                        otherwise

                            error ('Invalid DefaultOrientation: %s', ...
                                self.simInfo.DefaultOrientation);

                    end

                end
                   
                self.haveNetCDF = true;
                
            end
            
            self.nodeNames = fieldnames (self.nodes);
            
            self.resultsLoaded = true;
            
        end
        
        function clearData (self)
            % clears all the simulation data in the object
            
            self.nodes = [];
            self.nNodes = [];
            self.nodeNames = {};
            self.mBDynOutFileName = '';
            self.time = [];
            self.txtFile = '';
            self.movFile = '';
            self.ncFile = '';
            self.haveNetCDF = false;
            
            self.resultsLoaded = false;
            
        end
        
        function clearAll (self)
            % clears all the data in the object, resetting it
            
            self.clearData ();
            
            self.preProcSystem = [];
            
        end
        
        function [ hax, hplot, hfig ] = plotNodeTrajectories (self, varargin)
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
            options.View = [];
            options.PlotAxes = [];
            
            options = parse_pv_pairs (options, varargin);
            
            % don't repeat nodes
            options.OnlyNodes = unique (options.OnlyNodes(:));
           
            if ~self.resultsLoaded
                error ('No results have been loaded yet for plotting')
            end
            
            if isempty (options.PlotAxes)
                
                hfig = figure ();
                hax = axes ();

            else
                assert (self.isAxesHandle (options.PlotAxes), ...
                    'PlotAxes must be an axes handle' );
                
                hax = options.PlotAxes;
                hfig = get (hax, 'Parent');
                
            end
            
            hold on;
            legstrings = cell (1, numel(options.OnlyNodes));
            
            for ind = 1:numel(options.OnlyNodes)
                
                hplot = plot3 ( hax, ...
                                self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Position(:,1), ...
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
            
            if isempty (options.View)
                view(hax, 3);
            else
                view (hax, options.View);
            end
            
            drawnow;
            
        end
        
        function [ hax, hplot, hfig ] = plotNodePositions (self, varargin)
            % plot the positions of the nodes
            %
            % Syntax
            %
            % hax = plotNodePositions (mbp, 'Parameter', value, ...)
            %
            % Input
            %
            %  mbp - mbdyn.postproc object
            %
            % Additional optional arguments may be supplied as
            % parameter-value pairs. The available options are:
            %
            %  'Legend' - flag determining whether to add a legend to the
            %    position plot. Default is true.
            %
            %  'Title' - flag determining whether to add a title to the
            %    position plot. Default is true.
            %
            %  'OnlyNodes' - vector of indices of nodes to plot, only these
            %    nodes will have their positions plotted. By default all
            %    node trajectories are plotted.
            %
            % Output
            %
            %  hax - handle to plot axes created
            %
            %
            
            options.Legend = true;
            options.OnlyNodes = 1:self.nNodes;
            options.Title = true;
            options.PlotAxes = [];
            
            options = parse_pv_pairs (options, varargin);
            
            [ hax, hplot, hfig ] = plotNodeQuantity (self, ...
                            'Position', 'pos_', ...
                            'Legend', options.Legend, ...
                            'OnlyNodes', options.OnlyNodes, ...
                            'PlotAxes', options.PlotAxes );
            
            if options.Title
                title (sprintf ('Node position from results file:\n%s', strrep (self.mBDynOutFileName, '_', '\_')));
            end
            
        end
        
        function [ hax, hplot, hfig ] = plotNodeVelocities (self, varargin)
            % plot the velocities of the nodes
            %
            % Syntax
            %
            % hax = plotNodePositions (mbp, 'Parameter', value, ...)
            %
            % Input
            %
            %  mbp - mbdyn.postproc object
            %
            % Additional optional arguments may be supplied as
            % parameter-value pairs. The available options are:
            %
            %  'Legend' - flag determining whether to add a legend to the
            %    velocity plot. Default is true.
            %
            %  'Title' - flag determining whether to add a title to the
            %    velocity plot. Default is true.
            %
            %  'OnlyNodes' - vector of indices of nodes to plot, only these
            %    nodes will have their velocity plotted. By default all
            %    node velocities are plotted.
            %
            % Output
            %
            %  hax - handle to plot axes created
            %
            %
            
            options.Legend = true;
            options.OnlyNodes = 1:self.nNodes;
            options.Title = true;
            options.PlotAxes = [];
            
            options = parse_pv_pairs (options, varargin);
            
            [ hax, hplot, hfig ] = plotNodeQuantity (self, ...
                            'Velocity', 'v_', ...
                            'Legend', options.Legend, ...
                            'OnlyNodes', options.OnlyNodes, ...
                            'PlotAxes', options.PlotAxes );
            
            if options.Title
                title (sprintf ('Node velocities from results file:\n%s', strrep (self.mBDynOutFileName, '_', '\_')));
            end
            
        end
        
        function [ hax, hplot, hfig ] = plotNodeAngularVelocities (self, varargin)
            % plot the angular velocities of the nodes
            %
            % Syntax
            %
            % hax = plotNodePositions (mbp, 'Parameter', value, ...)
            %
            % Input
            %
            %  mbp - mbdyn.postproc object
            %
            % Additional optional arguments may be supplied as
            % parameter-value pairs. The available options are:
            %
            %  'Legend' - flag determining whether to add a legend to the
            %    angular velocity plot. Default is true.
            %
            %  'Title' - flag determining whether to add a title to the
            %    angular velocity plot. Default is true.
            %
            %  'OnlyNodes' - vector of indices of nodes to plot, only these
            %    nodes will have their angular velocity plotted. By default
            %    all node angular velocities are plotted.
            %
            % Output
            %
            %  hax - handle to plot axes created
            %
            %
            
            options.Legend = true;
            options.OnlyNodes = 1:self.nNodes;
            options.Title = true;
            options.PlotAxes = [];
            
            options = parse_pv_pairs (options, varargin);
            
            [ hax, hplot, hfig ] = plotNodeQuantity (self, ...
                            'AngularVelocity', 'omega_', ...
                            'Legend', options.Legend, ...
                            'OnlyNodes', options.OnlyNodes, ...
                            'PlotAxes', options.PlotAxes);
            
            if options.Title
                title (sprintf ('Node angular velocities from results file:\n%s', strrep (self.mBDynOutFileName, '_', '\_')));
            end
            
        end
        
        function [ names, info ] = availableNetCDFOutput (self)
            % gets information on the available output in the mbdyn netcdf file
            %
            % Syntax
            
            if ~self.haveNetCDF
                error ('No netcdf file is available')
            end
            
            info = ncinfo (self.ncFile);
            
            % get all the variable names into a cell array
            names = { info.Variables(:).Name };
            
        end
        
        function displayNetCDFVarNames (self, mbsys)
            % gets information on the available output in the mbdyn netcdf file
            %
            % Syntax
            
            [ names, info ] = availableNetCDFOutput (self);

            for ind = 1:numel (names)
                
                out = regexp (names{ind}, '[a-zA-Z_]+\.[a-zA-Z_]+\.(\d+)$', 'tokens');
                
                if ~isempty (out)
                    description = info.Variables(ind).Attributes(1).Value;
                    
%                     labelnum = str2num (out{1}{1});
                    
                    fprintf (1, '%s : type - %s\n', names{ind}, description );
                    
                    continue;
                end
                
                description = 'no description';
                for attind = 1:numel (info.Variables(ind).Attributes)
                    if strcmpi (info.Variables(ind).Attributes(attind).Name, 'description')
                        description = info.Variables(ind).Attributes(attind).Value;
                        fprintf (1, '%s : %s\n', names{ind}, description );
                        break;
                    end
                end
                
                
                
            end
            
        end
        
        function hax = plotNetCDFVar (self, component, quantity, varargin)
            % plot the output of a components (such as a joint or other element)
            %
            % Syntax
            %
            % hax = plotOutput (mbp, component_type, quantity, label) 
            % hax = plotOutput (..., 'Parameter', value);
            %
            % Description
            %
            % 
            %
            % Input
            %
            %  mbp - mbdyn.postproc object
            %
            %  Addtional optional arguments are supplied as paramter-value
            %  pairs. The available options are:
            %
            %  'Legend' - flag determining whether to add a legend to the
            %    position plot. Default is true.
            %
            %  'Title' - flag determining whether to add a title to the
            %    position plot. Default is true.
            %
            % Output
            %
            %  hfig - handle to figure created
            %
            %  hax - handle to plot axes created
            %
            %
            
            options.Legend = true;
            options.Title = true;
            
            options = parse_pv_pairs (options, varargin);
            
            if ~self.haveNetCDF
                error ('No netcdf file is available for plotting')
            end
            
            var = self.getNetCDFVariable (component, quantity);
            
            dims = size (var);
            
            plotdim = find (dims == length (self.time));
            
            if isempty (plotdim)
                error ('Variable not the right size to plot against time.');
            end
            
            hfig = figure;
            hax = axes;
            
            ColOrd = get(hax,'ColorOrder');
            [m,n] = size(ColOrd);

            legstrings = {};
            
            ndims = numel (dims);
            switch ndims
                
                case 2
                    
                    if plotdim == 1
                        plot (self.time, var);
                    else
                        plot (self.time, var.');
                    end
                    
%                 case 3
                    
                    
                    
                otherwise
                    
                    error ('Can only plot vars with 2 dimensions');
                    
            end
            
            xlabel ('Time [s]');
            
            y_label_str = '';
            if ischar (component)
                y_label_str = [y_label_str, component, ' : '];
            end
            y_label_str = [y_label_str, quantity];
            ylabel (y_label_str);
            set (hax, 'XLim', [self.time(1), self.time(end)]);
            
            if options.Legend
                legend (hax, legstrings, 'Interpreter', 'none', 'Location', 'BestOutside');
            end
            
        end
        
        
        function var = getNetCDFVariable (self, component, quantity)
            % gets the contents of a variable from the mbdyn netcdf output file
            %
            % Syntax
            %
            % var = getNetCDFVariable (mbp, component)
            % var = getNetCDFVariable (mbp, component, quantity)
            %
            % Description
            %
            % getNetCDFVariable loads data from the netcdf file created by
            % MBDyn. If not netcdf file is available then an error will be
            % thrown. The displayNetCDFVarNames method can be useful to
            % discover the names used for quantities in the netcdf file for
            % use with getNetCDFVariable.
            %
            % Input
            %
            %  mbp - mbdyn.postproc object
            %
            %  component - either a character vector, an mbdyn.pre.element
            %    object, or an mbdyn.pre.node object. If it is a character
            %    vector it must be the full name of a variable in the
            %    netcdf file, including the specifi quantity to be output
            %    e.g. 'elem.joint.12.F'. In this case contnets of
            %    "quantity" will be ignored. e.g.
            %
            %    var = getNetCDFVariable (mbp, 'elem.joint.12.F');
            %
            %    If component is an mbdyn.pre.element object, or an
            %    mbdyn.pre.node object then these must be objects which
            %    were used to create the system using mbdyn.pre.system, or
            %    at least have the correct labels in their "label"
            %    properties. "quantity" must then contain a string
            %    indicating what quantity should be loaded, e.g:
            %
            %    var = getNetCDFVariable (mbp, joint_obj, 'F')
            %
            %  quantity - character vecto indicating what variable is to be
            %    returned when "component" contains an mbdyn.pre.element
            %    object, or an mbdyn.pre.node object.
            %
            % Output
            %
            %  var - the contents of the variable loaded from the output
            %    netcdf file.
            %
            %
            % See also: mbdyn.postproc.displayNetCDFVarNames
            %
            
            if ~self.haveNetCDF
                error ('No netcdf file is available for plotting')
            end
            
            if nargin < 3
                quantity = [];
            end
            
            if ischar (component) && isempty (quantity)
                
                varname = component;
                
            elseif ischar (component) && ischar (quantity)
                
                varname = [ component, '.', quantity];
                
            elseif isa (component, 'mbdyn.pre.element')
                
                varname = [ 'elem.' component.netCDFName, '.' int2str(component.label), '.', quantity ];
                
            elseif isa (component, 'mbdyn.pre.node')
                
                varname = [ 'node.' component.netCDFName, '.' int2str(component.label), '.', quantity ];
                
            end
            
            ncid = self.netcdf_open (self.ncFile, 'NOWRITE');
            
            CC = onCleanup (@() self.netcdf_close (ncid));
            
            varid = self.netcdf_inqVarID (ncid, varname);
                            
            var = self.netcdf_getVar (ncid, varid);
            
        end
        
        
        function mov = animate (self, varargin)
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
            % Addtional optional arguments are supplied as parameter-value
            % pairs. The available options are:
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
            %  'ExternalDrawFcn' - function handle or string to function
            %    which  will be called after drawing the scene is complete.
            %    Intended to be used by external programs to add their own
            %    features to the scene. The funciton must take two
            %    arguments, the first is the handle to the axes contining
            %    the plot, the second is the time step index of the
            %    simulation.
            %
            %  'View' - viewpoint of plot in the same format as any of the
            %    single input styles of the 'view' function. See output of
            %    'help view' for more information. The contents of the
            %    'View' option is passed directly to the view function to
            %    set the axes view.
            %
            %  'FigPositionAndSize' - figure position and size as a four
            %    element vector in the same format as a figure's 'Position'
            %    property, i.e. [ x, y, w, h ] where x and y are the
            %    coordinates of the lower left corner of the figure and w
            %    and h are the height and width respectively.
            %
            %  'VideoFile' - optional file name into which the aimation
            %    will be saved as a video. Default is empty, in which case
            %    no video file will be produced unless a VideoWriter object
            %    is supplied instead through the 'VideoWriter' option (see
            %    below).
            %
            %  'VideoWriter' - optional VideoWriter object which will be
            %    used to create a video instead of creating one internally.
            %    This option cannot be used at the same time as the
            %    'VideoFile' option. If this option is used, the
            %    VideoProfile option (see below) is ignored. However,
            %    animate will still modify the FrameRate property orf the
            %    supplied VideoWriter object internally, and the VideoSpeed
            %    option (see below) will also be applied. Default is empty,
            %    in which case no video file will be produced unless a
            %    video file path is supplied instead through the
            %    'VideoFile' option (see above).
            %
            %  'VideoSpeed' - optional speed multiplier for the video when 
            %    using the 'VideoFile'or 'VideoWriter' option. It must be a
            %    scalar value greater than 0. The animation will play a
            %    speed multiplied by this factor (by changing the video
            %    frame rate) rather than real time. Default is 1 (no
            %    speedup).
            %
            %  'VideoProfile' - optional character vector indicating what
            %    type of video compression to use when using the
            %    'VideoFile' option. This option coresponds to the profile
            %    option for the VideoWrite class, soo see the help for
            %    VideoWriter to see what the possible options are. Default
            %    is 'Motion JPEG AVI'.
            %
            % Output
            %
            %  mov - if a video file is produced, this is VideoWriter
            %   object used to create it. If no video file was
            %   created/requested, it is an empty matrix.
            %
            % See also:  mbdyn.postproc.drawStep
            %
            %
            
            
            if ~self.resultsLoaded
                error ('No results have been loaded yet for plotting')
            end
            
            % the following options will be passed directly to drawStep,
            % and checked there
            options.PlotAxes = [];
            options.DrawLabels = false;
            options.AxLims = [];
            options.PlotTrajectories = false;
            options.NodeSize = [];
            options.DrawMode = 'wireghost';
            options.DrawNodes = true;
            options.DrawBodies = true;
            options.Skip = 1;
            options.Light = false;
            options.OnlyNodes = 1:self.nNodes;
            options.ExternalDrawFcn = [];
            options.View = [];
            options.FigPositionAndSize = [];
            options.Title = false;
            
            % the following options are specific to animate and must be
            % checked here
            options.VideoFile = [];
            options.VideoSpeed = 1;
            options.VideoProfile = 'Motion JPEG AVI';
            options.VideoWriter = [];
            options.VideoQuality = 75;
            
            options = parse_pv_pairs (options, varargin);
            
            mbdyn.pre.base.checkNumericScalar (options.VideoSpeed, true, 'VideoSpeed');
            if ~isempty (options.VideoFile)
                assert (ischar (options.VideoFile), ...
                    'VideoFile must be a character vector containing a valid file name for the video.');
            end
            if ~isempty (options.VideoProfile)
                assert (ischar (options.VideoProfile), ...
                    'VideoProfile must be a character vector containing a valid file video profile.');
            end
            if ~isempty (options.VideoWriter)
                assert (isa (options.VideoWriter, 'VideoWriter'), ...
                    'VideoWriter must be a VideoWriter object' );
            end
            if ~isempty (options.VideoFile) && ~isempty (options.VideoWriter)
                error ('The VideoFile and VideoWriter options cannot both be used at the same time');
            end
            mbdyn.pre.base.checkNumericScalar (options.VideoQuality, true, 'VideoQuality');
            assert (options.VideoQuality > 0 && options.VideoQuality <= 100, ...
                'VideoQuality msut be > 0 and <= 100' );
            
            if self.resultsLoaded == false
                error ('No results have been loaded yet');
            end
            
            
            make_video = ~isempty (options.VideoFile) || ~isempty (options.VideoWriter);
            
            plotdata.HAx = [];
            plotdata.htraj = [];
            
            for tind = 1:options.Skip:size(self.nodes.(self.nodeNames{1}).Position,1)
                
                if tind > 1
                    % clear old node label drawings
                    if options.DrawLabels
                        delete (plotdata.HNodeLabels);
                    end
                    
%                     if ~isempty (plotdata.htraj)
%                         delete (plotdata.htraj);
%                     end
                    
                    % clear the axes
                    cla (plotdata.HAx);
                    
                    % clear the FigPositionAndSize option so we don't keep
                    % setting the size every step
                    options.FigPositionAndSize = [];
                    
                    % get the current view, instead of continuing to force
                    % the initial view so we can rotate the thing during
                    % playback
%                     options.View = get (plotdata.HAx, 'View');
                    
                end
                
                plotdata = self.drawStep ( tind, ...
                              'NodeSize', options.NodeSize, ...
                              'DrawLabels', options.DrawLabels, ...
                              'AxLims', options.AxLims, ...
                              'DrawMode', options.DrawMode, ...
                              'DrawNodes', options.DrawNodes, ...
                              'DrawBodies', options.DrawBodies, ...
                              'Light', options.Light, ...
                              'PlotAxes', plotdata.HAx, ...
                              'Title', options.Title, ...
                              'OnlyNodes', options.OnlyNodes, ...
                              'ForceRedraw', true, ...
                              'ExternalDrawFcn', options.ExternalDrawFcn, ...
                              'View', options.View, ...
                              'FigPositionAndSize', options.FigPositionAndSize );
                          
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

                        axXlim(1) = axXlim(1) - stretchfactor*abs(axXlim(1)) - abs (min (xExcursion(1), 0));
                        axXlim(2) = axXlim(2) + stretchfactor*abs(axXlim(2)) + abs (max (xExcursion(2), 0));

                        axYlim(1) = axYlim(1) - stretchfactor*abs(axYlim(1)) - abs (min (yExcursion(1), 0));
                        axYlim(2) = axYlim(2) + stretchfactor*abs(axYlim(2)) + abs (max (yExcursion(2), 0));

                        axZlim(1) = axZlim(1) - stretchfactor*abs(axZlim(1)) - abs (min (zExcursion(1), 0));
                        axZlim(2) = axZlim(2) + stretchfactor*abs(axZlim(2)) + abs (max (zExcursion(2), 0));
                        
                        options.AxLims = [axXlim; axYlim; axZlim];
                    end
                    
                end
                
                if options.PlotTrajectories
                    
                    for indii = 1:self.nNodes
                    
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
%                         holdstatus = ishold (plotdata.HAx);
%                         hold (plotdata.HAx, 'on');
%                         plotdata.htraj(indii) = line ( plotdata.HAx, ...
%                                 self.nodes.(self.nodeNames{indii}).Position(1:options.Skip:tind,1), ...
%                                 self.nodes.(self.nodeNames{indii}).Position(1:options.Skip:tind,2), ...
%                                 self.nodes.(self.nodeNames{indii}).Position(1:options.Skip:tind,3) );
%                         if holdstatus
%                             hold (plotdata.HAx, 'on');
%                         else
%                             hold (plotdata.HAx, 'off');
%                         end
                        
                    end
                    
                end
                
                set (plotdata.HAx, 'XLim', options.AxLims(1,1:2));
                set (plotdata.HAx, 'YLim', options.AxLims(2,1:2));
                set (plotdata.HAx, 'ZLim', options.AxLims(3,1:2));
                
%                 drawnow;
                
                if make_video

                    if tind == 1
                        
                        if isempty (options.VideoWriter)
                            
                            mov = VideoWriter(options.VideoFile, options.VideoProfile); %#ok<TNMLP>
                            
                            mov.Quality = options.VideoQuality;
                            
                        else
                            mov = options.VideoWriter;
                        end
                        
                        FPS = options.VideoSpeed ...
                                   * numel (1:options.Skip:size(self.nodes.(self.nodeNames{1}).Position,1)) ...
                                            / (self.simInfo.FinalTime - self.simInfo.InitialTime);
                                        
                        mov.FrameRate = FPS;
                        
                        F = getframe (plotdata.HFig);
                        open (mov);
                        writeVideo (mov,F);
                    else
                        F = getframe (plotdata.HFig);
                        writeVideo (mov,F);
                    end
                    
                else
                    mov = [];
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
            %  'ExternalDrawFcn' - function handle or string to function
            %    which  will be called after drawing the scene is complete.
            %    Intended to be used by external programs to add their own
            %    features to the scene. The funciton must take two
            %    arguments, the first is the handle to the axes contining
            %    the plot, the second it the time step index of the
            %    simulation.
            %
            %  'View' - viewpoint of plot in the same format as any of the
            %    single input styles of the 'view' function. See output of
            %    'help view' for more information. The contents of the
            %    'View' option is passed directly to the view function to
            %    set the axes view.
            %
            %  'FigPositionAndSize' - figure position and size as a four
            %    element vector in the same format as a figure's 'Position'
            %    property, i.e. [ x, y, w, h ] where x and y are the
            %    coordinates of the lower left corner of the figure and w
            %    and h are the height and width respectively.
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
            options.ForceRedraw = false;
            options.ExternalDrawFcn = [];
            options.View = [];
            options.FigPositionAndSize = [];
            
            options = parse_pv_pairs (options, varargin);
            
            assert (tind <= numel(self.time), ...
                'Requested time index is greater than the maximum time index available.');
            
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
            
            if ~isempty (options.FigPositionAndSize)
                set (plotdata.HFig, 'Position', options.FigPositionAndSize);
            end
            
            if ~isempty (self.preProcSystem)
                
                defaultOrientationIsMatrix = strcmpi (self.simInfo.DefaultOrientation, 'orientation matrix');
                
                for indii = 1:self.nNodes
                    % set the node positions and orientations
                    self.preProcSystem.setNodePosition (self.nodes.(self.nodeNames{indii}).Label, self.nodes.(self.nodeNames{indii}).Position(tind,:)');
                    
                    if defaultOrientationIsMatrix
                        % get the transpose of the matrix as mbdyn writes
                        % it out row-wise, not columnwise
                        om = mbdyn.pre.orientmat (self.simInfo.DefaultOrientation, self.nodes.(self.nodeNames{indii}).Orientation(:,:,tind));
                    else
                        om = mbdyn.pre.orientmat (self.simInfo.DefaultOrientation, self.nodes.(self.nodeNames{indii}).Orientation(tind,:));
                    end
                    self.preProcSystem.setNodeOrientation (self.nodes.(self.nodeNames{indii}).Label, om);
                end
                
                % draw the system
                self.preProcSystem.draw ( 'AxesHandle', plotdata.HAx, ...
                    'Mode', options.DrawMode, ...
                    'Bodies', options.DrawBodies, ...
                    'StructuralNodes', options.OnlyNodes, ...
                    'Joints', false, ...
                    'Light', options.Light, ...
                    'ForceRedraw', options.ForceRedraw);
                
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
                                      self.time(tind), ...
                                      filestr ) ...
                          );
                else
                    title ( plotdata.HAx, ...
                            sprintf ( 'System plot at time index %d for MBDyn results file:\n%s', ...
                                      tind, filestr ) ...
                          );
                end
            end
                        
            if ~isempty (options.ExternalDrawFcn)
                feval (options.ExternalDrawFcn, plotdata.HAx, tind);
            end
            
            daspect (plotdata.HAx, [1,1,1]);
            
            if isempty (options.View)
                view (plotdata.HAx, 3);
            else
                view (plotdata.HAx, options.View);
            end

            drawnow;
            
        end

    end
    
    
    % internal methods
    methods (Access = protected)
        
        function fpath = checkFile (self, mbext)
            
            [pathstr, name, ext] = fileparts (self.mBDynOutFileName);
            
            fpath1 = fullfile (pathstr, [name, mbext]);
            fpath = fpath1;
            if ~exist (fpath1, 'file')
                
                if isempty (ext)
                    fpath = '';
                else
                    % in case file root name has dots in it
                    fpath2 = fullfile (pathstr, [name, ext, mbext]);
                    fpath = fpath2;
                    if ~exist  (fpath, 'file')
                        fpath = '';
                    end
                end
                
            end
            
        end

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
        
        function [ hax, hplot, hfig ] = plotNodeQuantity (self, fieldname, legprefix, varargin)
            % plot a variable associated with each node
            %
            % Syntax
            %
            % hax = plotNodeQuantity (mbp, fieldname, legprefix, 'Parameter', value, ...)
            %
            % Input
            %
            %  mbp - mbdyn.postproc object
            %
            %  fieldname - character vector containing the name of the
            %    variable associated with each node which is to be plotted
            %
            %  legprefix - prefix for legend entries
            %
            % Addtional optional arguments are supplied as parameter-value
            % pairs. The avaialable options are:
            %
            %  'Legend' - flag determining whether to add a legend to the
            %    position plot. Default is true.
            %
            %  'Title' - flag determining whether to add a title to the
            %    position plot. Default is true.
            %
            %  'OnlyNodes' - vector of indices of nodes to plot, only these
            %    nodes will have their positions plotted. By default all
            %    node trajectories are plotted.
            %
            %  'PlotAxes' - handle to the axes in which to perform the
            %    plot. Default is empty in which case a new figure and axes
            %    are created.
            %
            % Output
            %
            %  hfig - handle to figure created
            %
            %  hax - handle to plot axes created
            %
            %
            
            options.Legend = true;
            options.OnlyNodes = 1:self.nNodes;
            options.Title = true;
            options.PlotAxes = [];
            
            options = parse_pv_pairs (options, varargin);
            
            % don't repeat nodes
            options.OnlyNodes = unique (options.OnlyNodes(:));
            
            if ~self.resultsLoaded
                error ('No results have been loaded yet for plotting')
            end
            
            if isempty (options.PlotAxes)
                
                hfig = figure ();
                hax = axes ();

            else
                assert (self.isAxesHandle (options.PlotAxes), ...
                    'PlotAxes must be an axes handle' );
                
                hax = options.PlotAxes;
                hfig = get (hax, 'Parent');
                
            end
            
            ColOrd = get(hax,'ColorOrder');
            [m,n] = size(ColOrd);

            legstrings = {};
            
            hold on;
            for ind = 1:numel(options.OnlyNodes)
                
                ColOrd = circshift (ColOrd, 1, 1);
                
                Col = ColOrd(1,:);
                
                hplot = plot ( self.time, self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).(fieldname)(:,1), 'LineStyle', '-', 'Color', Col );
                hplot(2) = plot ( self.time, self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).(fieldname)(:,2), 'LineStyle', '--', 'Color', Col );
                hplot(3) = plot ( self.time, self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).(fieldname)(:,3), 'LineStyle', ':', 'Color', Col );
                
                node_label = self.nodes.(self.nodeNames{options.OnlyNodes(ind)}).Label;
                
                legstrings = [ legstrings, { sprintf('Node %d %sx', node_label, legprefix), ...
                                             sprintf('Node %d %sy', node_label, legprefix), ...
                                             sprintf('Node %d %sz', node_label, legprefix) ...
                                           } ...
                             ];
                    
            end
            hold off
            
            xlabel ('MBDyn Sim Time [s]');
            set (hax, 'XLim', [self.time(1), self.time(end)]);
            
            if options.Legend
                legend (hax, legstrings, 'Interpreter', 'none', 'Location', 'BestOutside');
            end
            
        end
        
%         function checkAxes (self, hax)
%             % checks if there is a valid set of axes to plot to, and if
%             % not, creates them
%             
%             % try to figure out if there is a valid set of axes to plot to
%             if isempty (hax)
%                 
%                 % plot in the existing axes if possible as no new axes
%                 % supplied
%                 if self.isAxesHandle (self.drawAxesH)
%                     % use existing axes
%                     
%                     if ~ishghandle (self.drawAxesH)
%                         % axes no longer exists, or figure has been closed
%                         self.drawAxesH = [];
%                         self.deleteAllDrawnObjects ();
%                         self.transformObject = [];
%                         self.needsRedraw = true;
%                         % make a new set of axes to draw to
%                         self.makeAxes ();
%                     end
%                     
%                 elseif isempty (self.drawAxesH)
%                     % make figure and axes
%                     self.makeAxes ();
%                     self.needsRedraw = true;
%                     
%                 else
%                     error ('drawAxesH property is not empty or an axes handle');
%                 end
%             
%             elseif self.isAxesHandle (hax)
%                 % an axes has been supplied, so we plot to this new axes
%                 
%                 if isoctave
%                     drawnow ();
%                 end
%                 
%                 if ~ishghandle (hax)
%                     error ('provided axes object is not valid');
%                 end
%                 self.drawAxesH = hax;
%                 % abandon existing shape objects and transform object
%                 % TODO: should we delet objects in old axes?
%                 self.shapeObjects = {};
%                 self.transformObject = [];
%                 % we need to redraw, as we're plotting in a different set
%                 % of axes
%                 self.needsRedraw = true;
%                 
%             else
%                 
%                 error ('hax must be an axes handle or empty');
%                 
%             end
%             
%         end
        
        function ret = isAxesHandle (self, hax)
            % test if variable is axes handle
            
            if isoctave ()
                ret = isaxes (hax);
            else
                ret = isa (hax, 'matlab.graphics.axis.Axes');
            end

        end

    end
    
    % octave and matlab compatible netcdf functions
    methods (Static)
        
        function ncid = netcdf_open (varargin)
            
            if isoctave
                if numel (varargin) > 1
                    switch varargin{2}
                        
                        case 'NOWRITE'
                            mode = netcdf_getConstant ('NC_NOWRITE');
                        case 'WRITE'
                            mode = netcdf_getConstant ('NC_WRITE');
                        case 'SHARE'
                            mode = netcdf_getConstant ('NC_SHARE');
                        otherwise
                            
                            error ('unrecognised mode');
                    end
                end
                
                ncid = netcdf_open (varargin{1}, mode);
            else
                ncid = netcdf.open (varargin{:});
            end
            
        end
        
        
        function netcdf_close (varargin)
            
            if isoctave
                netcdf_close (varargin{:});
            else
                netcdf.close (varargin{:});
            end
            
        end
        
        
        function varid = netcdf_inqVarID (varargin)
            
            if isoctave
                varid = netcdf_inqVarID (varargin{:});
            else
                varid = netcdf.inqVarID (varargin{:});
            end
            
        end
        
        
        function data = netcdf_getVar (varargin)
            
            if isoctave
                data = netcdf_getVar (varargin{:});
            else
                data = netcdf.getVar (varargin{:});
            end
            
        end
        
        
    end
    
end