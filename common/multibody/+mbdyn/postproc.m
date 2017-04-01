classdef postproc < handle
   
    properties 
        
        preProcSystem;
        resultsLoaded;
        nNodes;
        nodes;
        nodeNames;
        bodies;
        joints;
        simInfo;
        
    end
    
    methods
        
        function self = postproc (info)
            
            if nargin > 0
                if isa (info, 'mbdyn.pre.system')
                    self.preProcSystem = info;
                elseif isstruct (info) 
                    self.simInfo = info;
                else
                    error ('preprocsys must be a mbdyn.pre.system object or a structure')
                end
            else
                self.preProcSystem = [];
            end
            
        end
        
        function loadResultsFromFiles (self, mbdoutfilename)
            
%             if ~exist (mbdfilename, 'file')
%                 error ('input file not found');
%             end
            
            [pathstr, name, ~] = fileparts (mbdoutfilename);
            
            self.clear ();
            
            % load .log file
            %
            % The log file contains a description of the problem in it's
            % intial state including bodies and their labels
            logfile = fullfile (pathstr, [name, '.log']);
            
%             self.parseLogFile (logfile);
            
            % load .mov file
            movfile = fullfile (pathstr, [name, '.mov']);
            
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
            
            % load .ine file
            inefile = fullfile (pathstr, [name, '.ine']);
            
            % load .frc file
            frcfile = fullfile (pathstr, [name, '.frc']);
            
            self.nodeNames = fieldnames (self.nodes);
            
            self.resultsLoaded = true;
            
        end
        
        function clear (self)
            self.nodes = [];
            self.nNodes = [];
        end
        
        function [hfig, hax] = plotNodeTrajectories (self)
           
            if ~self.resultsLoaded
                error ('No results have been loaded yet for plotting')
            end
            
            hfig = figure;
            hax = axes;
            
            hold on;
            
            for ind = 1:numel (self.nodeNames)
                
                plot3 ( self.nodes.(self.nodeNames{ind}).Position(:,1), ...
                        self.nodes.(self.nodeNames{ind}).Position(:,2), ...
                        self.nodes.(self.nodeNames{ind}).Position(:,3), ...
                        '-' );
            
            end
            
            hold off;
            
            axis equal;
%             axis square;
            
            view(3);
            drawnow;
            
        end
        
        
        function animate (self, varargin)
            % animate the nodes and bodies of the system
            
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
            
            options = parse_pv_pairs (options, varargin);
            
            if self.resultsLoaded == false
                error ('No results have been loaded yet');
            end
            
            plotdata.HAx = [];
            
            for tind = 1:options.Skip:size(self.nodes.(self.nodeNames{1}).Position,1)
                
                if tind > 1
                    % clear old node label drawings
%                     delete (hnode);
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
                              'PlotAxes', plotdata.HAx);
                
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
                
                if ~isempty(options.VideoFile)

                    if tind == 1
                        FPS = options.VideoSpeed ...
                               * numel (1:options.Skip:size(self.nodes.(self.nodeNames{1}).Position,1)) ...
                                        / (self.preProcSystem.problems{1}.finalTime - self.preProcSystem.problems{1}.initialTime);
                        
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
            
            options = parse_pv_pairs (options, varargin);
            
            if isempty (options.PlotAxes)
                % plot the trajectories once to get appropriate axis limits
                [plotdata.HFig, plotdata.HAx] = plotNodeTrajectories (self);
                maximize (plotdata.HFig);
                cla (plotdata.HAx);
            else
                plotdata.HAx = options.PlotAxes;
                plotdata.HFig = get (plotdata.HAx, 'Parent');
            end
            
%             set(hfig,'WindowStyle','docked');
            
            if isempty (options.AxLims)
                stretchfactor = 0.25;
                
                axXlim = get (plotdata.HAx, 'XLim');
                axXlim(1) = axXlim(1) - stretchfactor*abs(axXlim(1));
                axXlim(2) = axXlim(2) + stretchfactor*abs(axXlim(2));
                
                axYlim = get (plotdata.HAx, 'YLim');
                axYlim(1) = axYlim(1) - stretchfactor*abs(axYlim(1));
                axYlim(2) = axYlim(2) + stretchfactor*abs(axYlim(2));
                
                axZlim = get (plotdata.HAx, 'ZLim');
                axZlim(1) = axZlim(1) - stretchfactor*abs(axZlim(1));
                axZlim(2) = axZlim(2) + stretchfactor*abs(axZlim(2));
                
            else
                axXlim = options.AxLims (1,1:2);
                axYlim = options.AxLims (2,1:2);
                axZlim = options.AxLims (3,1:2);
            end
            
            if isempty (options.NodeSize)
                % choose 
                plotdata.nodeSize = repmat (0.03 * mean ( diff ([axXlim; axYlim; axZlim ], 1, 2) ), ...
                             1, self.nNodes);
            else
                if isscalar (options.NodeSize)
                    plotdata.nodeSize = repmat (options.NodeSize, 1, self.nNodes);
                elseif numel (options.NodeSize) ~= self.nNodes
                    error ('You supplied a vector of node sizes which is not the same size as the number of nodes.')
                end
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
                    'StructuralNodes', options.DrawNodes, ...
                    'Joints', false, ...
                    'Light', options.Light );
                
            end
            
            for indii = 1:self.nNodes
                %                     plot3 ( self.nodes.(self.nodeNames{indii}).Position(tind,1), ...
                %                             self.nodes.(self.nodeNames{indii}).Position(tind,2), ...
                %                             self.nodes.(self.nodeNames{indii}).Position(tind,3), ...
                %                             'o' );
                
                
                %                     hnode(indii) = drawnode (self, self.nodes.(self.nodeNames{indii}).Position(tind,:), ...
                %                               -self.nodes.(self.nodeNames{indii}).Orientation(tind,:), ...
                %                               nodesz(indii), nodesz(indii), nodesz(indii) ...
                %                               ... % 0.5, 0.5, 0.5
                %                              );
                
                if options.DrawLabels
                    plotdata.HNodeLabels(indii) ...
                        = text ( plotdata.HAx, ...
                                 self.nodes.(self.nodeNames{indii}).Position(tind,1) + plotdata.nodeSize(indii), ...
                                 self.nodes.(self.nodeNames{indii}).Position(tind,2), ...
                                 self.nodes.(self.nodeNames{indii}).Position(tind,3) + plotdata.nodeSize(indii), ...
                                 self.nodeNames{indii}, ...
                                 'Interpreter', 'none');
                end
                
%                 if options.PlotTrajectories
%                     if tind == 1
%                         htraj(indii) = animatedline ( ...
%                             self.nodes.(self.nodeNames{indii}).Position(tind,1), ...
%                             self.nodes.(self.nodeNames{indii}).Position(tind,2), ...
%                             self.nodes.(self.nodeNames{indii}).Position(tind,3) ...
%                             );
%                         %                             drawnow;
%                         
%                     else
%                         
%                         addpoints ( htraj(indii), ...
%                             self.nodes.(self.nodeNames{indii}).Position(tind,1), ...
%                             self.nodes.(self.nodeNames{indii}).Position(tind,2), ...
%                             self.nodes.(self.nodeNames{indii}).Position(tind,3) );
%                         %                             drawnow;
%                     end
%                 end
                
            end
            
            if ~isempty (self.preProcSystem)
                title (plotdata.HAx, sprintf ('t = %.2fs', (tind-1)*self.preProcSystem.problems{1}.timeStep));
            end
            
            set (plotdata.HAx, 'XLim', axXlim);
            set (plotdata.HAx, 'YLim', axYlim);
            set (plotdata.HAx, 'ZLim', axZlim);
            
            %                 axis equal;
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