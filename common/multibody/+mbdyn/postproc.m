classdef postproc < handle
   
    properties 
        
        preProcSystem;
        resultsLoaded;
        nNodes;
        nodes;
        bodies;
        joints;
        
    end
    
    methods
        
        function self = postproc (preprocsys)
            
            if nargin > 0
                if isa (preprocsys, 'mbdyn.pre.system')
                    self.preProcSystem = preprocsys;
                else
                    error ('preprocsys must be a mbdyn.pre.system object')
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
            
            self.resultsLoaded = true;
            
        end
        
        function clear (self)
            self.nodes = [];
            self.nNodes = [];
        end
        
        function [hfig, hax] = plotNodeTrajectories (self)
           
            nodenames = fieldnames (self.nodes);
            
            hfig = figure;
            hax = axes;
            
            hold on;
            
            for ind = 1:numel (nodenames)
                
                plot3 ( self.nodes.(nodenames{ind}).Position(:,1), ...
                        self.nodes.(nodenames{ind}).Position(:,2), ...
                        self.nodes.(nodenames{ind}).Position(:,3), ...
                        '-' );
            
            end
            
            hold off;
            
            axis equal;
            
            view(3);
            
        end
        
        
        function animate (self, varargin)
            
            options.DrawLabels = false;
            options.AxLims = [];
            options.PlotTrajectories = false;
            options.NodeSize = [];
            options.DrawMode = 'wireghost';
            options.DrawNodes = true;
            options.DrawBodies = true;
            options.Skip = 1;
            
            options = parse_pv_pairs (options, varargin);
            
            if self.resultsLoaded == false
                error ('No results have been loaded yet');
            end
            
            % plot the trajectories once to get appropriate axis limits
            [hfig, hax] = plotNodeTrajectories (self);
            
            set(hfig,'WindowStyle','docked');
            
            if isempty (options.AxLims)
                stretchfactor = 0.25;
                
                axXlim = get (hax, 'XLim');
                axXlim(1) = axXlim(1) - stretchfactor*abs(axXlim(1));
                axXlim(2) = axXlim(2) + stretchfactor*abs(axXlim(2));
                
                axYlim = get (hax, 'YLim');
                axYlim(1) = axYlim(1) - stretchfactor*abs(axYlim(1));
                axYlim(2) = axYlim(2) + stretchfactor*abs(axYlim(2));
                
                axZlim = get (hax, 'ZLim');
                axZlim(1) = axZlim(1) - stretchfactor*abs(axZlim(1));
                axZlim(2) = axZlim(2) + stretchfactor*abs(axZlim(2));
            else
                axXlim = options.AxLims (1,1:2);
                axYlim = options.AxLims (2,1:2);
                axZlim = options.AxLims (3,1:2);
            end
            
            nodenames = fieldnames (self.nodes);
            
            nnodes = numel (nodenames);
            
            if isempty (options.NodeSize)
                % choose 
                nodesz = repmat (0.03 * mean ( diff ([axXlim; axYlim; axZlim ], 1, 2) ), ...
                             1, nnodes);
            else
                if isscalar (options.NodeSize)
                    nodesz = repmat (options.NodeSize, 1, nnodes);
                elseif numel (options.NodeSize) ~= nnodes
                    error ('You supplied a vector of node sizes which is not the same size as the number of nodes.')
                end
            end
            
            % clear the previous trajectory plot
            cla;
            
            % TODO: the sped of all this could be improved using hgtransform
            %
            % see https://uk.mathworks.com/help/matlab/creating_plots/transform-objects.html
            %
            for tind = 1:options.Skip:size(self.nodes.(nodenames{1}).Position,1)
                
                if tind > 1
                    % clear old node drawings
%                     delete (hnode);
                    if options.DrawLabels
                        delete (hlabel);
                    end
                end
                    
%                 hold on;

                if ~isempty (self.preProcSystem)
                    
                    for indii = 1:nnodes
                        % set the node positions and orientations
                        self.preProcSystem.setNodePosition (self.nodes.(nodenames{indii}).Label, self.nodes.(nodenames{indii}).Position(tind,:)');
                        om = mbdyn.pre.orientmat ('orientation', reshape (self.nodes.(nodenames{indii}).Orientation(tind,:), 3, 3));
                        self.preProcSystem.setNodeOrientation (self.nodes.(nodenames{indii}).Label, om);
                    end

                    % draw the system 
                    self.preProcSystem.draw ( 'AxesHandle', hax, ...
                                              'Mode', options.DrawMode, ...
                                              'Bodies', options.DrawBodies, ...
                                              'StructuralNodes', options.DrawNodes, ...
                                              'Joints', false );
                                      
                end
                
                for indii = 1:nnodes    
%                     plot3 ( self.nodes.(nodenames{indii}).Position(tind,1), ...
%                             self.nodes.(nodenames{indii}).Position(tind,2), ...
%                             self.nodes.(nodenames{indii}).Position(tind,3), ...
%                             'o' );

                    
%                     hnode(indii) = drawnode (self, self.nodes.(nodenames{indii}).Position(tind,:), ...
%                               -self.nodes.(nodenames{indii}).Orientation(tind,:), ...
%                               nodesz(indii), nodesz(indii), nodesz(indii) ...
%                               ... % 0.5, 0.5, 0.5 
%                              );
                        
                    if options.DrawLabels
                        hlabel(indii) = text ( self.nodes.(nodenames{indii}).Position(tind,1) + nodesz(indii), ...
                                self.nodes.(nodenames{indii}).Position(tind,2), ...
                                self.nodes.(nodenames{indii}).Position(tind,3) + nodesz(indii), ...
                                nodenames{indii}, ...
                                'Interpreter', 'none');
                    end
                    
                    if options.PlotTrajectories
                        if tind == 1
                            htraj(indii) = animatedline ( ...
                                        self.nodes.(nodenames{indii}).Position(tind,1), ...
                                        self.nodes.(nodenames{indii}).Position(tind,2), ...
                                        self.nodes.(nodenames{indii}).Position(tind,3) ...
                                                        );
%                             drawnow;
                            
                        else
                            
                            addpoints ( htraj(indii), ...
                                        self.nodes.(nodenames{indii}).Position(tind,1), ...
                                        self.nodes.(nodenames{indii}).Position(tind,2), ...
                                        self.nodes.(nodenames{indii}).Position(tind,3) );
%                             drawnow;
                        end
                    end

                end
                
%                 hold off;

                title (hax, sprintf ('t = %.2fs', (tind-1)*self.preProcSystem.problems{1}.timeStep));
                
                set (hax, 'XLim', axXlim);
                set (hax, 'YLim', axYlim);
                set (hax, 'ZLim', axZlim);
                
%                 axis equal;
                view(3);
                
                drawnow;
%                 pause (1e-4)
            end
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