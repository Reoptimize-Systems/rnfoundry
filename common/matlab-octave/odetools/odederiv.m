classdef odederiv < handle
    % calculate numerical derivative of quantity for use in ode solver
    % routine
    %
    % Description
    %
    % odederiv is intended for use with the ode solver routines such as
    % ode45, ode15s et. It calulates the numerical derivative of one or
    % more solution components at each time step. 
    %
    % The typical use 
    %
    % odederiv Methods:
    %
    %  derivative - returns the derivative of quantities at a given time
    %  update - updates the previous values of quanties
    %  outputfcn - provides function for use in ode OutputFcn 
    %
    
   properties (GetAccess = public, SetAccess = private)
       ulast;
       tlast;
       dudtlast;
       strict;
       numberOfStepsToFit;
       fastMethod;
   end
   
   properties (GetAccess = private, SetAccess = private)
       doPlot;
       TARGET_FIGURE;
       TARGET_AXIS;
       nt;
       t0;
       u0;
       du0dt0;
   end
   
   methods
       
       function self = odederiv (t0, u0, varargin)
           % constructor for odederiv class
           %
           % Syntax
           %
           % od = odederiv ()
           % od = odederiv ('Parameter', value)
           %
           % Input 
           %
           % Some optional arguments can be supplied as parameter-value
           % pairs. The following options are available:
           %
           %  'Strict' - true or false flag. If true, odederiv will check
           %    it's inputs more carefully for the correct sizes etc.
           %    during simulation and issue more helpful error messages
           %    when something is not right. Defaults to false (for speed).
           %
           %  'DoPlot' - true or false flag. If true, odederiv will plot
           %    the functions and their derivatives at every output step.
           %    Requires that you call the outputfcn method in your ode
           %    OutputFcn. Defaults to false.
           %
           %  'du0dt0' - optional initial value for derivatives at time t0,
           %    if not supplied a vector of zeros will be used.
           %
           %  'UseFastMethod' - flag determining how the derivatives are
           %    computed. 
           % Output
           %
           %   od - a new odederiv class. Note that odederiv derives from
           %     the handle class not value class.
           %
           
           options.Strict = false;
           options.DoPlot = false;
           options.du0dt0 = zeros (size (u0));
           options.NumberOfStepsToFit = 8;
           options.UseFastMethod = true;
           
           options = parse_pv_pairs (options, varargin);
           
           % check inputs
           check.multicheck ( @(x) (isscalar (x) && (isnumeric(x) || islogical (x))), ...
    'Inputs ''t0'' (1) and options ''Strict'' (2), ''DoPlot'' (3), ''NumberOfStepsToFit'' (4) and ''UseFastMethod'' (5) should all be scalar logical or numeric values', ...
               '', ...
               t0, ... % 1
               options.Strict, ... % 2
               options.DoPlot, ... % 3
               options.NumberOfStepsToFit, ... % 4
               options.UseFastMethod ... % 5
                            ); 
           
           assert (isnumeric (options.NumberOfStepsToFit) ...
                    && isint2eps (options.NumberOfStepsToFit) ...
                    && options.NumberOfStepsToFit >= 3, ...
               'If supplied ''NumberOfStepsToFit'' should be an integer >= 3.');
           
           self.strict = options.Strict;
           self.doPlot = options.DoPlot;
           self.numberOfStepsToFit = options.NumberOfStepsToFit;
           self.fastMethod = options.UseFastMethod;
           
           self.t0 = t0;
           self.u0 = u0;
           
           if samesize (options.du0dt0, u0)
               self.du0dt0 = options.du0dt0;
           else
               error ('ODEDERIV:baddu0Vdt0oru0', ...
                   'Supplied du0dt0 was not the same size as u0')
           end
           
           reset (self);
           
       end
       
       function reset (self)
           self.tlast = self.t0 ;
           self.ulast = self.u0;
           self.dudtlast = self.du0dt0;
           self.nt = 0;
       end
       
       function dudt = derivative (self, t, u)
           % get the current value of the derivative
           %
           % call in ode function
           %
           % Syntax
           %
           % odederiv.derivative (t, u)
           %
           % Description
           %
           % derivative returns the numerical derivative of the supplied
           % quantity with respect to it's previous value such that
           %
           % du/dt = (u - prev u) / (t - prev t)
           %
           % derivative is vectorised so u may be a column of several
           % quantities for which the derivative is to be calculated, in
           % this case t is a scalar and u is a column vector of separate
           % quantities all at the point in time.
           %
           % In addition, derivative can take a row vector of n time values
           % and a matrix of values of u, e.g.
           %
           % t = [ t1, t2, t3 ]
           % u = [ u11, u12, u13;
           %       u21, u22, u23;
           %       u31, u32, u33 ]
           %
           % In this case the columns of u are the values of the quantities
           % at each time point in t. derivative will return a matrix
           %
           % dudt = [ d u11 / d t1, d u12 / d t2, d u13 / d t3;
           %          d u21 / d t1, d u22 / d t2, d u23 / d t3;
           %          d u31 / d t1, d u32 / d t2, d u33 / d t3 ]
           %
           % This facilitates vectorisation of ode solution functions
           %
           
           dudt = zeros (size (u));
           
           if self.nt == 0 || t == self.tlast(end)
               dudt = repmat (self.dudtlast(:,end), 1, size (u,2));
           elseif self.fastMethod || self.nt == 1
               % do simple calc
               dudt = bsxfun (@rdivide, bsxfun (@minus, u, self.ulast(:,end) ), (t - self.tlast(end)));
           elseif self.nt == 2
               % use slm
               for compind = 1:size (u,1)       
                   for ucolind = 1:size (u,2)
                       uslm = slmengine ([self.tlast, t], [self.ulast(compind,:), u(compind,ucolind)], ...
                           'knots', 3, ...
                           'LeftSlope', self.dudtlast(compind,1), ...
                           'RightValue', u(compind,ucolind), ...
                           'xy', [self.tlast(2:end)', self.ulast(compind,(2:end))'] );
                       
                       dudt(compind,ucolind) = slmeval (t, uslm, 1, false);
                   end
               end
               
           else
               
               for compind = 1:size (u,1)
                   for ucolind = 1:size (u,2)
                       uslm = slmengine ([self.tlast, t], [self.ulast(compind,:), u(compind,ucolind)], ...
                           'Order', 2, ...
                           'knots', max([3,ceil(numel(self.tlast)/1.5)]), ...
                           ...'LeftSlope', self.dudtlast(compind,1), ...
                           'RightValue', u(compind,ucolind) ...
                           ...'xy', [self.tlast(end-1:end)', self.ulast(compind,(end-1:end))'] ...
                           ... 'xyp', [self.tlast(2:end)', self.dudtlast(compind,(2:end))']
                           );
                       
                       dudt(compind,ucolind) = slmeval (t, uslm, 1, false);
                   end
               end
               
           end
           
       end
       
       function update (self, t, u)
           % update the previous values of t and u
           %
           % call in Outputfcn
           
           if self.strict
               if numel(t) > 1
                   error ('When updating, t must be a scalar');
               end
               if ~samesize(u, self.ulast)
                   error ('When updating, new value of u must be the same size as the previous value');
               end
           end
           
           fprintf (1, 't: %f\n', t);
           
           if self.fastMethod
               self.ulast = u;
               self.tlast = t;
           else
               if self.nt == 0
                   
               elseif self.nt > self.numberOfStepsToFit
                   self.dudtlast  = [self.dudtlast(:,end-self.numberOfStepsToFit:end), derivative(self, t, u)];
                   self.tlast = [self.tlast(end-self.numberOfStepsToFit:end), t];
                   self.ulast = [self.ulast(:,end-self.numberOfStepsToFit:end), u];
                   
               elseif self.nt >= 1
                   self.dudtlast  = [self.dudtlast, derivative(self, t, u)];
                   self.tlast = [self.tlast, t];
                   self.ulast = [self.ulast, u];
                   
               end
           end
           
           % first get derivative at this time step
           %            self.dudtlast = derivative (self, t, u);
           % then update the previous values of t and u
           %            self.ulast = u;
           %            self.tlast = t;
           
           
           self.nt = self.nt + 1;
           
       end
       
       function status = outputfcn (self, t, u, flag, varargin)
           
           status = 0;
           
           callbackDelay = 1;  % Check Stop button every 1 sec.
           
           if isempty (flag)
               
               if self.doPlot
                   
                   dudt = derivative (self, t, u(:,end));
                   
                   y = [u; dudt];
                   
                   if (isempty(self.TARGET_FIGURE) || isempty(self.TARGET_AXIS))
                       
                       error(message('MATLAB:odeplot:NotCalledWithInit'));
                       
                   elseif (ishghandle(self.TARGET_FIGURE) && all(ishghandle(self.TARGET_AXIS)))  % figure still open
                       
                       try
                           ud = get(self.TARGET_FIGURE,'UserData');
                           if ud.stop == 1  % Has stop button been pushed?
                               status = 1;
                           else
                               for i = 1 : length(ud.anim)
                                   addpoints (ud.anim(i), t, y(i,:));
                               end
                               if etime(clock,ud.callbackTime) < callbackDelay
                                   drawnow update;
                               else
                                   ud.callbackTime = clock;
                                   set(self.TARGET_FIGURE,'UserData',ud);
                                   drawnow;
                               end
                           end
                       catch ME
                           error(message('MATLAB:odeplot:ErrorUpdatingWindow', ME.message));
                       end
                   end
               end
               
               self.update (t(end), u(:,end));
               
           elseif strcmp (flag, 'init')
               
               % input is tspan and y0
               self.update (t(1), u(:));
               
               if self.doPlot
                   
                   dudt = derivative (self, t(1), u(:));
                   
                   ncomps = size (dudt, 1);
                   
                   f = figure;
                   self.TARGET_FIGURE = f;
                   ud = get(f,'UserData');
                   
                   % Initialize lines
                   if ~ishold || ~isfield(ud,'lines')
                       [self.TARGET_AXIS,ulines,dudtlines] = plotyy (t(1), u(:), t(1), dudt);
                       self.TARGET_AXIS(1).YLimMode = 'auto';
                       self.TARGET_AXIS(2).YLimMode = 'auto';
                       ud.lines = [ulines, dudtlines];
                   end
                   for i = 1 : ncomps
                       ulines(i).LineStyle = '--';
                       ulines(i).Color = 'b';
                       dudtlines(i).LineStyle = ':';
                       dudtlines(i).Color = 'r';
                       
                       ud.anim(i) = animatedline(t(1),u(i),'Parent',self.TARGET_AXIS(1),...
                           'Color',get(ud.lines(i),'Color'),...
                           'Marker',get(ud.lines(i),'Marker'));
                   end
                   for i = 1 : ncomps
                       ud.anim(i+ncomps) = animatedline(t(1),dudt(i),'Parent',self.TARGET_AXIS(2),...
                           'Color',get(ud.lines(i+ncomps),'Color'),...
                           'Marker',get(ud.lines(i+ncomps),'Marker'));
                   end
                   
                   if ~ishold
                       set(self.TARGET_AXIS(1),'XLim',[min(t) max(t)]);
                       set(self.TARGET_AXIS(2),'XLim',[min(t) max(t)]);
                   end
                   
                   % The STOP button
                   h = findobj(f,'Tag','stop');
                   if isempty(h)
                       pos = get(0,'DefaultUicontrolPosition');
                       pos(1) = pos(1) - 15;
                       pos(2) = pos(2) - 15;
                       uicontrol( ...
                           'Style','pushbutton', ...
                           'String',getString(message('MATLAB:odeplot:ButtonStop')), ...
                           'Position',pos, ...
                           'Callback',@StopButtonCallback, ...
                           'Tag','stop');
                       ud.stop = 0;
                   else
                       % make sure it's visible
                       set(h,'Visible','on');
                       % don't change old ud.stop status
                       if ~ishold || ~isfield(ud,'stop')
                           ud.stop = 0;
                       end
                   end
                   
                   % Set figure data
                   ud.callbackTime = clock;
                   set(f,'UserData',ud);
                   
                   % fast update
                   drawnow update;
               end
               
               
           elseif strcmp (flag, 'done')
               
               if self.doPlot
                   f = self.TARGET_FIGURE;
                   self.TARGET_FIGURE = [];
                   ta = self.TARGET_AXIS;
                   self.TARGET_AXIS = [];
                   
                   if ishghandle(f)
                       ud = get(f,'UserData');
                       if ishghandle(ta)
                           ncomps = size (self.ulast, 1);
                           
                           for i = 1 : size (self.ulast, 1)
                               [tt,yy] = getpoints(ud.anim(i));
                               np = get(ta(1),'NextPlot');
                               set(ta(1),'NextPlot','add');
                               ud.lines(i) = plot(tt,yy,'Parent',ta(1),...
                                   'Color',get(ud.anim(i),'Color'),...
                                   'Marker',get(ud.anim(i),'Marker'));
                               set(ta(1),'NextPlot',np);
                               delete(ud.anim(i));
                               
                               [tt,yy] = getpoints(ud.anim(i+ncomps));
                               np = get(ta(2),'NextPlot');
                               set(ta(2),'NextPlot','add');
                               ud.lines(i+ncomps) = plot(tt,yy,'Parent',ta(2),...
                                   'Color',get(ud.anim(i+ncomps),'Color'),...
                                   'Marker',get(ud.anim(i+ncomps),'Marker'));
                               set(ta(2),'NextPlot',np);
                               delete(ud.anim(i+ncomps));
                           end
                       end
                       set(f,'UserData',rmfield(ud,{'anim','callbackTime'}));
                       if ~ishold
                           set(findobj(f,'Tag','stop'),'Visible','off');
                           
                           if ishghandle(ta(1))
                               set(ta(1),'XLimMode','auto');
                           end
                           if ishghandle(ta(2))
                               set(ta(2),'XLimMode','auto');
                           end
                       end
                   end
                   
                   % full update
                   drawnow;
               end
               
           end
           
       end
       
   end
    
end

% Test scraps
%
% od = odederiv (0, [1;2;3])
% 
% od.derivative (1, [4;4;4])
% 
% ans =
% 
%      3.0000e+000
%      2.0000e+000
%      1.0000e+000
% 
% od.derivative ([1, 2], [4, 5; 4, 5; 4, 5])
% 
% ans =
% 
%      3.0000e+000     2.0000e+000
%      2.0000e+000     1.5000e+000
%      1.0000e+000     1.0000e+000
% 
% od.derivative ([1, 2], [4, 4; 4, 4; 4, 4])
% 
% ans =
% 
%      3.0000e+000     1.5000e+000
%      2.0000e+000     1.0000e+000
%      1.0000e+000   500.0000e-003

