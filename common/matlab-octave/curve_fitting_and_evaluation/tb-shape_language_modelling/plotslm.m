function hfig = plotslm(slm, plotstyle, varargin)
% plotslm: plots a Shape Language Model (slm) or its derivatives
% usage: plotslm(slm)
% usage: hfig = plotslm(slm)
%
%
% arguments: (input)
%  slm - slm struct, as returned by slmengine or slmfit, or evaluated by 
%    slmeval
%
%  plotstyle - char or cellarrays of strings determining the style of plot
%    to be made
%
% Additional arguments may be supplied a Parameter-Value pairs
%
%  'Title' - string containing the figure title, default is 'Function plot'
%  'XLabel' - string containing the label for the first x-axis, default is 
%      'f(X)'
%  'YLabel' - string containing the y-axis label, default is 'X'
%
% arguments: (output)
%  hfig - figure handle. Only returned if an output argument
%        is requested
%
%
% Example:
%  Plot of a fit to exponential data
%  x = 0:.1:1;
%  y = exp(x) + randn(size(x))/10;
%  slm = slmengine(x,y);
%  plotslm(slm)
%
%
% See also: slmset, slmengine, slmeval, ppval, slmfit
%
%
% Author: John D'Errico
% E-mail: woodchips@rochester.rr.com
% Release: 1.0
% Release date: 2/6/07

% ------------------           licence           --------------------
% Copyright (c) 2007, John D'Errico
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
% ------------------           licence           --------------------

if nargin < 2
    plotstyle = {''};
elseif ~iscell(plotstyle)
    plotstyle = {plotstyle};
end

options.Title = 'Function plot';
options.XLabel = 'X';
options.YLabel = 'f(X)';
options.ShowKnots = true;
options.ShowData = true;
options.DYLabel = 'dy';
options.DY2Label = 'dy2';
options.DY3Label = 'dy3';

options = parse_pv_pairs (options, varargin);

% initial plot parameters
params.LineStyle = '-';
params.LineColor = [0.4660, 0.6740, 0.1880]; %'r';
params.LineWidth = 0.5;

params.ShowData = options.ShowData;
params.DataMarker = 'o';
params.DataColor = [0, 0.4470, 0.7410]; %'b';
params.DataSize = 3;

params.ShowKnots = options.ShowKnots;
params.KnotStyle = '--';
params.KnotColor = 'g';

params.Grid = 'off';

params.NumberOfPoints = 1001;

params.Title = options.Title;
params.XLabel = options.XLabel;
params.YLabel = options.YLabel;

% is it a slm or pp or old style hermite form?
if isstruct(slm)
  % this covers both the pp or slm forms
  params.model = slm;
  
  if strcmp(slm.form,'slm')
    params.Degree = slm.degree;
  else
    % its a pp
    params.Degree = slm.order - 1;
  end
elseif isa(slm,'double') && ismember(size(slm,2),[2 3])
  % It must have been a linear or cubic Hermite array form.
  % Be friendly, and convert to a slm.
  params.model = hermite2slm(slm);
  
  % define linear as first order, cubic is third
  params.Degree = params.model.degree;
end

% do we have any data?
if isfield(params.model,'x') && ~isempty(params.model.x)
  params.DataExists = true;
else
  params.DataExists = false;
end

% which style will we use initially
if params.DataExists
  params.PlotStyle = 'fundata';
else
  params.PlotStyle = 'fun';
end

% generate plot figure
figuregen

% and plot the model. The options will be
% 'fundata', 'fun', 'residuals', 'dy', 'dy2', 'dy3', 'integral'
if isempty(plotstyle)
    plotcurve(params.PlotStyle)
else
    plotcurve([params.PlotStyle, plotstyle]);
end

% any further interactions will be by callbacks.

% do we return the figure handle?
if nargout > 0
  hfig = params.handles.figure;
end

% ===============================================
%    nested functions - initial plot figure
% ===============================================
function figuregen
  % generates/initializes the plot gui
  
  % a new figure
  params.handles.figure = figure;
  
  % generate the plot menu
  h=uimenu(params.handles.figure,'Label','Plot');
  uimenu(h,'Label','Curve only','Callback',@(s,e) plotcurve('fun'));
  if params.DataExists
    uimenu(h,'Label','Curve & Data','Callback',@(s,e) plotcurve('fundata'));
    uimenu(h,'Label','Residuals','Callback',@(s,e) plotcurve('residuals'));
  end
  if params.Degree > 0
    uimenu(h,'Label','First derivative','Callback',@(s,e) plotcurve('dy'));
  end
  if params.Degree > 1
    uimenu(h,'Label','Second derivative','Callback',@(s,e) plotcurve('dy2'));
    uimenu(h,'Label','Third derivative','Callback',@(s,e) plotcurve('dy3'));
  end
  uimenu(h,'Label','Integral','Callback',@(s,e) plotcurve('integral'));
  
  
  if params.DataExists
    hh=uimenu(h,'Label','Data points','Separator','on');
    hhh=uimenu(hh,'Label','Data markers');
    uimenu(hhh,'Label','o','Callback',@(s,e) setplotparam('DataMarker','o'));
    uimenu(hhh,'Label','+','Callback',@(s,e) setplotparam('DataMarker','+'));
    uimenu(hhh,'Label','x','Callback',@(s,e) setplotparam('DataMarker','x'));
    uimenu(hhh,'Label','.','Callback',@(s,e) setplotparam('DataMarker','.'));
    uimenu(hhh,'Label','Triangle','Callback',@(s,e) setplotparam('DataMarker','v'));
    uimenu(hhh,'Label','Square','Callback',@(s,e) setplotparam('DataMarker','s'));
    uimenu(hhh,'Label','Pentagram','Callback',@(s,e) setplotparam('DataMarker','p'));
    uimenu(hhh,'Label','Hexagram','Callback',@(s,e) setplotparam('DataMarker','h'));
    uimenu(hhh,'Label','None','Callback',@(s,e) setplotparam('DataMarker',''));
    
    hhh=uimenu(hh,'Label','Symbol color');
    uimenu(hhh,'Label','Red','Callback',@(s,e) setplotparam('DataColor','r'));
    uimenu(hhh,'Label','Green','Callback',@(s,e) setplotparam('DataColor','g'));
    uimenu(hhh,'Label','Blue','Callback',@(s,e) setplotparam('DataColor','b'));
    uimenu(hhh,'Label','Cyan','Callback',@(s,e) setplotparam('DataColor','c'));
    uimenu(hhh,'Label','Magenta','Callback',@(s,e) setplotparam('DataColor','m'));
    uimenu(hhh,'Label','Yellow','Callback',@(s,e) setplotparam('DataColor','y'));
    uimenu(hhh,'Label','Black','Callback',@(s,e) setplotparam('DataColor','k'));
    
    hhh=uimenu(hh,'Label','Symbol size');
    uimenu(hhh,'Label','1','Callback',@(s,e) setplotparam('DataSize',1));
    uimenu(hhh,'Label','2','Callback',@(s,e) setplotparam('DataSize',2));
    uimenu(hhh,'Label','3','Callback',@(s,e) setplotparam('DataSize',3));
    uimenu(hhh,'Label','4','Callback',@(s,e) setplotparam('DataSize',4));
    uimenu(hhh,'Label','6','Callback',@(s,e) setplotparam('DataSize',6));
    uimenu(hhh,'Label','10','Callback',@(s,e) setplotparam('DataSize',10));
    uimenu(hhh,'Label','20','Callback',@(s,e) setplotparam('DataSize',20));
  end
  
  hh=uimenu(h,'Label','Knots/Breaks');
  hhh=uimenu(hh,'Label','Line style');
  uimenu(hhh,'Label','Solid','Callback',@(s,e) setplotparam('KnotStyle','-'));
  uimenu(hhh,'Label','Dotted','Callback',@(s,e) setplotparam('KnotStyle',':'));
  uimenu(hhh,'Label','Dash/dot','Callback',@(s,e) setplotparam('KnotStyle','-.'));
  uimenu(hhh,'Label','Dashed','Callback',@(s,e) setplotparam('KnotStyle','--'));
  uimenu(hhh,'Label','None','Callback',@(s,e) setplotparam('KnotStyle','none'));
  
  hhh=uimenu(hh,'Label','Line color');
  uimenu(hhh,'Label','Red','Callback',@(s,e) setplotparam('KnotColor','r'));
  uimenu(hhh,'Label','Green','Callback',@(s,e) setplotparam('KnotColor','g'));
  uimenu(hhh,'Label','Blue','Callback',@(s,e) setplotparam('KnotColor','b'));
  uimenu(hhh,'Label','Cyan','Callback',@(s,e) setplotparam('KnotColor','c'));
  uimenu(hhh,'Label','Magenta','Callback',@(s,e) setplotparam('KnotColor','m'));
  uimenu(hhh,'Label','Yellow','Callback',@(s,e) setplotparam('KnotColor','y'));
  uimenu(hhh,'Label','Black','Callback',@(s,e) setplotparam('KnotColor','k'));
  
  hh=uimenu(h,'Label','Curve');
  hhh=uimenu(hh,'Label','Line style');
  uimenu(hhh,'Label','Solid','Callback',@(s,e) setplotparam('LineStyle','-'));
  uimenu(hhh,'Label','Dotted','Callback',@(s,e) setplotparam('LineStyle',':'));
  uimenu(hhh,'Label','Dash/dot','Callback',@(s,e) setplotparam('LineStyle','-.'));
  uimenu(hhh,'Label','Dashed','Callback',@(s,e) setplotparam('LineStyle','--'));
  uimenu(hhh,'Label','None','Callback',@(s,e) setplotparam('LineStyle','none'));
  
  hhh=uimenu(hh,'Label','Line color');
  uimenu(hhh,'Label','Red','Callback',@(s,e) setplotparam('LineColor','r'));
  uimenu(hhh,'Label','Green','Callback',@(s,e) setplotparam('LineColor','g'));
  uimenu(hhh,'Label','Blue','Callback',@(s,e) setplotparam('LineColor','b'));
  uimenu(hhh,'Label','Cyan','Callback',@(s,e) setplotparam('LineColor','c'));
  uimenu(hhh,'Label','Magenta','Callback',@(s,e) setplotparam('LineColor','m'));
  uimenu(hhh,'Label','Yellow','Callback',@(s,e) setplotparam('LineColor','y'));
  uimenu(hhh,'Label','Black','Callback',@(s,e) setplotparam('LineColor','k'));
  
  hh=uimenu(h,'Label','Interrogate data by mouse','Separator','on');
  uimenu(hh,'Label','Label a point by mouse','Callback',@(s,e) selectp('remove'));
  uimenu(hh,'Label','Permanent Label','Callback',@(s,e) selectp('permanent'));
  
  hh=uimenu(h,'Label','Toggle grid','Separator','on','Callback',@(s,e) togglegrid);
  
  
end % figuregen
  
% ===============================================
%    nested functions - set plot parameters
% ===============================================
function setplotparam(parameter,val)
% sets plot parameters from callbacks
switch parameter
  case 'DataMarker'
    params.DataMarker = val;
  case 'DataColor'
    params.DataColor = val;
  case 'DataSize'
    params.DataSize = val;

  case 'KnotStyle'
    params.KnotStyle = val;
  case 'KnotColor'
    params.KnotColor = val;
    
  case 'LineStyle'
    params.LineStyle = val;
  case 'LineColor'
    params.LineColor = val;
    
end

% and regenerate the figure to reflect the change
plotcurve(params.PlotStyle)

end


% ===============================================
%    nested functions - plots function/data/derivatives
% ===============================================
function selectp(action)
  % uses selectdata to flag a point
  
  switch action
    case 'remove'
      removelabel = 'on';
    case 'permanent'
      removelabel = 'off';
  end
  
  selectdata('selectionmode','c','label','on','removelabel',removelabel);
  
end % 

% ===============================================
%    nested functions - plots function/data/derivatives
% ===============================================
function togglegrid
  
  if strcmp(params.Grid,'on')
    params.Grid = 'off';
    grid off
  else
    params.Grid = 'on';
    grid on
  end
  
end % 

% ===============================================
%    nested functions - plots function/data/derivatives
% ===============================================
function plotcurve(plotstyle)
% plots the data and model. The options will be
% 'fundata', 'fun', 'residuals', 'dy', 'dy2', 'dy3', 'integral'

% make sure the figure is on top
figure(params.handles.figure)

% if a single string has been passed in, convert to a cell array of size 1
if ~iscell(plotstyle)
    plotstyle = {plotstyle};
end

% We will check if we need to plot the knots during the loop
knotsplotted = false;

if any(strcmp(plotstyle{1}, {'fundata' 'fun'}))
    
    % plot the function itself. we will add the
    % data later if appropriate
    if strcmpi(params.model.form,'slm')
        xrange = params.model.knots([1,end]);
    else
        xrange = params.model.breaks([1,end]);
    end

    xev = linspace(xrange(1),xrange(2),params.NumberOfPoints);

    % evaluate
    if strcmpi(params.model.form,'slm')
        ypred = slmeval(xev,params.model);
    else
        ypred = ppval(params.model,xev);
    end

    % plot the curve
    h = plot(xev,ypred);
    set(h,'Marker','none', ...
        'Color',params.LineColor, ...
        'LineStyle',params.LineStyle, ...
        'LineWidth',params.LineWidth)
    
    if ~knotsplotted && params.ShowKnots

        hold on

        % plot the vertical knot marker lines
        axlim = axis;
        yrange = axlim(3:4);

        if strcmpi(params.model.form,'slm')
            knots = params.model.knots(:);
        else
            knots = params.model.breaks(:);
        end

        h = plot(repmat(knots',2,1),yrange(:));
        set(h,'Marker','none', ...
            'Color',params.KnotColor, ...
            'LineStyle',params.KnotStyle)

        % release the hold
        hold off

        knotsplotted = true;

    end
    
    xlabel (params.XLabel)
    ylabel (params.YLabel)
    title  (params.Title)
    
end

% now for each plot requested in the cell array, do the plot
for i = 2:numel(plotstyle)

    % which plot was called for?
    switch plotstyle{i}

        case 'residuals'
            % only a residual plot, no function at all

            if isfield(params.model,'x') && ~isempty(params.model.x)
                % get the residuals
                if strcmpi(params.model.form,'slm')
                    % a slm form
                    yhat = slmeval(params.model.x,params.model);
                else
                    % it was a pp form
                    yhat = ppval(params.model,params.model.x);
                end
                resid = yhat - params.model.y;

                % plot them
                hplts = addaxis(params.model.x,resid);

                addaxislabel(i, 'Residual errors f(X)');

                setaddaxisprops(i, 'YColor', params.DataColor);
                
                set(hplts, 'Marker', '+', ...
                           'Color', params.DataColor, ...
                           'MarkerSize', params.DataSize, ...
                           'LineStyle', 'none')

            end

        case {'dy' 'dy2' 'dy3'}

            % plot the k'th derivative
            switch plotstyle{i}
                case 'dy'
                    k = 1;
                case 'dy2'
                    k = 2;
                case 'dy3'
                    k = 3;
            end

            if strcmpi(params.model.form,'slm')
                xrange = params.model.knots([1,end]);
            else
                xrange = params.model.breaks([1,end]);
            end

            xev = linspace(xrange(1),xrange(2),params.NumberOfPoints);

            % evaluate
            if strcmpi(params.model.form,'slm')
                ypred = slmeval(xev,params.model,k);
            else
                pp = params.model;
                % differentiate the model
                order = size(pp.coefs,2);
                for j = 1:k
                    pp.coefs = pp.coefs*diag((order-1):-1:1,1);
                end
                ypred = ppval(pp,xev);
            end

            % plot the curve
            h = addaxis(xev,ypred);

            derivLineColor = [0.8500, 0.3250, 0.0980]; %'m';
            
            setaddaxisprops(i, 'YColor', derivLineColor);
            
            switch plotstyle{i}
                case 'dy'
                    addaxislabel(i,options.DYLabel);
                    derivLineStyle = '--';
                case 'dy2'
                    addaxislabel(i,options.DY2Label)
                    derivLineStyle = ':';
                case 'dy3'
                    addaxislabel(i,options.DY3Label)
                    derivLineStyle = '-.';
            end

            set(h,'Marker','none', ...
                   'Color', derivLineColor, ...
                   'LineStyle', derivLineStyle, ...
                   'LineWidth', params.LineWidth)


        case 'integral'

            % plot the cumulative integral

            % if its a slm, easiest is to convert to a
            % pp form, then integrate that
            if strcmpi(params.model.form,'slm')
                pp = slm2pp(params.model);
            else
                % it already was a pp
                pp = params.model;
            end

            xrange = pp.breaks([1,end]);
            xev = linspace(xrange(1),xrange(2),params.NumberOfPoints);

            % integrate (once)
            order = size(pp.coefs,2);
            op = diag(1./(order:-1:1),-1);
            op(1,:) = [];
            pp.coefs = pp.coefs*op;
            pp.order = pp.order + 1;
            % offset the constant terms to make the integral cumulative
            % Just use horner's rule to evaluate at the upper endpoint
            % of each interval
            knots = pp.breaks';
            dx = diff(knots); %#ok
            offset = pp.coefs(:,1);
            for j = 2:(order+1)
                offset = offset.*dx + pp.coefs(:,j);
            end
            pp.coefs(2:end,end) = pp.coefs(2:end,end) + ...
                cumsum(offset(1:(end-1)));

            % evaluate
            ypred = ppval(pp,xev);

            % plot the curve
            h = plot(xev,ypred);

            addaxislabel(i, 'Cumulative integral');

            set(h,'Marker','none', ...
                'Color',params.LineColor, ...
                'LineStyle',params.LineStyle, ...
                'LineWidth',params.LineWidth);

    end

end % end plotstyle cell array loop


if any(strcmp('fundata', plotstyle))
  % we must also plot the data
  if params.ShowData ...
          && isfield(params.model,'x') ...
          && ~isempty(params.model.x)

    h = addaxisplot(params.model.x,params.model.y,1);
    
    set(h,'LineStyle','none', ...
          'Marker',params.DataMarker, ...
          'Color',params.DataColor, ...
          'MarkerSize',params.DataSize)
%     hold off
  end
  
  if isfield(params.model,'prescription') && ~isempty(params.model.prescription.ErrorBar)
    hold on
    errorbar(params.model.x,params.model.y, ...
      params.model.prescription.ErrorBar(:,1), ...
      params.model.prescription.ErrorBar(:,2), ...
      'LineStyle','none', ...
      'Marker',params.DataMarker, ...
      'Color',params.DataColor, ...
      'MarkerSize',params.DataSize)
    hold off
  end
end

% do we need the grid plotted?
if strcmp(params.Grid,'on')
  grid on
end  

end % plotcurve

% ===============================================
%    end mainline
% ===============================================
end

