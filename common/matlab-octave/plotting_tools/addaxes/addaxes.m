function HA = addaxes(varargin)
%ADDAXES   Adds a new linked axis related by any monotonic function.
%
%   SYNTAX:
%          addaxes(...,'PropertyName',PropertyValue)
%          addaxes('off')
%          addaxes(AX,...)
%     HA = addaxes(...);
%
%   INPUT:
%     AX      - Uses given axes handle instead of current one.
%               DEFAULT: gca
%     'PN'/PV - Paired property/value inputs to define the new axes, at
%               least one of the following (see NOTE below):
%
%                --------------|---------------------|--------------------
%                   'NAME'             VALUE                DEFAULT
%                --------------|---------------------|--------------------
%                 'XFunction'     Function name or     @(x)interp1(...
%                     or              handle            get(AX,'XLim'),...
%                   'xfun'                              get(AX,'XLim'),...
%                                                       x,'linear',...,
%                                                         'extrap')
%
%                 'XInverse'     Function name or       Inverse function
%                     or              handle          of 'xfun' estimated
%                   'xinv'                                with FZERO.
%
%                 'XLegend'        XLABEL legend          '' (none)
%                     or
%                   'xleg'
%
%                   'XDate'     true, false, 'tlabel'    false (TLABEL not 
%                     or           or 'datetick'            used)
%                   'xdat'
%                --------------|---------------------|--------------------
%
%     'off'   - Deletes any axes previously added. See NOTE below.
%
%   OUTPUT:
%     HA - Handle of new axes. Not recommended to modify it. See NOTE
%          below.
%
%   DESCRIPTION:
%     This functions adds axes to an existent axes, and links them to work
%     with ZOOM/PAN and DATETICK (or TLABEL) if required.
%
%     The big difference with the PLOTYY function (besides of working with
%     the x-axis as well) is that the added axis is not forced to be
%     linearly related with the old one. See the EXAMPLE below.
%
%     The idea behind this function, is the creation of temporal invisible
%     axes to get the Ticks and TickLabels acoordingly to the specifyed
%     relationship. Then other axes are drawn above the current one with
%     the same limits but with the Ticks mapped by inverting the
%     relationships using the FZERO formulae or the given inverse function.
%     See NOTE below for some warnings.
%
%   NOTE:
%     * Optional inputs use its DEFAULT value when not given or [].
%     * Optional outputs may or not be called.
%     * Besides from the paired optional inputs in Table above, the
%       following normal axes properties may be used: 
%         'XColor'  'XTick'  'XTickLabel'
%       for  axis customization and 
%         'CLim'        'FontAngle'   'FontName'  'FontSize'
%         'FontUnits'   'FontWeight'  'Layer'     'LineWidth'   
%         'TickDir'     'TickLength' 
%       for axes customization.
%     * Of course, 'Y' axis properties as 'YFunction', 'YColor', etc., may
%       be used as well, but not for 'Z'.
%     * 'xFun', 'yFun', 'xInv' and 'yInv' must return an array of the same
%       size as the input and must be monotonically increasing or
%       decreasing.
%     * When the relationships 'xFun' and 'yFun' are not linear and given
%       by an interpolation (INTERP1 for example) the user must use them
%       carefully taking into account extrapolations and the approximations
%       results of the invertion formulae provided by FZERO.
%     * If the new axis are dates, use the 'XDate'or 'YDate' options like
%       in the following EXAMPLE.
%     * By now, the function only can be used to add another axes once.
%       That is, multiple axes cannot be added, and to include both x- and
%       y-axis give both 'xFun' and 'yFun' functions.
%
%   EXAMPLE:
%     % DATA
%      t = sort(unique(round(rand(50,1)*100)/100),'descend');
%      [X,Y,Z] = peaks(length(t));
%      xFun = @(x)interp1(X(1,:)',exp(t),x,'linear','extrap');  % Numerical
%      yFun = @(x)x.*abs(x)+datenum(date);                      % Cuadratic
%      xInv = @(x)interp1(exp(t),X(1,:)',x,'linear','extrap');  % Numerical
%     % PLOT
%      figure
%       imagesc(X(1,[1 end]),Y([1 end],1),Z), set(gca,'Layer','top')
%       addaxes(...                           % <- THE CLUE
%        'XFun'   ,xFun,...
%        'YFun'   ,yFun,...
%        'XInv'   ,xInv,...
%        'XLeg'   ,func2str(xFun),...
%        'YDat'   ,'tlabel',...     % date ticks!
%        'XColor' ,'b',...
%        'YColor' ,'r');
%        zoom on
%     
%   SEE ALSO:
%     AXES, PLOT, PLOTYY, FUNCTION_HANDLE, FZERO, DATETICK
%     and
%     TLABEL by Carlos Vargas
%     at http://www.mathworks.com/matlabcentral/fileexchange
%
%
%   ---
%   MFILE:   addaxes.m
%   VERSION: 1.1 (Sep dd, 2009) (<a href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/authors/11258')">download</a>) 
%   MATLAB:  7.7.0.471 (R2008b)
%   AUTHOR:  Carlos Adrian Vargas Aguilera (MEXICO)
%   CONTACT: nubeobscura@hotmail.com

%   REVISIONS:
%   1.0      Released. (Jul 29, 2009)
%   1.1      Fixed bug with function inversion, thanks to Allen Hall. Added
%            'Inverse' options and 'Position' link. Fixed small bug related
%            with 'off' option. (Sep dd, 2009)

%   DISCLAIMER:
%   addaxes.m is provided "as is" without warranty of any kind, under the
%   revised BSD license.

%   Copyright (c) 2009 Carlos Adrian Vargas Aguilera

% INPUTS CHECK-IN
% -------------------------------------------------------------------------

% Parameters.
myAppName   = 'addAxes'; % Application data name.
zoomAppName = 'zoom_zoomOrigAxesLimits';

% Checks number of inputs and outputs.
if     nargin<1
 error('CVARGAS:addaxes:notEnoughInputs',...
  'At least 1 input is required.')
elseif nargout>1
 error('CVARGAS:addaxes:tooManyOutputs',...
  'At most 1 output is allowed.')
end


% -------------------------------------------------------------------------
% MAIN
% -------------------------------------------------------------------------
 
if ~((nargin==2) && isstruct(varargin{2}))
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % ADDAXES called from command window or a M-file
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % Default.
 AX   = gca;
 xFun = @(x) interp1(get(AX,'XLim'),get(AX,'XLim'),x,'linear','extrap');
 yFun = @(x) interp1(get(AX,'YLim'),get(AX,'YLim'),x,'linear','extrap');
 xInv = [];
 yInv = [];
 xLeg = '';
 yLeg = '';
 dFun = 'tlabel';      % Function to use on date axis.
 nSeg = 0.25;          % Second to wait for double-click.
 
 % Program data. To be saved as application data named appname.
 data.AX   = AX;
 data.HA   = [];
 data.HF   = [];
 data.xNew = 0;
 data.yNew = 0;
 data.xFun = xFun;
 data.yFun = yFun;
 data.xInv = xInv;
 data.yInv = yInv;
 data.xLeg = xLeg;
 data.yLeg = yLeg;
 data.xDat = false;
 data.yDat = false;
 data.aOpt = {};
 data.aNam = myAppName;
 data.nSeg = nSeg;
 
 % Parse inputs.
 [data,xTic,yTic,xTLa,yTLa,flag] = parseInputs(data,dFun,varargin{:});
 clear varargin
 
 % Checks if 'off' option.
 if flag
  data = getappdata(data.AX,myAppName);
  if ~isempty(data)
   if ishandle(data.AX) % Fixed bug, Sep 2009
    zh = zoom(data.AX);
    ph = pan(data.HF);
    set(zh,'ActionPostCallback',[])
    set(ph,'ActionPostCallback',[])
    rmappdata(data.AX,myAppName)
    if ishandle(data.HA)
     linkaxes([data.AX data.HA],'off')
    end
   end
   if ishandle(data.HA)
    delete(data.HA)
   end
  end
  if nargout<0, HA = []; end
  return
 end
 
 % Gets axis parent.
 data.HF = ancestor(data.AX,{'figure','uipanel'});

 % Changes some aspects on old axes.
 view(data.AX,2)
 box (data.AX,'off')
 grid(data.AX,'off')
 
 % Gets new axis locations in front of old one.
 pxLoc = get(data.AX,'XAxisLocation');
 pyLoc = get(data.AX,'YAxisLocation');
 xLoc  = 'top';   if strcmp(pxLoc,xLoc), xLoc = 'bottom'; end
 yLoc  = 'right'; if strcmp(pyLoc,yLoc), yLoc = 'left';   end
 
 % Generates new axes.
 tempF = get(0    ,'CurrentFigure');
 tempA = get(tempF,'CurrentAxes');
 data.HA = axes(...
  'Parent'          ,data.HF,...
  'Color'           ,'none',...
  'Box'             ,'off',...
  'Units'           ,get(data.AX,'Units'),...
  'Position'        ,get(data.AX,'Position'),...
  'View'            ,get(data.AX,'View'),...
  'FontUnits'       ,get(data.AX,'FontUnits'),...
  'FontSize'        ,get(data.AX,'FontSize'),...
  'FontName'        ,get(data.AX,'FontNam'),...
  'FontWeight'      ,get(data.AX,'FontWeight'),...
  'FontAngle'       ,get(data.AX,'FontAngle'),...
  'Layer'           ,get(data.AX,'Layer'),...
  'LineWidth'       ,get(data.AX,'LineWidth'),...
  'Projection'      ,get(data.AX,'Projection'),...
  'TickLength'      ,get(data.AX,'TickLength'),...
  'TickDir'         ,get(data.AX,'TickDir'),...
  'XColor'          ,get(data.AX,'XColor'),...
  'YColor'          ,get(data.AX,'YColor'),...
  'XDir'            ,get(data.AX,'XDir'),...
  'YDir'            ,get(data.AX,'YDir'),...
  'XAxisLocation'   ,xLoc,...
  'YAxisLocation'   ,yLoc,...
  'XTick'           ,[],...
  'YTick'           ,[],...
  'Tag'             ,'addaxes',...
  data.aOpt{:});
 set(0    ,'CurrentFigure',tempF)
 set(tempF,'CurrentAxes'  ,tempA)
 
 % Links properties.
 theLinks = linkprop([data.AX data.HA],{'Units','Position','View'});
 setappdata(data.AX,[data.aNam 'Links'],theLinks);
 setappdata(data.HA,[data.aNam 'Links'],theLinks);
 
 % Fixes title.
 if strcmp(xLoc,'top')
  title(data.HA,get(get(data.AX,'Title'),'String'))
  title(data.AX,'')
  set(data.HA,'Position',get(data.AX,'Position')) % In case it moves
 end
 
 % Updates Ticks's and TickLabel's.
 updateTickAndTickLabel(data,xTic,yTic,xTLa,yTLa)
 
 % Sets ZOOM and PAN functionalities.
 zoom(data.AX,'reset')
 zh = zoom(data.AX);
 ph = pan(data.HF);
 set(zh,'ActionPostCallback',@addaxes)
 set(ph,'ActionPostCallback',@addaxes)

 % Links axes.
 linkaxes([data.AX data.HA])
 
 % Saves data.
 setappdata(data.HA,data.aNam,data)
 setappdata(data.AX,data.aNam,data)
 
else
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % ADDAXES called after ZOOM or PAN
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % Loads data.
 data = getappdata(varargin{2}.Axes,myAppName);
 if isempty(data)
  if nargout==1, HA = data.HA; end
  return
 end
 
 % Checks axes.
 if any(~ishandle([data.AX data.HA]))
  if nargout==1, HA = data.HA; end
  return
 end
 
 % Finds out if double-click.
 pause(data.nSeg)
 drawnow
 doubleclick = strcmp(get(data.HF,'SelectionType'),'open');
 if doubleclick
  axis(data.AX,getappdata(data.AX,zoomAppName))
 end
 
 % Initializes.
 xTic = [];  xTLa = '';
 yTic = [];  yTLa = '';
 
 % Updates Ticks's and TickLabel's.
 updateTickAndTickLabel(data,xTic,yTic,xTLa,yTLa);
 
end
 

% OUTPUTS CHECK-OUT
% -------------------------------------------------------------------------

if nargout==1, HA = data.HA; end


% =========================================================================
% SUBFUNCTIONS
% -------------------------------------------------------------------------

function updateTickAndTickLabel(data,xTic,yTic,xTLa,yTLa)
% Updates Ticks's and TickLabel's.

% Gets axis limits.
pxLim = get(data.AX,'XLim');
pyLim = get(data.AX,'YLim');
set(data.HA,...
 'XLim'         ,pxLim,...
 'YLim'         ,pyLim)
xLim  = sort(data.xFun(pxLim));
yLim  = sort(data.yFun(pyLim));
if any(~isfinite(xLim)) || ~(diff(xLim)>eps(xLim(1)))
 warning('CVARGAS:addaxes:incorrectXfunctionLim',...
  ['"' func2str(data.xFun) '" is not a monotonical function ' ...
   'or generates NaNs within current X-axis range.'])
 return
end
if any(~isfinite(yLim)) || ~(diff(yLim)>eps(yLim(1)))
 warning('CVARGAS:addaxes:incorrectYfunctionLim',...
  ['"' func2str(data.yFun) '" is not a monotonical function ' ...
   'or generates NaNs within current Y-axis range.'])
 return
end

% Generates a temporal axes with new limits.
tempF = get(0    ,'CurrentFigure');
tempA = get(tempF,'CurrentAxes');
fTemp = figure(...
 'Visible'         ,'off',...
 'Units'           ,get(data.HF,'Units'),...
 'Position'        ,get(data.HF,'Position'));
aTemp = axes(...
 'Parent'          ,fTemp,...
 'Units'           ,get(data.AX,'Units'),...
 'Position'        ,get(data.AX,'Position'),...
 'View'            ,get(data.AX,'View'),...
 'XLim'            ,xLim,...
 'YLim'            ,yLim,...
 data.aOpt{:});
set(0    ,'CurrentFigure',tempF);
set(tempF,'CurrentAxes'  ,tempA);

% Sets date TickLabel's.
if data.xDat
 if strcmpi(data.xDat,'tlabel')
  try
   tlabel(aTemp,'x','keeplimits')
   if isempty(data.xLeg)
    data.xLeg = get(get(aTemp,'XLabel'),'String');
   end
  catch
   % A possible error is to zoom in beyond seconds!
   warning('CVARGAS:addaxes:zoomBeyondSecons',...
    'TLABEL does not work beyond seconds.')
   datetick(aTemp,'x','keeplimits') 
  end
 else
  datetick(aTemp,'x','keeplimits')
 end
end
if data.yDat
 if strcmpi(data.yDat,'tlabel')
  try
   tlabel(aTemp,'y','keeplimits')
   if isempty(data.yLeg)
    data.yLeg = get(get(aTemp,'YLabel'),'String');
   end
  catch
   % A possible error is to zoom in beyond seconds!
   warning('CVARGAS:addaxes:zoomBeyondSecons',...
    'TLABEL does not work beyond seconds.')
   datetick(aTemp,'y','keeplimits')
  end
 else
  datetick(aTemp,'y','keeplimits')
 end
end

% Changes Tick's and TickLabel's.
if (data.xNew~=0) && isempty(xTic)
 % Fixed bug, Sep 2009
 [xTic,xTLa] = changeTickAndTicklabel(...
  get(aTemp,'XTick'),get(aTemp,'XTickLabel'),data.xFun,data.xInv,pxLim,...
                                                                data.xDat);
 if (length(xTic)<2)
  warning('CVARGAS:addaxes:zoomXOut',...
   'Too depth ZOOM on x-axis!')
 else
  dTic = diff(xTic);
  if any(dTic<=eps(min(dTic)))
   warning('CVARGAS:addaxes:incorrectXfunctionTick',...
    ['"' func2str(data.xFun) '" is not a monotonical function ' ...
     'within current X-axis range.'])
   return
  end
 end
end
if (data.yNew~=0) && isempty(yTic)
 % Fixed bug, Sep 2009
 [yTic,yTLa] = changeTickAndTicklabel(...
  get(aTemp,'YTick'),get(aTemp,'YTickLabel'),data.yFun,data.yInv,pyLim,...
                                                                data.yDat);
 if (length(yTic)<2)
  warning('CVARGAS:addaxes:zoomYOut',...
   'Too depth ZOOM on y-axis!')
 else
  dTic = diff(yTic);
  if any(dTic<=eps(min(dTic)))
   warning('CVARGAS:addaxes:incorrectYfunctionTick',...
    ['"' func2str(data.yFun) '" is not a monotonical function ' ...
     'within current Y-axis range.'])
   return
  end
 end
end

% Deletes temporal axes.
delete(aTemp), delete(fTemp)
 
% Updates axes
set(data.HA,...
 'XTick'        ,xTic,...
 'YTick'        ,yTic,...
 'XTickLabel'   ,xTLa,...
 'YTickLabel'   ,yTLa);
 
% Sets legends.
if ~isempty(data.xLeg), set(get(data.HA,'XLabel'),'String',data.xLeg), end
if ~isempty(data.yLeg), set(get(data.HA,'YLabel'),'String',data.yLeg), end

function [Tic,TLa] = changeTickAndTicklabel(Tic,TLa,FUN,INV,LIM,DAT)
% Changes Tick's and TickLabel's.
if ~DAT
 ind = 1; if (str2double(TLa(ind,:))*Tic(ind))==0, ind = 2; end
 if str2double(TLa(ind,:))~=Tic(ind)
  expstr = int2str(log10(Tic(ind)/str2double(TLa(ind,:))));
  if ~strcmp(expstr,'0') && ~strcmp(expstr,'-0')
   ind = double(isspace(TLa));
   ini = arrayfun(@(x)repmat(' ',1,x),sum(cumprod(ind,2),2),...
                                                    'Uniformoutput',false);
   ind = fliplr(ind);
   fin = arrayfun(@(x)repmat(' ',1,x),sum(cumprod(ind,2),2),...
                                                    'Uniformoutput',false);
   TLa = strcat(ini,strtrim(cellstr(TLa)),['e' expstr],fin);
   if iscellstr(TLa), TLa = char(TLa); end
  end
 end
end

for k = 1:length(Tic)
 % Fixed Bug (thanks to Allen Hall), Sep 2009
 if isempty(INV)
  temp = fzero(@(x) FUN(x)-Tic(k),LIM(1));
  if ~((temp>=LIM(1) && (temp<=LIM(2))))
   temp = fzero(@(x) FUN(x)-Tic(k),LIM(2))
   if ~((temp>=LIM(1) && (temp<=LIM(2))))
    temp = fzero(@(x) FUN(x)-Tic(k),mean(LIM))
    datenum(temp)
    if ~((temp>=LIM(1) && (temp<=LIM(2))))
     error('CVARGAS:addaxes:unableToGetInverse',...
      ['Unable to get ''x'' from ' func2str(FUN) ' = ' num2str(Tic(k)) ...
       '. Better try given the inverse formulae.'])
    end
   end
  end
 else
  try
   temp = INV(Tic(k));
  catch
   error('CVARGAS:addaxes:errorWithInverse',...
     ['Unable to get ''x'' from ' func2str(INV) ' = ' num2str(Tic(k)) ...
     '. Error while evaluating.'])
  end
  if ~((temp>=LIM(1) && (temp<=LIM(2))))
     error('CVARGAS:addaxes:incorrectInverse',...
      ['Unable to get ''x'' from ' func2str(INV) ' = ' num2str(Tic(k)) ...
       '. Value out of bounds.'])
  end
 end
 Tic(k) = temp;
end
ind        = ~isfinite(Tic);
Tic(ind)   = [];
TLa(ind,:) = [];
if (length(Tic)>1) && (Tic(1)>Tic(2))
 Tic = sort(Tic);
 TLa = flipud(TLa);
end

function [data,xTic,yTic,xTLa,yTLa,flag] = parseInputs(data,dFun,varargin)
% Parses inputs.

% Defaults.
xTic = [];  xTLa = '';
yTic = [];  yTLa = '';

% Gets axes.
if ~isempty(varargin) && (length(varargin{1})==1) && ishandle(varargin{1})
 data.AX = varargin{1};
 if ~strcmp(get(data.AX,'Type'),'axes')
  error('CVARGAS:addaxes:incorrectHandleInput',...
   'First input must ve a valid axes handle.')
 end
 varargin(1) = [];
 % Updates functions.
 AX        = data.AX;
 data.xFun = eval(func2str(data.xFun));
 data.yFun = eval(func2str(data.yFun));
end

% Checks 'off' option.
flag = false;
if length(varargin)==1
 if isempty(varargin{1})
  % continue
 elseif strcmp(varargin{1},'off')
  flag = true;
  return
 end
 varargin(1) = [];
end

% Checks already used function.
if isappdata(data.AX,data.aNam)
 error('CVARGAS:addaxes:invalifFunctionUse',...
  'By now ADDAXES only can be used once. Use ''off'' option before.')
end

% Checks if paired options.
if rem(length(varargin),2)~=0
 error('CVARGAS:addaxes:incorrectPairedOptions',...
  'Option(s) must be paired property/value(s).')
end

% Saves defaults temporarly in case the given one is empty.
dxFun = data.xFun;
dyFun = data.yFun;

% Loop.
while ~isempty(varargin)
 if isempty(varargin{1})
  error('CVARGAS:addaxes:emptyPropertyName',...
   'Input property name can not be empty.')
 elseif ischar(varargin{1})
  n = length(varargin{1});
  if     strncmpi(varargin{1},'XFunction' ,max(4,n)), data.xFun = varargin{2}; data.xNew = 1;
  elseif strncmpi(varargin{1},'YFunction' ,max(4,n)), data.yFun = varargin{2}; data.yNew = 1;
  elseif strncmpi(varargin{1},'XInverse'  ,max(4,n)), data.xInv = varargin{2};
  elseif strncmpi(varargin{1},'YInverse'  ,max(4,n)), data.yInv = varargin{2};
  elseif strncmpi(varargin{1},'XLegend'   ,max(4,n)), data.xLeg = varargin{2};
  elseif strncmpi(varargin{1},'YLegend'   ,max(4,n)), data.yLeg = varargin{2};
  elseif strncmpi(varargin{1},'XDate'     ,max(4,n)), data.xDat = varargin{2};
  elseif strncmpi(varargin{1},'YDate'     ,max(4,n)), data.yDat = varargin{2};
  elseif strncmpi(varargin{1},'XTick'     ,max(5,n)), xTic      = varargin{2};
  elseif strncmpi(varargin{1},'YTick'     ,max(5,n)), yTic      = varargin{2};
  elseif strncmpi(varargin{1},'XTickLabel',max(6,n)), xTLa      = varargin{2}; 
  elseif strncmpi(varargin{1},'YTickLabel',max(6,n)), yTLa      = varargin{2}; 
  elseif strncmpi(varargin{1},'XColor'            ,max(2,n)) || ...
         strncmpi(varargin{1},'YColor'            ,max(2,n)) || ...
         strncmpi(varargin{1},'CLim'              ,max(2,n)) || ...
         strncmpi(varargin{1},'FontAngle'         ,max(5,n)) || ...
         strncmpi(varargin{1},'FontName'          ,max(5,n)) || ...
         strncmpi(varargin{1},'FontSize'          ,max(5,n)) || ...
         strncmpi(varargin{1},'FontUnits'         ,max(5,n)) || ...
         strncmpi(varargin{1},'FontWeight'        ,max(5,n)) || ...
         strncmpi(varargin{1},'Layer'             ,max(2,n)) || ...
         strncmpi(varargin{1},'LineWidth'         ,max(5,n)) || ...
         strncmpi(varargin{1},'TickDir'           ,max(3,n)) || ...
         strncmpi(varargin{1},'TickLength'        ,max(5,n))
    data.aOpt = {data.aOpt{:},varargin{1},varargin{2}};
  else
    error('CVARGAS:addaxes:invalidProperty',...
     ['Invalid property ''' varargin{1} '''.'])
  end
 else
  error('CVARGAS:addaxes:invalidPropertyType',...
   'Property name must be a valid one.')
 end
 varargin(1:2) = [];
end

% Checks empty functions.
if isempty(data.xFun), data.xFun = dxFun; end
if isempty(data.yFun), data.yFun = dyFun; end

% Checks string functions.
if ischar(data.xFun), data.xFun = str2func(data.xFun); end
if ischar(data.yFun), data.yFun = str2func(data.yFun); end

% Checks ticks and labels.
if ((size(xTLa,1)~=0) && (length(xTic)~=size(xTLa,1))) || ...
   ((size(yTLa,1)~=0) && (length(yTic)~=size(yTLa,1)))
 error('CVARGAS:addaxes:invalidLabels',...
  ['When ''TickLabel''s are specified, '...
   'the respective ''Tick''s should be given too.'])
end

% Checks date axis.
if data.xDat
 if ~ischar(data.xDat)
  data.xDat = dFun;
 else
  data.xDat = lower(data.xDat);
  if ~strcmp(data.xDat,'tlabel') && ~strcmp(data.xDat,'datetick')
   error('CVARGAS:addaxes:invalidDateFunction',...
    ['Unrecognized ''' data.xDat ''' function. ' ...
     'Must be one of ''tlabel'' or ''datetick''.'])
  end
 end
end
if data.yDat
 if ~ischar(data.yDat)
  data.yDat = dFun;
 else
  data.yDat = lower(data.yDat);
  if ~strcmp(data.yDat,'tlabel') && ~strcmp(data.yDat,'datetick')
   error('CVARGAS:addaxes:invalidDateFunction',...
    ['Unrecognized ''' data.yDat ''' function. ' ...
     'Must be one of ''tlabel'' or ''datetick''.'])
  end
 end
end


% [EOF]   addaxes.m