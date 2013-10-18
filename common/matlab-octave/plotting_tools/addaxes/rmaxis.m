function rmaxis(n, axishandle)
% Removes an axis added using addaxis (and all plots associated with this
% this axis)
%
% Syntax
%
% rmaxis(n)
%
% 
%

    if nargin < 2
        axishandle = gca;
    end
    
    aad = getappdata(axishandle, 'addaxis_data');
    
    if ischar(n)
        % some day we wil have named axes we can search for
    end
    
    for plotlineind = 2:numel(aad.axisdata{n-1})
        delete(aad.axisdata{n-1}(plotlineind));
    end
    
    delete(aad.axisdata{n-1}(1));
    
    % remove the axis from the list of axes
    aad.axisdata(n-1) = [];
    
    setappdata(axishandle,'addaxis_data',aad);
  
end
  