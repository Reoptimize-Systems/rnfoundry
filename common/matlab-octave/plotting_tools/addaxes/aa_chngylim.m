function aa_chngylim(schemaProp, evnt)

    hax = evnt.AffectedObject;
    
    origylim = get(hax, 'YLim');
    
    newylim = evnt.newValue;
    
    normnewylim1 = newylim ./ origylim;
    
    normylimshift = normnewylim1 - ([1,1]);
    
    aad = getappdata(hax, 'addaxis_data');
    
    for i = 1:numel(aad.axisdata)
        
        thisaxylim = get(aad.axisdata{i}(1), 'YLim');
        
        thisaxnewylim = thisaxylim + (thisaxylim .* normylimshift);
        
        set(aad.axisdata{i}(1), 'YLim', thisaxnewylim);
        
        % now change the scaling of all the plots associated with the sub
        % axis
        for ii = 2:numel(aad.axisdata{i})
            % The information necessary to reconstruct the actual data from
            % the position is stored in a field named 'plotdata' in the
            % plot's application data. This is set by aa_splot.
            plotdata = getappdata(aad.axisdata{i}(ii), 'plotdata');
            
            y = get(aad.axisdata{i}(ii), 'YData');
            
            plotdata.yl = newylim;
            plotdata.yl2 = thisaxnewylim;
            
            y = (y-plotdata.yl(1)) ./ (plotdata.yl(2)-plotdata.yl(1)) ...
                .* (plotdata.yl2(2)-plotdata.yl2(1)) + plotdata.yl2(1);
            
            y = (y-plotdata.yl2(1)) ./ (plotdata.yl2(2)-plotdata.yl2(1)) .* (plotdata.yl(2)-plotdata.yl(1)) + plotdata.yl(1);
            
            set(aad.axisdata{i}(ii), 'YData', y);
            
            setappdata(aad.axisdata{i}(ii), 'plotdata', plotdata);
            
        end
        
    end

end