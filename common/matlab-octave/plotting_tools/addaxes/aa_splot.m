function varargout = aa_splot(plotdata, varargin)
% AA_SPLOT replaces plot command and automatically changes color
%
%  This function is a replacement for the plot command.  
%  It automatically changes the color with subsequent
%  uses of aa_splot.  
%  
%  The first argument to aa_splot is a structure containing the field yl
%  and yl2. 
%
% Subsequent arguments to aa_splot are exactly the same as the arguments
% for plot.m.

    np = get(gca,'nextplot');
    oldplots = get(gca,'children');

    cord = get(gca,'colorord');

    if ~isempty(oldplots)
        lastcolor = get(oldplots(1),'color');
        if lastcolor == cord(1,:),
            set(gca,'colorord',cord(mod([0:6]+1,7)+1,:));
        end
    end

    hold on;

    h = plot(varargin{:});

    % Store the reconstruction data structure in a field named 'plotdata' in
    % the plot's application data
    
    for i = 1:numel(h)
        setappdata(h(i), 'plotdata', plotdata);
    end

    try % datacursormode is not always available, e.g. in octave

        % Now change figure datacursormode to use our special datatip callback if
        % it has not already been set to do so
        dcm_obj = datacursormode(gcf);

        if strcmpi(get(dcm_obj, 'Enable'), 'on')
        
            % get the handle to the update function of the plot. If this is
            % empty, we may already be in datacursor mode before calling
            % addaxis or addaxisplot. 
            updatefhandle = get(dcm_obj,'UpdateFcn');
            
            if isempty(updatefhandle)
                % we are in datacursormode but no plot update function has been
                % set. This is probably due to the user creating a datatip on a
                % normal plot not created by addaxis or addaxisplot before
                % calling one of these functions. In this case, go ahead and
                % change the UpdateFcn field of our plot.
                set(dcm_obj,'UpdateFcn',@addaxisdtips)
            else
                % there is a function handle, so check if it is the one we
                % expect
                fhandleinfo = functions(updatefhandle);
    
                % warn the user if the function is not addaxisdtips
                if ~strcmpi(fhandleinfo.function, 'addaxisdtips')
                    warning('Figure update function has already been set, but is not set to addaxisdtips. Data tips will not display the correctly scaled data');
                end
            
            end
        else
            % toggle the datacursormode to on and set UpdateFcn to addaxisdtips
            %hdt = datacursormode;
            set(dcm_obj, 'Enable', 'on');
            set(dcm_obj,'UpdateFcn',@addaxisdtips)
        end

    catch
        warning('Unable to set datacursormode, values reported by data cursor will not be correct.')
    end

    for IND = 1:nargout
       varargout(IND) = {h};
    end

    %if nargout > 0, varargout{:} = h; end;

    set(gca,'colorord',cord(mod([0:6]+1,7)+1,:));
    set(gca,'nextplot',np);
    set(gca,'box','on');

end
