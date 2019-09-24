function [b, hax, hfig] = stackedForcePlot (logger, bodynum, axis, varargin)

    options.HydroForceVarNames = { 'ForceExcitationRamp', ...
                                   'ForceRadiationDamping', ...
                                   'ForceMorrison', ...
                                   'ForceViscousDamping', ...
                                   'ForceAddedMassUncorrected', ...
                                   'ForceRestoring' };
                                
	options.OtherForces = [];
    options.Times = [ logger.data.Time(1), logger.getPrevLoggedVal('Time', 0)];
    options.Skip = 1;
    
    options = parse_pv_pairs (options, varargin);
    
    time_inds = find (logger.data.Time >= options.Times(1) & logger.data.Time <= options.Times(2));
    
    time_inds = time_inds(1:options.Skip:end);
    
    n_stacked_vars = numel(options.HydroForceVarNames);
    var_names = options.HydroForceVarNames;
    
    data = nan ( numel(time_inds), n_stacked_vars );
    
    for ind = 1:n_stacked_vars
        data(:,ind) = squeeze (logger.data.(var_names{ind})(axis,bodynum,time_inds));
    end
    
    hfig = figure;
    hax = axes;
    
    if n_stacked_vars > 1
        x_offset = mean (diff (logger.data.Time(time_inds))) / (2 * n_stacked_vars);
    else
        x_offset = 0;
    end
    this_x_offset = -x_offset * n_stacked_vars;
    
    
    b = [];
    hold (hax, 'all');
    
%     cord = get (hax,'colororder');
    
    for ind = 1:n_stacked_vars
        
        if ind > 1
            b = [ b, errorbar( hax, ...
                               logger.data.Time(time_inds) + this_x_offset, ...
                               data(:,ind) - data(:,ind)./2 + sum(data(:,1:ind-1), 2), ...
                               abs(data(:,ind))./2 ) ...
                ];
            
        else
            b = errorbar(hax, logger.data.Time(time_inds) + this_x_offset, data(:,ind) - data(:,ind)./2, abs(data(:,ind))./2);
        end
        
        set(b(ind), 'LineStyle', 'none');
        set(b(ind), 'LineWidth', 2);
%         set(b(ind), 'Color', cord(1,:));
        set(b(ind), 'CapSize', 0);
        
%         cord = circshift (cord, [1,0]);
        
        this_x_offset = this_x_offset + x_offset;
    end
    
    plot (hax, logger.data.Time(time_inds), squeeze (logger.data.ForceHydro(axis,bodynum,time_inds)));
    
    hold (hax, 'off');
    
    legend ( var_names{:}, 'ForceHydro' );
    
    
%                 PTO_1_InternalForce: [1035×1 double]
%          PTO_1_RelativeDisplacement: [1035×1 double]
%              PTO_1_RelativeVelocity: [1035×1 double]
%      PTO_1_ControllerTargetVelocity: [1035×1 double]
%               PTO_1_ForcePTODamping: [1035×1 double]
%                PTO_1_ForcePTOSpring: [1035×1 double]
%         PTO_1_ForceControllerTcAmax: [1035×1 double]
%     PTO_1_ForceRateControllerTcAmax: [1035×1 double] };
    
    

end