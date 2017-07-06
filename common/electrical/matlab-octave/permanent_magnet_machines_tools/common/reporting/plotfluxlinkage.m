function plotfluxlinkage (design, varargin)
% plots the flux linkage for a machine design with relative displacement
%
% Syntax
%
%  design - machine design structure. The flux linkage data must have been
%    generated for this machine
%
% Additional options are avaialable as Parameter-Value pairs
%
%  'ShowEMF' - flag determining whether to show the derivative of the flux
%    linkage curve to give an idea of the shape of the EMF waveform,
%    default is true
%
%  'ShowKnots' - flag determining whether to show the 'knot' positions in
%    the pieceswise curve fit. The flux linkage is stored in a piecewise 
%
%  'ShowRawData' - flag determining whether to plot the raw data from which
%    the slm curve fit was produced, default is true.
%
%

    options.ShowEMF = true;
    options.ShowKnots = true;
    options.ShowRawData = true;
    
    options = parse_pv_pairs (options, varargin);
    
    if ~isfield (design, 'slm_fluxlinkage')
        error ('This design does not appear to have the flux linkage information present, you may need to run its associated simulation funcitons')
    end
    
    if options.ShowEMF
        plotstyle = 'dy';
    else
        plotstyle = '';
    end
    
    plotslm( design.slm_fluxlinkage, plotstyle, ...
         'ShowKnots', options.ShowKnots, ...
         'ShowData', options.ShowRawData, ...
         'Title', 'Flux Linkage With Position', ...
         'XLabel', 'Relative Position [Poles]', ...
         'YLabel', 'Flux Linkage [WbT]')

end