function [design, simoptions] = simfun_AM(design, simoptions)
% performs common preprocessing tasks common to all machines prior to a
% dynamic simulation
%
% Syntax
%
% [design, simoptions] = simfun_AM(design, simoptions)
%
% Input
%
% design - a structure intended to hold the machine design parameters.
%   simfun_AM looks for the following fields in the sturcture, and if not
%   present adds them with the specified default.
%
%   Field                  Default
%  
%   FEMMTol                1e-5
%   WireResistivityBase    1.7e-8
%   AlphaResistivity       3.93e-3
%   TemperatureBase        20
%   CoilLayers             1
%   NStrands               1
%   NStages                1
%   sides                  design.NStages
%   HcMag                  979000
%   MagnetSkew             0
%
%   The properties of the machine coils are also checked, see
%   checkcoilprops_AM.m for details of what must be present for this, and
%   what will be added or modified in the structure.
%
%   
%
% See also: checkcoilprops_AM.m
%

% Copyright Richard Crozier 2012-2013

    % store a tolerance for use in FEMM simulations in the design structure
    design = setfieldifabsent(design, 'FEMMTol', 1e-5);
    
    % use resistivity of copper at 20 degrees C
    design = setfieldifabsent(design, 'WireResistivityBase', 1.7e-8);
    
    design = setfieldifabsent(design, 'AlphaResistivity', 3.93e-3);
    
    design = setfieldifabsent(design, 'TemperatureBase', 20);
    
    % set the magnet coercivity if not present
    design = setfieldifabsent( design, 'HcMag', 979000);
    
    % set the magnet skew if not present
    design = setfieldifabsent( design, 'MagnetSkew', 0);
    
    % set numeber of magnets per pole, to create the magnet skew if not
    % already present
    design = setfieldifabsent( design, 'NSkewMagnetsPerPole', 10);
    
    % default temperature (used for resistance calculations is 20 degrees)
    simoptions = setfieldifabsent(simoptions, 'Temperature', 20);
    
    % some machines can optionally skip the FEA if it has already been
    % done. We here set the default option to not skip any FEA step
    simoptions = setfieldifabsent(simoptions, 'SkipFEA', false);
    
    % further to this, a more fine choice can sometimes be desired to skip
    % only the main field FEA step, or to skip only an inductance
    % simulation, the following options by default set these options to
    % false
    simoptions = setfieldifabsent(simoptions, 'SkipMainFEA', false);
    
    simoptions = setfieldifabsent(simoptions, 'SkipInductanceFEA', false);

    simoptions = setfieldifabsent(simoptions, 'GetVariableGapForce', false);
    
    simoptions = setfieldifabsent(simoptions, 'usefemm', false);
	simoptions = setfieldifabsent(simoptions, 'quietfemm', true);
    
    simoptions = setfieldifabsent(simoptions, 'SkipCheckCoilProps', false);
    
    if ~simoptions.SkipCheckCoilProps
        % Determine the area of the coil, and the find either the number of
        % turns, coil wire diameter if required
        design = checkcoilprops_AM(design);
    end
    
    % If loss data is supplied, fit the loss function coefficients to it,
    % this loss data is expected to be in the form of three vectors
    % containing a set of frequencies, max B values and the resulting
    % losses at each combination
    if isfield(design, 'CoreLoss')
        % check the required data fields are present before attempting the
        % fit
        if all(isfield(design.CoreLoss, {'fq', 'Bq', 'Pq'})) ...
             && ~all(isfield (design.CoreLoss, {'kh', 'kc', 'ke', 'beta'}))
            
            for clind = 1:numel(design.CoreLoss)
                
                %     kh      kc     ke   beta
                Xo = [0.02, 0.0001, 0.001, 1.5];

                options = LMFnlsq('default');
                options = LMFnlsq(options, 'Display', 0, 'XTol', 1e-9);

                % fit values of kh, kc, ke and beta using ironlossfitfcn
                xf = LMFnlsq(@(fitvars) ironlossfitfcn(fitvars,design.CoreLoss(clind).fq,design.CoreLoss(clind).Bq,design.CoreLoss(clind).Pq), Xo, options);

                xf = abs(xf);
                
                design.CoreLoss(clind).kh = xf(1);
                design.CoreLoss(clind).kc = xf(2);
                design.CoreLoss(clind).ke = xf(3);
                design.CoreLoss(clind).beta = xf(4);
            
            end
            
        end
        
    end
    
end