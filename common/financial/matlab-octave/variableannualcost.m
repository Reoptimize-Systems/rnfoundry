function vac = variableannualcost(varargin)
% calculates the amortised annual cost of capital at a given discount rate
% over a project lifetime
%
% Syntax
%
% vac = variableannualcost(fixedcost, drate, years)
%
% Description
%
% The formula used to calculate the annual cost is that given in 'A Brief
% Review of Wave Energy', T W Thorpe, A report produced for The UK
% Department of Trade and Industry:
%
%     vac = fixedcost * drate / ( 1 - (1 + drate)^(-years)); 
%
% however, if the discount rate is zero the following simpler formula is
% used instead:
%
%     vac = fixedcost / years
%
% Inputs can be any combination of scalars and matrices of the same size,
% scalar values will be expanded to the same size as the other matrices.
%

% Copyright Richard Crozer, The University of Edinburgh

    error(nargchk(3,3,nargin,'struct'))
    
    % check for nans and inf values (although inf values are allowed for
    % the fixed cost and number of years
    if any(isnan(varargin{1})) || any(isnan(varargin{2})) || any(isnan(varargin{3}))
        error('Input contained NaN values')
    elseif any(isinf(varargin{2})) || any(isinf(varargin{3}))
        error('Either drate or years contained Inf values')
    elseif any(varargin{3} < 1)
        error('All values of years must be greater than or equal to one.')
    end
    
    % check for a combination of scalars and same size matrices
    [tf, scalarinds, nonscalarinds] = samesizeorscalars(varargin{:});
    
    if ~tf
        error('All inputs must be scalar or matrices of the same size')
    end
    
    % expand the scalars to the same size as the other matrices
    if any(scalarinds) && ~isempty(nonscalarinds)
        
        varargin(scalarinds) = cellfun(@(x) repmat(x, size(varargin{nonscalarinds(1)})), varargin(scalarinds), 'UniformOutput', false);
                        
    end
    
    fixedcost = varargin{1};
    drate = varargin{2};
    years = varargin{3};
    
    % formula used in 'A Brief Review of Wave Energy' A report produced for
    % The UK Department of Trade and Industry, T W Thorpe
    vac = fixedcost .* drate ./ ( 1 - (1 + drate).^(-years)); 
   
    % formula seen elsewhere (gives exactly the same result as Thorpe)
%     vac = ( fixedcost .* drate .* (1 + drate).^years ) ./ ...
%           ((1 + drate).^years - 1);
    
    % replace zero dicount rate values with the fixed cost divided by the
    % years (otherwise there is zero/zero which is NaN)
    vac(drate==0) = fixedcost(drate==0) ./ years(drate==0);
 
end