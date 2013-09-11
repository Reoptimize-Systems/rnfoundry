function psidot = dpsidx_linear(design, x)

    if isfield(design, 'psipoly') 
        % First check if a polynomial giving psi is available
        [psidot, err] = derivest(@(var) psipolypsi_linear(design, var), x, 'NominalStep', 0.99, 'FixedStep', 0.5);
        
        % psidot is corrected here for the unitless polynomial fitting dimension
        psidot = psidot ./ design.Taup;
        
    elseif isfield(design, 'psilookup') 
        
        % If no polynomial is available check if there is a lookup table of
        % values
        [psidot, err] = derivest(@(var) psilookuppsi_linear(design, var), x,  'NominalStep', 0.99);
        
    elseif isfield(design, 'Apoly')
        
        % if no psi polynomial but there is a polynomial giving the vector
        % potential, use this
        [psidot, err] = derivest(@(var) psi_linear(design, var), x,  'NominalStep', 0.99);
     
    else
        error(['You must supply either a flux linkage polynomial, or ',...
            'lookup table, or a vector potential polynomial for ', ...
            'the design'])
    end

end