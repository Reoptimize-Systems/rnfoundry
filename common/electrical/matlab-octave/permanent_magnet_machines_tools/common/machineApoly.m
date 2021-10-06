function Apoly = machineApoly(xycoords, order, A)
% machineApoly: gets the values of the potential at the corrdinates in
% xycoords and returns a 2d polynomial fitted to the values
%
%
    if nargin == 1
        order = 8;
    end

    if nargin < 3
        A = mo_geta2(xycoords);
    end
    
    if ~isoctave
        warning('off', 'MATLAB:nearlySingularMatrix');
    end
    
    Apoly = polyfitn(xycoords, A, order);
    
    if ~isoctave
        warning('on', 'MATLAB:nearlySingularMatrix');
    end

end