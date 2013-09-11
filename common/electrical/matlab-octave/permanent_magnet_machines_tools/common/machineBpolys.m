function [Bxpoly, Bypoly] = machineBpolys(xycoords, order, B)
% machineApoly: gets the values of the potential at the corrdinates in
% xycoords and returns a 2d polynomial fitted to the values
%
%
    if nargin == 1
        order = 8;
    end

    if nargin < 3
        B = mo_getb2(xycoords);
    end

    if ~isoctave
        warning('off', 'MATLAB:nearlySingularMatrix');
    end

    Bxpoly = polyfitn(xycoords, B(:,1), order);

    Bypoly = polyfitn(xycoords, B(:,2), order);

    if ~isoctave
        warning('on', 'MATLAB:nearlySingularMatrix');
    end

end