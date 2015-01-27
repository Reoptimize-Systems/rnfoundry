function docs

    % both ppmval and ppual take two arguments: Evalsites X, and polynomial description pp
    % that must be in pp-form. Error is printed, if polynomial is not in pp-form.
    
    % ppual is tuned to handle functions, that are mapping from R->R. Input argument 
    % can be of any shape and return value is of same dimension as X
    
    % ppmval is designed to handle all other cases. That is functions can be mappings from 
    % R^m -> R^n.
    % If function is a map from R^m -> R^n, then X must have m rows, in
    % which case each column of X is interpreted to be independent input variable
    % (x1,...,xm). Multiple evaluation sites run in columns.
    % If X has dimensions m X k, then return matrix is of size n X k
    
    % Both functions perform extrapolation in a such way, that they extend
    % polynomial out of range values. 
    % If function is out of range from left
    % in dimension i, then first piece of polynomial in dimension i is used
    % to extrapolate in that dimension
    
    % On the other hand, if function is out of range in dimension i, then
    % last piece of the polynomial in dimension i is used to extrapolate
    
    % Examples
    
    % Univariate case
    %X = 10:100;
    %Y = sqrt(X);
    %spline = interp1(X,Y,'cubic','pp');
    
    %Xi = 0:110;
    %Yi = ppuval(Xi,spline);
    
    % Case in higher dimension
    %X = 0:10;
    %Y = 0:10;
    %Z = 0:10;
    %[XX,YY,ZZ] = meshgrid(X,Y,Z);
    %W = sin(sqrt(XX.^2+YY.^2+ZZ.^2));
    %spline = fn2fm(spapi({aptknt(X,3),aptknt(Y,3),aptknt(Z,3)},{X,Y,Z},W),'pp');
    
    %Xi = [1.2,5;6,2;3,4];
    %Yi = ppmval(Xi,spline);
    
end