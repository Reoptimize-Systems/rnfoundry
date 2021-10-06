function [Qc, k] = concwindingdesign(p, ndesigns, n, maxi)
% generates conventional winding designs for nonoverlapping concentrated
% windings
% 
% Syntax
%
% [Qc, k] = concwindingdesign(p, ndesigns, n)
%
% Input
%
%  p - the number of Poles in the machine (must be an even number)
%
%  ndesigns - the maximum number of winding designs to return
%
%  n - (optional )the number of coils side-by-side forming a phase group
%    (conventionally 1, with 1 also being the default if not supplied). 
%
%  maxi - (optional) integer controlling the search range for the winding
%    designs, see the description for further information. maxi is set to
%    1000 if not supplied.
%
% Output
%
%  Qc - vector of length ndesigns or less containing the total number of
%    stator coils in each design
%
%  k - vector of k numbers, one for each winding design (see [1])
%
% Description
%
% According to [1], the procedure to determine valid layouts for
% three-phase concentrated-coil windings is as follows: 
%
% 1) Select the number of Poles divisible by two.
%
% 2) Identify those i’s, where i is a
% positive integer that meet the following: 
%
%       36 * (p/6i) ?TRUNC( p / 6i) = k, where k =6, 12, 24, or 30
%
% 3) For n =1, 2,... (i.e., the number of phase coils side by side forming
% a coil phase group), calculate the possible number of stator coils as 
% Qc=3ni.
% 
% concwindingdesign performs this procedure on increasing values of i until
% either ndesigns have been found or maxi has been reached
%
% [1] M. J. Kamper and F. G. Rossouw, “Analysis and Performance of Axial
% Flux Permanent-Magnet Machine With Air-Cored Nonoverlapping Concentrated
% Stator Windings,” IEEE Transactions on Industry Applications, vol. 44,
% no. 5, pp. 1495-1504, Sep. 2008.
%

    if ~iseven(p)
        error('ALLMACHINES:concwindingdesign:oddpoles', ...
              'The number of Poles must be an even number, you supplied: %d', p)
    end
    
    if nargin < 4
        maxi = 1000;
    end
    
    if nargin < 3
        n = 1;
    end
    
    if nargin < 2
        ndesigns = 1;
    end
    
    designcount = 0;
    i = 1;
    Qc = [];
    k = [];
    
    while designcount < ndesigns && i < maxi
        
        % [1] uses the formula 36 * [ p/6i - TRUNC(p/61)]
        [w,num,denom] = rat(fr(p,6*i));

        thisk = 36 * fr(num,denom);

        if thisk == 6 || thisk == 12 || thisk == 24 || thisk == 30
            
            designcount = designcount + 1;
            
            k(designcount) = double(thisk);
            
            Qc(designcount) = 3*n*i;
            
        end
        
        i = i + 1;
    
    end


end