function A = UnequalLegAngle (vars)
% Calculates the cross-sectional area of a beam with an unequal legged
% angle cross-section, as calculated in 'Roark's Formulas Stress & Strain'
%
%
%            t     y
%         ->   <-  .
%            _     .
%        ^  | |    .
%        :  | |    .
%        :  | |    .
%     d  :  |.|.................. x
%        :  | |    .
%        :  | |    .         :
%        :  | |____________  v
%        v  |______._______|   t
%           <------.------>  ^
%               b  .         :
%                  .
%
%
% Input: 
%
%   vars - (n x 3) matrix of values as below
%           vars(:,1) - d, vertical leg length
%           vars(:,2) - b, horizontal leg length
%           vars(:,3) - t, leg thickness
%
% Output:
%
%   A - (n x 1) column vector of values of A, the cross-sectional area of a
%       beam with a unequal-legged angle profile
%

    d = vars(:,1);
    b = vars(:,2);
    t = vars(:,3);
    
    A = t .* (b + d - t);

end