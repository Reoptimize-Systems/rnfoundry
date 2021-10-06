function Force = gapclosingforce_ACPMSM(dgVlm, bpVlm, taupVbp, lsVbp, dbiVlm, lm, pClosing)
% gapclosingforce_ACPMSM: a function to find the force between the two
% magnets in the linear air-cored permanent magnet synchronous machine
%
% Input: 
%
%   dgVlm - distance from back of magnet to machine centre versus the
%           magnet thickness ratio
%
%   bpVlm - magnet width to magnet depth ratio
%
%   taupVbp - pole pitch to magnet width ratio
%
%   lsVbp - machine depth to magnet width ratio
%
%   lm - The magnet width in the same direction as the force
%
%   pClosing - (optional) A structure containing the polynomial to be
%           used to evaluate the air-gap closing force in the machine
%
% Output:
%
%   Force - The total per-pole force between each side of the field.
%

    % First calculate the Force / Area from the polynomial
    if nargin > 6
        Force = polyvaln(pClosing, [dgVlm, bpVlm, taupVbp, lsVbp, dbiVlm, lm]);
    else
        ACPMSMdir = fileparts(which('ratios2dimensions_ACPMSM'));
        load(fullfile(ACPMSMdir, 'Polynomials', 'ClosingForcePoly_ACPMSM.mat'));
        Force = polyvaln(pClosing, [dgVlm, bpVlm, taupVbp, lsVbp, dbiVlm, lm]);
    end
    
    % Now calculate the actual total force
    bp = bpVlm * lm;
    
    Force = Force * (taupVbp * bp) * (lsVbp * bp);

end