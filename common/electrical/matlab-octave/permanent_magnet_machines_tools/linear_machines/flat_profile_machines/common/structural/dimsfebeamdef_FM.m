function [x, y, z, n, znodes, zguidecon, zsupp, zframe, zextreme, guideyshift, tolerance] = dimsfebeamdef_FM(BeamInfo, tolerance)
% determine the locations of importent positions on the flat machine frame
% ensuring nodes within a certain distance of each other will be merged.
%
% Syntax
%
% [x, y, z, n, znodes, zguidecon, zsupp, zframe, ...
%       zextreme, guideyshift, tolerance] = dimsfebeamdef_FM(BeamInfo)
%
% [x, y, z, n, znodes, zguidecon, zsupp, zframe, ...
%       zextreme, guideyshift, tolerance] = dimsfebeamdef_FM(BeamInfo, tolerance)
%
% Input
%
%   BeamInfo - A structure containing information on the frame design
%
%   tolerance - optional scalar value to use as the tolerance for merging
%     nearby node locations
%
% Output
%
%   x
% 
%   y 
% 
%   z 
% 
%   n 
% 
%   znodes 
% 
%   zguidecon
% 
%   zsupp 
% 
%   zframe
% 
%   zextreme
% 
%   guideyshift
% 
%   tolerance - tolerance used in determining if nodes should be merged or
%     kept separate
%

    % Dimensions and parameters
    if nargin < 2
        tolerance = 100 * eps;
    end

    % The centre of the machine is located a (0,0,0) so it is convenient to
    % half the dimensions to find the appropriate coordinates for the model
    x = BeamInfo.width / 2;
    y = BeamInfo.depth / 2;
    z = BeamInfo.height / 2;

    % calculate appropriate points for the members keeping the two field parts
    % apart
    if BeamInfo.OuterWebs.NoPerSide == 1
        % if there's only one outer web, it has to go in the middle
        zsupp = 0;
    else
        % otherwize arrange the supports evenly spaced from the extremities
        % inward
        zsupp = (-z:2*z/(BeamInfo.OuterWebs.NoPerSide-1):z)';
    end

    % each support beam will be divided into n elements
    n = BeamInfo.OuterPoleSupports.NoPerSide - 1;
    
    % calculate appropriate node heights for the beams on the outer back iron
    zframe = (-z:2*z/n:z)';

    % there will also be connections to the linear guide rails
    zguidecon = (-z+(2*z/(BeamInfo.GuideBearings.NoPerGuide+1)):2*z/(BeamInfo.GuideBearings.NoPerGuide+1):z-(2*z/(BeamInfo.GuideBearings.NoPerGuide+1)))';

    % get the unique members of all the nodes, put them into a column vector
    % and sort it in ascending order
    %     znodes = sortrows(consolidator([zsupp; zframe; zguidecon],[],[],tolerance));
    %     znodes = sortrows([zsupp; zframe; zguidecon]);
    %     znodes = znodes(abs(znodes(2:end) - znodes(1:end-1)) > tolerance);
    znodes = sortrows(uniquetol([zsupp; zframe; zguidecon],tolerance));

    % we must shift the support beams in the positive and negative y directions
    % to accomodate the connections to the linear guides
    guideyshift = 50/1000;
    
    % zextreme = (design.Poles(1) .* design.Wp);
    zextreme = BeamInfo.GuideRails.length / 2;
    
end