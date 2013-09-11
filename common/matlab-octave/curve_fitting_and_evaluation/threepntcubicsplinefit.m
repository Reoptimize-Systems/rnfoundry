function [ p1p2, p1p2der1, p1p2der2 ] = threepntcubicsplinefit(pnts)
% fits two cubic spline polynomials to one or more sets of three
% coordinates
%
% Syntax
%
% [ p1p2, p1p2der ] = threepntcubicsplinefit(pnts)
%
% Input
%
%   pnts - (3 x 2 x n) matrix of n sets of three points. Two smoothly
%     joined spline polynomials will be fitted to each of the three point
%     sets. The first column is the x coordinate of the points, the second
%     the y coordinate.
%
% Output
%
%   p1p2 - (2 x 4 x n) matrix, the top slice (1,:,:) contains the
%     coefficients of the first spline between points 1 and 2 in each of
%     the n coordinate sets. The bottom slice (2,:,:) contains the
%     coefficients of the second polynomial between points 2 and three 
%
%   p1p2der1 - (2 x 3 x n), same as the p1p2, but the coefficients of the
%     first derivative of the splines in this case.
%
%   p1p2der2 - (2 x 3 x n), same as the p1p2der1, but the coefficients of 
%     the second derivative of the splines in this case.
%
% Example 1
%
% % simple example
% pnts = [ 1, 1; 2, 5; 3, 4 ];
% p1p2 = threepntcubicsplinefit(pnts);
% figure; plot(pnts(:,1), pnts(:,2), 'ok');
% hold on;
% x = linspace(pnts(1,1), pnts(2,1), 100);
% y = p1p2(1,1)*x.^3 + p1p2(1,2)*x.^2 + p1p2(1,3)*x + p1p2(1,4);
% plot(x, y, 'r');
% x = linspace(pnts(2,1), pnts(3,1), 100);
% y = p1p2(2,1)*x.^3 + p1p2(2,2)*x.^2 + p1p2(2,3)*x + p1p2(2,4);
% plot(x, y, 'b');
% hold off
% 
% Example 2
%
% % multiple point sets, 2 in this case
% pnts = cat(3, [ 1, 1; 2, 5; 3, 4 ], [ 0, 0; 2, 3; 3, 1 ]);
% p1p2 = threepntcubicsplinefit(pnts);
% figure; plot(pnts(:,1,1), pnts(:,2,1), 'ob');
% hold on; 
% x = linspace(pnts(1,1,1), pnts(2,1,1), 100);
% y = p1p2(1,1,1)*x.^3 + p1p2(1,2,1)*x.^2 + p1p2(1,3,1)*x + p1p2(1,4,1);
% plot(x, y, 'b');
% x = linspace(pnts(2,1,1), pnts(3,1,1), 100);
% y = p1p2(2,1,1)*x.^3 + p1p2(2,2,1)*x.^2 + p1p2(2,3,1)*x + p1p2(2,4,1);
% plot(x, y, 'b');
% plot(pnts(:,1,2), pnts(:,2,2), 'or');
% x = linspace(pnts(1,1,2), pnts(2,1,2), 100);
% y = p1p2(1,1,2)*x.^3 + p1p2(1,2,2)*x.^2 + p1p2(1,3,2)*x + p1p2(1,4,2);
% plot(x, y, 'r');
% x = linspace(pnts(2,1,2), pnts(3,1,2), 100);
% y = p1p2(2,1,2)*x.^3 + p1p2(2,2,2)*x.^2 + p1p2(2,3,2)*x + p1p2(2,4,2);
% plot(x, y, 'r');
% hold off
% 
% Example 3
%
% % the first example, but with derivatives too
% pnts = [ 1, 1; 2, 5; 3, 4 ];
% [ p1p2, p1p2der1, p1p2der2 ] = threepntcubicsplinefit(pnts);
% figure; plot(pnts(:,1), pnts(:,2), 'ok');
% hold on; 
% x = linspace(pnts(1,1), pnts(2,1), 100);
% y = p1p2(1,1)*x.^3 + p1p2(1,2)*x.^2 + p1p2(1,3)*x + p1p2(1,4);
% plot(x, y, 'b');
% x = linspace(pnts(2,1), pnts(3,1), 100);
% y = p1p2(2,1)*x.^3 + p1p2(2,2)*x.^2 + p1p2(2,3)*x + p1p2(2,4);
% plot(x, y, 'b');
% x = linspace(pnts(1,1), pnts(2,1), 100);
% y = p1p2der1(1,1)*x.^2 + p1p2der1(1,2)*x + p1p2der1(1,3);
% plot(x, y, 'r');
% x = linspace(pnts(2,1), pnts(3,1), 100);
% y = p1p2der1(2,1)*x.^2 + p1p2der1(2,2)*x + p1p2der1(2,3);
% plot(x, y, 'r');
% x = linspace(pnts(1,1), pnts(2,1), 100);
% y = p1p2der2(1,1)*x + p1p2der2(1,2);
% plot(x, y, 'm');
% x = linspace(pnts(2,1), pnts(3,1), 100);
% y = p1p2der2(2,1)*x + p1p2der2(2,2);
% plot(x, y, 'm');
% hold off

% Created by Richard Crozier 2013

    p = size(pnts,3);
    zels = zeros(1,1,p);
    oneels = zels + 1;
    twoels = oneels + 1;
    
    A = ...
    [ pnts(1,1,:).^3,   pnts(1,1,:).^2, pnts(1,1,:),  oneels,  zels,              zels,            zels,         zels; ...
      pnts(2,1,:).^3,   pnts(2,1,:).^2, pnts(2,1,:),  oneels,  zels,              zels,            zels,         zels; ...
      zels,             zels,           zels,         zels,    pnts(2,1,:).^3,    pnts(2,1,:).^2,  pnts(2,1,:),  oneels; ...
      zels,             zels,           zels,         zels,    pnts(3,1,:).^3,    pnts(3,1,:).^2,  pnts(3,1,:),  oneels; ...
      3*pnts(2,1,:).^2, 2*pnts(2,1,:),  oneels,       zels,    -3*pnts(2,1,:).^2, -2*pnts(2,1,:),  -oneels,      zels; ...
      6*pnts(2,1,:),    twoels,         zels,         zels,    -6*pnts(2,1,:) ,   -twoels,         zels,         zels; ...
      6*pnts(1,1,:),    twoels,         zels,         zels,    zels,              zels,            zels,         zels; ...
      zels,             zels,           zels,         zels,    6*pnts(3,1,:) ,    twoels,          zels,         zels ];

    b = [ pnts(1,2,:);
          pnts(2,2,:);
          pnts(2,2,:);
          pnts(3,2,:)
          zels
          zels
          zels
          zels ];
      
    
    % preallocate space for the solutions
    x = zeros(8,1,p);
    
    % solve the system, wish bsxfun worked for this
    for ind = 1:p
        x(1:8,1,ind) = A(:,:,ind) \ b(:,:,ind);
    end
    
    % transpose along dimension 3
    x = permute(x, [2 1 3]);
    
    % spline coefficients
    p1p2 = [ x(1,1:4,:); x(1,5:8,:) ];
    
    % spline first derivative coefficients
    if nargout > 1
        p1p2der1 = [ 3 .* p1p2(:,1,:), 2 .* p1p2(:,2,:), p1p2(:,3,:) ] ;
    end
    
    % spline second derivative coefficients
    if nargout > 2
        p1p2der2 = [ 2 .* p1p2der1(:,1,:), p1p2der1(:,2,:) ] ;
    end

end