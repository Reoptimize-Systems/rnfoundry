function eul = orient2euler123 (om)
% convert full orientation matrix to euler123 angles
%
% Syntax
%
% eul = orient2euler123 (om)
%
% Description
%
% orient2euler123 converts one or more full 3x3 orientation matrices to
% euler123 angle representations.
%
% Input
%
%  om - either a (3 x 3 x n) matrix of n orientation matrices or a (n x 9)
%   matrix containing the orientiation matrices on each row with the format
%   
%   [ A11, A21, A31, A12, A22, A32, A13, A23, A33 ]
%
%   i.e. the matrix traversed column-wise.
%
% Output
%
%  eul - (n x 3) matrix af euler angles for each input matrix
%
%

    if size (om, 2) == 3
        
       alpha = -atan2 ( om(2,3,:) , om(3,3,:) );
           
       beta = atan2 ( om(1,3,:) ...
                        , ( cos (alpha).*om(3,3,:) - sin (alpha).*om(2,3,:) ) );

       gamma = atan2 ( ( cos (alpha).*om(2,1,:) - sin (alpha).*om(3,1,:) ) ...
                       , ( cos (alpha).*om(2,2,:) - sin (alpha).*om(3,2,:) ) );
                       
                       
%     
%         alpha = atan2 ( om(1,3,:), -om(2,3,:) );
%         
%         beta = atan2 ( -om(2,3,:) .* cos (alpha) + om(1,3,:) .* sin (alpha), om(3,3,:) );
%         
%         gamma = atan2 ( om(3,1,:), -om(3,2,:) );
    
    elseif size (om, 2) == 9
        
       alpha = -atan2 ( om(:,8) , om(:,9) );
           
       beta = atan2 ( om(:,7) ...
                        , ( cos (alpha).*om(:,9) - sin (alpha).*om(:,8) ) );

       gamma = atan2 ( ( cos (alpha).*om(:,2) - sin (alpha).*om(:,3) ) ...
                       , ( cos (alpha).*om(:,5) - sin (alpha).*om(:,6) ) );
        
%         alpha = atan2 ( om(:,7), -om(:,8) );
%         
%         beta = atan2 ( -om(:,8) .* cos (alpha) + om(:,7) .* sin (alpha), om(:,9) );
%         
%         gamma = atan2 ( om(:,3), -om(:,6) );        
        
    end

    eul = [ alpha(:), beta(:), gamma(:) ];
    
end