function length = isotrapzcoilmtl(w1, w2, h, t)
% estimates the mean turn length in an isosceles trapezoidal coil winding 
%
% Syntax
% 
% length = isotrapzcoilwirelen(w1, w2, h, t)
% 
% Input
%                                      w1
%                        :<-------------------------->:
%                        :                            :
%                   ***************************************                 
%                 *                                         *               
%                 *                                          *              
%   ............. *       *****************************      *              
%           ^      *     *                             *     *              
%           �       *    *                             *    *               
%           �        *    *                           *    *                
%           �         *    *                         *    *                 
%           �          *    *                       *    *                  
%        h  �           *    *               ----> *    *<---- t   
%           �            *    *                   *    *                    
%           �             *    *                 *    *                     
%           �              *    *               *    *                      
%           �               *    *             *    *                       
%           �                *    *           *    *                        
%           v                 *    *         *    *                         
%    ........................  *     *******     *                          
%                               *    :     :    *                           
%                                *   :     :   *                                                     
%                                  **:     :**                              
%                                    *******                                
%                                    :     :                                 
%                                    :<--->:
%                                      w2
%                                                       
%   w1 - the internal length of the top parallel edge
%
%   w1 - the internal length of the top parallel edge
%
%   h - the distance from the internal bottom edge to the internal top edge
%
%   t - the coil winding thickness, i.e. the cross-sectional length of the
%       coil
% 
% Output
%
%   length - the mean turn length of the wire in the coil
%

    % calculate the length of the non-parallel sides
    sidelen = sqrt( (h+t).^2 - (((w1+t) - (w2+t)) ./ 2).^2 );
    
    % now calculate the turn length
    length = (2 .* sidelen) + w1 + w2 + (2*t) + (pi .* t);
                                         
end