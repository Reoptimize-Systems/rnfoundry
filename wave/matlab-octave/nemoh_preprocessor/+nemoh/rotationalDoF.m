classdef rotationalDoF < nemoh.degreeOfFreedom
    
   properties 

       rotationCentre;
       
   end
   
   
   methods
       
       function  self = rotationalDoF (direction, rot_centre, varargin)
           
           options.Comment = 'Rotational degree of freedom (first three values are direction vector, other three are centre of rotation)';
           
           options = parse_pv_pairs (options, varargin);
       
           self = self@nemoh.degreeOfFreedom(direction, 'Comment', options.Comment);
           
           self.rotationCentre = rot_centre;
           
       end
       
       function str = generateDoFStr (self)

           str = sprintf ( '2 %s %s %s %s %s %s\t\t! %s\n', ...
                           self.formatNumber (self.directionVector(1)), ...
                           self.formatNumber (self.directionVector(2)), ...
                           self.formatNumber (self.directionVector(3)), ...
                           self.formatNumber (self.rotationCentre(1)), ...
                           self.formatNumber (self.rotationCentre(2)), ...
                           self.formatNumber (self.rotationCentre(3)), ...
                           self.comment );           
       end
       
       function set.rotationCentre (self, rot_centre)
           
           assert ( isnumeric (rot_centre) ...
                      && numel (rot_centre) == 3 ...
                      && isreal (rot_centre) ...
                    , 'rotationCentre must be a real numeric vector of length 3' );
                
            self.rotationCentre = rot_centre;
            
       end
       
       
   end
    
end