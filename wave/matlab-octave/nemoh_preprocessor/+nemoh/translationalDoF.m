classdef translationalDoF < nemoh.degreeOfFreedom
    
   properties 

   end
   
   
   methods
       
       function  self = translationalDoF (direction, varargin)
           
           options.Comment = 'Translational degree of freedom (first three values are direction vector, other three can be ignored)';
           
           options = parse_pv_pairs (options, varargin);
       
           self = self@nemoh.degreeOfFreedom(direction, 'Comment', options.Comment);
           
       end
       
       function str = generateDoFStr (self)

           str = sprintf ( '1 %s %s %s 0. 0. 0.\t\t! %s\n', ...
                           self.formatNumber (self.directionVector(1)), ...
                           self.formatNumber (self.directionVector(2)), ...
                           self.formatNumber (self.directionVector(3)), ...
                           self.comment );
           
       end
       
       
   end
    
end