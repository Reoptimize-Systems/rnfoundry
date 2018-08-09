classdef linearElasticIsotropicConstituativeLaw < mbdyn.pre.constituativeLaw 
    
   properties
       
       stiffness;
       
   end
   
   methods
      
       function self = linearElasticIsotropicConstituativeLaw (stiffness)
           
           
           self.checkNumericScalar (stiffness, true, 'stiffness');
           
           self.type = 'linear elastic isotropic';
           self.stiffness = stiffness;
           
       end
       
        function str = generateMBDynInputString (self)
                       
            str = self.commaSepList (self.type, self.stiffness);     
            
        end
       
   end
    
end