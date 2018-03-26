classdef scalarFunctionElasticIsotropicConstituativeLaw < mbdyn.pre.constituativeLaw 
    
   properties
       
       scalarFunction;
       
   end
   
   methods
      
       function self = scalarFunctionElasticIsotropicConstituativeLaw (scalar_function)
           
           
           assert (isa (scalar_function, 'mbdyn.pre.scalarFunction'), ...
               'scalar_function must be an mbdyn.pre.scalarFunction object' );
           
           self.type = 'scalar function elastic isotropic';
           self.scalarFunction = scalar_function;
           
       end
       
        function str = generateMBDynInputString (self)
                       
            str = self.commaSepList (self.type, self.scalarFunction.generateMBDynInputString ());     
            
        end
       
   end
    
end