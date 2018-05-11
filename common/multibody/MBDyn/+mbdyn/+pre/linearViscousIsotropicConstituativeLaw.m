classdef linearViscousIsotropicConstituativeLaw < mbdyn.pre.constituativeLaw 
    
   properties
       
       viscosity; % the linear viscous coefficient
       
   end
   
   methods
      
       function self = linearViscousIsotropicConstituativeLaw (viscosity)
           % linearViscousIsotropicConstituativeLaw constructor
           %
           %
           
           self.checkNumericScalar (viscosity, true, 'viscosity');
           
           self.type = 'linear viscous isotropic';
           self.viscosity = viscosity;
           
       end
       
        function str = generateMBDynInputString (self)
            % generates MBDyn input string for linearViscousIsotropicConstituativeLaw
            % 
            % Syntax
            %  
            % str = generateMBDynInputString (lvi)
            %  
            % Description
            %  
            % generateMBDynInputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %  
            % Input
            %  
            %  lvi - mbdyn.pre.linearViscousIsotropicConstituativeLaw
            %   object
            %  
            % Output
            %  
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %
            
            str = self.commaSepList (self.type, self.viscosity);     
            
        end
       
   end
    
end