function ObjVal = objactiam_prescribed_mot(Chrom, rtn_type, simoptions, multicoredir)
% objactm: evaluates actiam score
%
% Syntax:  ObjVal = objactm(Chrom, rtn_type, simoptions, options)
%
% Input parameters:
%    Chrom     - Matrix containing the chromosomes of the current
%                population. Each row corresponds to one individual's
%                string representation.
%                if Chrom == [], then special values will be returned
%    rtn_type  - if Chrom == [] and
%                rtn_type == 1 (or []) return boundaries
%                rtn_type == 2 return title
%                rtn_type == 3 return value of global minimum
%
% Output parameters:
%    ObjVal    - Column vector containing the objective values of the
%                individuals in the current population.
%                if called with Chrom == [], then ObjVal contains
%                rtn_type == 1, matrix with the boundaries of the function
%                rtn_type == 2, text for the title of the graphic output
%                rtn_type == 3, value of global minimum
%                
%
% Author:     Richard Crozier
% History:    11.05.2010
%

   FieldBounds = [ 0.05,  0.95;   %  1. WmVWp
                   0.05,  2.0;    %  2. WpVRm   
                   1.005, 3.0;    %  3. RoVRi  
                   1.05,  1.5;    %  4. RaVRo
                   0.01,  0.975;  %  5. RsoVRm  
                   1.001, 1.1;    %  6. RiVRm 
                   0.1,   1/3;    %  7. WcVWp
                   0.01,  1.0;    %  8. Rm
                   0.5,   15.0;   %  9. RgVRc
                   0.2,   0.65;   % 10. kfill
                   0,     1;      % 11. DcAreaFac
                   10,    200;    % 12. TransPoles 
                   0,     1;      % 13. BranchFac                     
                   0,     5;   ]; % 14. bpoints
                
   Dim = size(FieldBounds, 1);
   
% Compute population parameters
   [Nind,Nvar] = size(Chrom);

% Check size of Chrom and do the appropriate thing
   % if Chrom is [], then define size of boundary-matrix and values
   if Nind == 0
      % return text of title for graphic output
      if rtn_type == 2
         ObjVal = ['ACTIAM MACHINE-' int2str(Dim)];
      % return value of global minimum
      elseif rtn_type == 3
         ObjVal = 0;
      % define size of boundary-matrix and values
      else
          
         %numbuoys = buoynum2buoydata(simoptions.buoylibdir);
         
         % return the bounds of the variables
         ObjVal = FieldBounds';
      end

   elseif Nvar == Dim
       
       if nargin < 4
           
           rootpath = fileparts(which('GetEMF_ACTIAM'));

           multicoredir = fullfile(rootpath, 'Temp');
           
       end
       
       ObjVal = objelectricalmachine(simoptions, Chrom, ...
                                     'chrom2design_prescribed_mot_ACTIAM' , ...
                                     'designandevaluate_ACTIAM', ...
                                     multicoredir);
       
   else
      error('size of matrix Chrom is not correct for function evaluation');
   end   

end