function ObjVal = objcorelesstorus(Chrom, rtn_type, simoptions, multicoredir)
% evaluates coreless torus score
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

   FieldBounds = [ 0.1,   0.9;    %  1. taucmiVtaucmo 
                   0.1,   6.0;    %  2. tcVtm
                   0.0,   0.5;    %  3. gVtm
                   0.1,   1.0;    %  4. tbiiVtbio
                   0.1,   2.0;    %  5. tbioVtm
                   1.005, 4.0;    %  6. tmVtaumm
                   0.5,   0.95;   %  7. taummVtaupm
                   0.5,   0.99;   %  8. RsVRbi 
                   0.8,   0.99;   %  9. RbiVRmi
                   0.1,   0.9;    % 10. RmiVRmo
                   0.5,   3.0;    % 11. Rmo
                   0,     0.3;    % 12. NBasicWindings
                   0.5,   15.0;   % 13. RgVRc
                   0.2,   0.85;   % 14. kfill
                   0,     1;      % 15. DcAreaFac
                   5,     50;     % 16. pole pairs
                   0,     1;      % 17. branchfac                     
                   0.5,   1;      % 18. modulefac
                   1,     10 ];   % 19. NStages
                
   Dim = size(FieldBounds, 1);
   
% Compute population parameters
   [Nind,Nvar] = size(Chrom);

% Check size of Chrom and do the appropriate thing
   % if Chrom is [], then define size of boundary-matrix and values
   if Nind == 0
      % return text of title for graphic output
      if rtn_type == 2
         ObjVal = ['TORUS CORELESS-' int2str(Dim)];
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
           
           multicoredir = fullfile(rootpath, 'Temp');
           
       end
       
       ObjVal = objelectricalmachine(simoptions, Chrom, ...
                                     'chrom2design_TORUS_CORELESS' , ...
                                     'evaluatedesign_TORUS_CORELESS', ...
                                     multicoredir);
       
   else
      error('size of matrix Chrom is not correct for function evaluation');
   end   

end