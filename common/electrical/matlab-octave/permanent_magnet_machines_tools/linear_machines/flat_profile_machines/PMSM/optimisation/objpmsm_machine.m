function ObjVal = objpmsm_machine(Chrom, rtn_type, simoptions, multicoredir)
% objpmsm_machine: evaluates the score for a chromosome of variables
% describing PMSM machine designs
%
% Syntax:  ObjVal = objactm(Chrom, rtn_type, simoptions, evaloptions)
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

%   WmVWp - magnet width to pole width ratio
%
%   WtVWc - tooth width to combined tooth and slot width ratio
%
%   hmVWm - magnet height to magnet width ratio
% 
%   htVWt - tooth height to tooth width ratio
% 
%   hbaVht - armature back iron height to tooth height ratio
% 
%   hbfVhm - field back iron height to magnet height ratio
% 
%   lsVWp - stack length to pole width ratio
% 
%   gVhm - air gap to magnet height ratio
% 
%   Wp - pole width
    
    maxbeams = size(beamvars(simoptions.Evaluation.GuideRailIMethod), 1);

    FieldBounds = [ 0.6,   0.95;       %  1. WmVWp 
                    0.1,   0.8;        %  2. WtVWc
                    0.05,  3.0;        %  3. hmVWm
                    0.1,   10.0;       %  4. htVWt
                    0.05,  3.0;        %  5. hbaVht
                    0.1,   2.0;        %  6. hbfVhm
                    1.0,   20.0;       %  7. lsVWp
                    0,     3.0;        %  8. gVhm
                    0.025, 0.4;        %  9. Wp
                    0.5,   15.0;       % 10. RlVRp
                    0.2,   0.65;       % 11. kfill
                    0,     1;          % 12. DcAreaFac
                    5,     200;        % 13. TransPoles 
                    0,     1;          % 14. BranchFac                      
                    0,     2;          % 15. BeamSpreadFactor
                    1,     10;         % 16. NoOfMachines
                    1,     maxbeams;   % 17. Beam to use for field guide rails
                    0.1,   3;          % 18. Outer Webs per m on field
                    2,     20 ];       % 19. maxAllowedxT Factor 
                
    Dim = size(FieldBounds, 1);
   
    % Compute population parameters
    [Nind,Nvar] = size(Chrom);

    % Check size of Chrom and do the appropriate thing
    % if Chrom is [], then define size of boundary-matrix and values
    if Nind == 0
        % return text of title for graphic output
        if rtn_type == 2
            ObjVal = ['PMSM MACHINE-' int2str(Dim)];
            % return value of global minimum
        elseif rtn_type == 3
            ObjVal = 0;
            % define size of boundary-matrix and values
        else
            % return the bounds of the variables
            ObjVal = FieldBounds';
        end

    elseif Nvar == Dim
        
        if nargin < 4
            
            multicoredir = fullfile(pmsmrootpath(), 'Temp');
            
        end
       
        ObjVal = objelectricalmachine(simoptions, Chrom, ...
                                      'preprocsystemdesign_PMSM' , ...
                                      'designandevaluate_PMSM', ...
                                      multicoredir, ...
                                      false);
    else
        error('size of matrix Chrom is not correct for function evaluation');
    end

end