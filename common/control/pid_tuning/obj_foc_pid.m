function ObjVal = obj_foc_pid(Chrom, rtn_type, FieldBounds, design, simoptions)
% objactm: evaluates actm score as part of a heaving buoy wec system
%
% Syntax:  ObjVal = objactm_machine(Chrom, rtn_type, simoptions, multicoredir)
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

% Copyright Richard Crozer, The University of Edinburgh
                
   Dim = size(FieldBounds, 1);
   
% Compute population parameters
   [Nind,Nvar] = size(Chrom);

   % Check size of Chrom and do the appropriate thing
   % if Chrom is [], then define size of boundary-matrix and values
   if Nind == 0
       % return text of title for graphic output
       if rtn_type == 2
           ObjVal = ['FOC PID Optimisation-' int2str(Dim)];
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
       
       parfor ind = 1:size (Chrom,1)
           ObjVal(ind,1) = dosimandscoreit (Chrom(ind,:), design, simoptions);
       end
       
   else
       error('size of matrix Chrom is not correct for function evaluation');
   end

end

function score = dosimandscoreit (chrom, design, simoptions)

    Kp_d = chrom(1);
    Ki_d = chrom(2);
    Kd_d = 0;
    Kp_q = chrom(3);
    Ki_q = chrom(4);
    Kd_q = 0;

    design.FOControl.PI_d = pidController ( Kp_d, Ki_d, Kd_d, 'MaxOut', 1000, 'MinOut', -1000, 'InitialTime', simoptions.FOCTApp );
    design.FOControl.PI_q = pidController ( Kp_q, Ki_q, Kd_q, 'MaxOut', 1000, 'MinOut', -1000, 'InitialTime', simoptions.FOCTApp );
    
    [T, Y, results, design, simoptions] = simulatemachine_AM ( design, ...
                                                               simoptions );
                                                           
	TStart = simoptions.FOCTApp;
    TStartind = find (T>=TStart, 1);
    ITAE = trapz ( T(TStartind+1:end), T(TStartind+1:end) .* abs (results.PI_vals(TStartind+1:end,1:2)), 1 );
    
    score = sum (ITAE);
    
end
