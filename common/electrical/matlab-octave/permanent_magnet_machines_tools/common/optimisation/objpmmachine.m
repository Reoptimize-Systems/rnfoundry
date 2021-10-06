function ObjVal = objpmmachine(Chrom, rtn_type, simoptions, multicoredir)
% objective function for permanent magnet machines suitible for use with
% the GA toolbox from the university of Sheffield UK
%
% Syntax
% 
% ObjVal = objpmmachine(Chrom, rtn_type, simoptions, multicoredir)
%
% Input parameters:
%
%    Chrom - Matrix containing the chromosomes of the current
%      population. Each row corresponds to one individual's string
%      representation. if Chrom == [], then special values will be returned
%
%    rtn_type - if Chrom == [], and rtn_type == 1 (or []), the boundaries
%      of the variables to be optimised will be returned. If rtn_type == 2
%      returns a string, used for the title in optimisation progress
%      graphs. If rtn_type == 3 the value of global minimum is returned. 
%
%    simoptions - structure that must contain at least the fields:
%
%      FieldBounds: (n x 2) matrix of lower and upper bounds for the
%        variables subject to the optimisation.
%
%      Chrom2DesignFcn: String or function handle for the function to be
%        used to convert the chromosomal representation of the design into
%        a full design. Note that this function is never evaluated in
%        parallel, always in a loop for the full chromosome locally.
%        Therefore you should reserve any intensive computation for the
%        DesignEvaluationFcn (see below). Chrom2DesignFcn must have the
%        calling syntax:
%
%      DesignEvaluationFcn: String or function handle for the function to
%        be used to evaluate the designs created by the Chrom2DesignFcn.
%
%    multicoredir - (optional) The evaluation can be performed in parallel
%      using the Multicore package. This argument sets the directory for
%      Multicore to use for temporary files.
%
% Output parameters:
%
%    ObjVal    - Column vector containing the objective values of the
%                individuals in the current population.
%                if called with Chrom == [], then ObjVal contains
%                rtn_type == 1, matrix with the boundaries of the function
%                rtn_type == 2, text for the title of the graphic output
%                rtn_type == 3, value of global minimum
%                
%
% Author:     Richard Crozier

    Dim = size(simoptions.FieldBounds, 1);

    % Compute population parameters
    [Nind,Nvar] = size(Chrom);

    % Check size of Chrom and do the appropriate thing
    % if Chrom is [], then define size of boundary-matrix and values
    if Nind == 0
        % return text of title for graphic output
        if rtn_type == 2
            ObjVal = ['Optimisation Progress-' int2str(Dim)];
            % return value of global minimum
        elseif rtn_type == 3
            ObjVal = 0;
            % define size of boundary-matrix and values
        else
            % return the bounds of the variables
            ObjVal = simoptions.FieldBounds';
        end

    elseif Nvar == Dim

        if nargin < 4
            % if a multicore temporary files directory is not supplied, use
            % one in the optimisation directory
            multicoredir = fullfile ( pm_machines_tools_rootdir (), ...
                                     'common', ...
                                     'optimisation', ...
                                     'temp' );
            
            mkdir (multicoredir);

        end

        % run the optimisation on the current population
        ObjVal = objelectricalmachine ( simoptions, ...
                                        Chrom, ...
                                        simoptions.Chrom2DesignFcn, ...
                                        simoptions.DesignEvaluationFcn, ...
                                        multicoredir, ...
                                        false );

    else
        error('size of matrix Chrom is not correct for function evaluation');
    end

end