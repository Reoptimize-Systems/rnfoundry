function [score, design, simoptions, T, Y, results] = regenfrommpgaoptfile_AM (optfile)
% regenerates a design and simulation results from an mpga output file
%
% Syntax
%
% [score, design, simoptions, T, Y, results] = regenfrommpgaoptfile_AM (optfile)
%
% Input
%
%  optfile - optimisation progress file created by the 'mpga' function.
%    This is the same file used to resume progress of an optimisation from a
%    previous run.
%
%    regenfrommpgaoptfile_AM reconstructs and evaluates the design using
%    the 'mpga' function, which re-runs the objective function, but returns
%    all results returned by this function.
%
% Output
%
%  score - the score for the design
%
%  design - the reconstructed design structure
%
%  simoptions - the simulation options used in the evaluation fucntion
%
%  T - ode solution time points
%
%  Y - ode solution outputs
%
%  results - other results and outputs from the ode solution
%

    mpgaoptions.resumefile = optfile;
    mpgaoptions.evaluate = true;

    out = mpga (mpgaoptions);

    score = out{1}{1};
    design = out{1}{2};
    simoptions = out{1}{3};
    T = out{1}{4};
    Y = out{1}{5};
    results = out{1}{6};

end