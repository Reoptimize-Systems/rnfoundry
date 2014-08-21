function [kh, kc, ke, beta] = corelosscoeffs (grade, gage, varargin)
% returns Steinmitz core loss coefficients for various grades of lamination
% materials
%
% Syntax
%
% [kh, kc, ke, beta] = corelosscoeffs (grade, gage)
% [kh, kc, ke, beta] = corelosscoeffs (..., 'Parameter', Value, ...)
% 
%
% Possible parameter-value pairs:
%
%  LimitFreqs - 2 element vector containing an upper and lower limit on the
%    frequencies to be used in the fit.
%
%  Processing - 
%
%  InterpolateMissing - 

    Inputs.LimitFreqs = [-inf, inf];
    Inputs.AddZeros = false;
    % following are inputs for corelossdata.m
    Inputs.Processing = 'AsSheared';
    Inputs.InterpolateMissing = true;
    
    Inputs = parse_pv_pairs (Inputs, varargin);

    [fq, Bq, Pq] = corelossdata(grade, gage, ...
                      'Processing', Inputs.Processing, ...
                      'InterpolateMissing', Inputs.InterpolateMissing );

    % strip frequencies of no interest
    Bq = Bq(fq > Inputs.LimitFreqs(1) & fq < Inputs.LimitFreqs(2));
    Pq = Pq(fq > Inputs.LimitFreqs(1) & fq < Inputs.LimitFreqs(2));
    fq = fq(fq > Inputs.LimitFreqs(1) & fq < Inputs.LimitFreqs(2));

    %     kh      kc     ke   beta
    Xo = [0.02, 0.0001, 0.001, 1.5];

    options = LMFnlsq('default');
    options = LMFnlsq(options, 'Display', 0, 'XTol', 1e-9);

    % fit values of kh, kc, ke and beta using ironlossfitfcn
    xf = LMFnlsq( @(fitvars) ironlossfitfcn(fitvars, fq, Bq, Pq), Xo, options );

    xf = abs(xf);

    kh = xf(1);
    kc = xf(2);
    ke = xf(3);
    beta = xf(4);
                
end
