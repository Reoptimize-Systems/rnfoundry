function [kh, kc, ke, beta] = corelosscoeffs (grade, gage, varargin)

    [fq, Bq, Pq] = corelossdata(grade, gage, varargin{:});

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
