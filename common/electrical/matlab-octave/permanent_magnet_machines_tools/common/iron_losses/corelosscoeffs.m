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
%
%  Plot - display a plot of the data. Set to 'log' for a log-log plot, and
%    'linear'
%
%  FitMode - 

    Inputs.LimitFreqs = [-inf, inf];
    Inputs.LimitInduction = [-inf, inf];
    Inputs.AddZeros = false;
    % following are inputs for corelossdata.m
    Inputs.Processing = 'AsSheared';
    Inputs.InterpolateMissing = true;
    Inputs.Plot = 'none';
    Inputs.FitMode = 'simple';
    
    Inputs = parse_pv_pairs (Inputs, varargin);

    [fq, Bq, Pq] = corelossdata(grade, gage, ...
                      'Processing', Inputs.Processing, ...
                      'InterpolateMissing', Inputs.InterpolateMissing );

    % strip frequencies of no interest
    Bq = Bq(fq > Inputs.LimitFreqs(1) & fq < Inputs.LimitFreqs(2));
    Pq = Pq(fq > Inputs.LimitFreqs(1) & fq < Inputs.LimitFreqs(2));
    fq = fq(fq > Inputs.LimitFreqs(1) & fq < Inputs.LimitFreqs(2));
    
    % strip inductions of no interest
    Pq = Pq(Bq > Inputs.LimitInduction(1) & Bq < Inputs.LimitInduction(2));
    fq = fq(Bq > Inputs.LimitInduction(1) & Bq < Inputs.LimitInduction(2));
    Bq = Bq(Bq > Inputs.LimitInduction(1) & Bq < Inputs.LimitInduction(2));

    if strcmpi ( Inputs.FitMode, 'simple')
        
        %     kh      kc     ke   beta
        Xo = [0.02, 0.0001, 0.001, 1.5];

        options = LMFnlsq('default');
        options = LMFnlsq(options, 'Display', 5000, 'XTol', 1e-12, 'FunTol', 1e-12, 'MaxIter', 200000);

        % fit values of kh, kc, ke and beta using ironlossfitfcn
        xf = LMFnlsq( @(fitvars) ironlossfitfcn(fitvars, fq, Bq, Pq), Xo, options );

        xf = abs(xf);

        kh = xf(1);
        kc = xf(2);
        ke = xf(3);
        beta = xf(4);
        if beta > 10, beta = 10; end
    
    else
        
        Xo = [ 0.2, 0.0001, 0.001, 0,0,1.5 ];
       
        options = LMFnlsq('default');
        options = LMFnlsq(options, 'Display', 5000, 'XTol', eps, 'FunTol', eps, 'MaxIter', 200000);

        % fit values of kh, kc, ke and beta using ironlossfitfcn
        xf = LMFnlsq( @(fitvars) ironlossfitfcn_complex(fitvars, fq, Bq, Pq), Xo, options );

        xf = abs(xf);

        kh = xf(1);
        kc = xf(2);
        ke = xf(3);
        beta = xf(4:end);
       
    end
    
    if strcmpi(Inputs.Plot, '3DLog')
        
        scatter3 (fq, Bq, Pq, 'kx');
        hold on
        scatter3 (fq, Bq, steinmitz(fq, Bq, kh, kc, ke, beta));
        hold off
        
        set (gca, 'XScale', 'linear');
        set (gca, 'YScale', 'linear');
        set (gca, 'ZScale', 'log');
        xlabel ('Frequency [Hz]');
        ylabel ('Induction [T]');
        zlabel ('Power Loss [Wm^{-2}]');
        legend ('Data', 'Fitted');
        set(gcf, 'Color', 'w');
        
    elseif strcmpi(Inputs.Plot, '3DLinear')
        
        scatter3 (fq, Bq, Pq, 'kx');
        hold on
        scatter3 (fq, Bq, steinmitz(fq, Bq, kh, kc, ke, beta));
        hold off
        
        set (gca, 'XScale', 'linear');
        set (gca, 'YScale', 'linear');
        xlabel ('Frequency [Hz]');
        ylabel ('Induction [T]');
        zlabel ('Power Loss [Wm^{-2}]');
        legend ('Data', 'Fitted');
        set(gcf, 'Color', 'w');
        
    elseif strcmpi(Inputs.Plot, '2DLog')
        
        cols = numel(unique (fq));
        
        fq = reshape (fq, [], cols);
        Bq = reshape (Bq, [], cols);
        Pq = reshape (Pq, [], cols);
        
        plot (Bq, Pq, ':');
        hold on
        plot (Bq, steinmitz(fq, Bq, kh, kc, ke, beta), 'x');
        hold off
        
        set (gca, 'XScale', 'linear');
        set (gca, 'YScale', 'log');
        xlabel ('Induction [T]');
        ylabel ('Power Loss [Wm^{-2}]');
        
        legstrs = {};
        for ind = 1:size(fq, 2)
            legstrs = [legstrs, {sprintf('Data %d Hz', fq(1,ind))}];
        end
        for ind = 1:size(fq, 2)
            legstrs = [legstrs, {sprintf('Fit %d Hz', fq(1,ind))}];
        end
        legend (legstrs{:}, 'Location', 'EastOutside');
        set(gcf, 'Color', 'w');
    end
                
end

function res = ironlossfitfcn_complex (fitvars, f, B, coreloss)

    fitvars = [abs(fitvars(1:3)); fitvars(4:end)];

    kh = fitvars(1);
    kc = fitvars(2);
    ke = fitvars(3);
    beta = fitvars(4:end);
    
    % res = FUN(x) - y
    fitcoreloss = steinmitz (f, B, kh, kc, ke, beta);
    
    res = real(fitcoreloss - coreloss(:));
    
end

function val = steinmitz (f, B, kh, kc, ke, beta)
    
    if numel(beta) ~= 1

        beta = beta(1).*f(:) + beta(1).*B(:) + beta(3);

    end
    
    if beta > 10, beta = 10; end
    
    val = abs(kh) .* f .* B.^abs(beta) + abs(kc) .* (f.*B).^2 + abs(ke) .* (f.*B).^1.5;

end


