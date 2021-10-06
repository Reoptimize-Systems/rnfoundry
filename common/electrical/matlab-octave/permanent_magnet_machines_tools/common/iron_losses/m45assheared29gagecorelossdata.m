function [fq, Bq, Pq] = m45assheared29gagecorelossdata(dointerp)
% returns table of losses per kg in M-45 AS Sheared 29 Gage (0.014") Steel
% laminations at a number of frequencies and field strengths in Tesla
%
% Syntax
%
% [fq, Bq, Pq] = m45assheared26gagecorelossdata()
% [fq, Bq, Pq] = m45assheared26gagecorelossdata(dointerp)
%
% 

    if nargin < 1
        dointerp = true;
    end
    
    matdensity = 7700; % kg / m^3
    freqs = [10,20,30,50,60,100,150,200,300,400,600,1000,1500,2000;];
    tesla = [0.1000;0.2000;0.4000;0.7000;1;1.200;1.300;1.400;1.500;1.550;1.600;1.650;1.700;];
    % loss data is provided in W/lb, convert to W/kg using conversion
    % factor of 2.2 lb/kg
    lossperkg = 2.2 * [ 0.001520,0.003200,0.005010,0.008930,0.01106,0.02070,0.03460,0.05080,0.08940,0.1370,0.2510;
                        0.006320,0.01330,0.02100,0.03750,0.04660,0.08650,0.1450,0.2110,0.3670,0.5520,1;
                        0.02130,0.04530,0.07140,0.1300,0.1613,0.3050,0.5170,0.7600,1.330,2.010,3.640;
                        0.05250,0.1120,0.1780,0.3250,0.4070,0.7810,1.350,2,3.580,5.480,10.10;0.09370,0.2000,0.3180,0.5840,0.7330,1.420,2.480,3.740,6.810,10.70,20;
                        0.1290,0.2750,0.4370,0.8040,1.009,1.960,3.450,5.230,9.650,15.20,29.30;
                        0.1500,0.3210,0.5100,0.9380,1.177,2.290,4.010,6.100,11.30,17.90,34.60;0.1760,0.3760,0.5980,1.099,1.377,2.670,4.690,7.120,13.20,20.90,40.70;
                        0.2090,0.4450,0.7070,1.298,1.626,3.150,5.500,8.360,15.50,24.40,NaN;
                        0.2270,0.4840,0.7680,1.408,1.763,3.420,5.970,9.040,16.70,29.40,NaN;
                        0.2460,0.5240,0.8300,1.523,1.904,3.690,6.450,9.730,17.90,28.40,NaN;
                        0.2920,0.5590,0.8860,1.625,2.030,3.930,6.860,10.40,NaN,NaN,NaN;
                        0.2780,0.5890,0.9330,1.708,2.130,4.130,7.240,10.90,NaN,NaN,NaN; ];
                       
    losspervol = lossperkg .* matdensity; % (P / kg) * (kg / m^3) = P / m^3
    
    if dointerp

        lossperkgnonans = infillpower(freqs, losspervol);

        % add a row of zeros along the top
        losspervolnonans = [zeros(1, size(losspervolnonans, 2)); losspervolnonans];
        losspervolnonans = [zeros(size(losspervolnonans, 1), 1), losspervolnonans];

        % ensure row vectors
        freqs = freqs(:);
        tesla = tesla(:);

        % add zero values for freq and tesla
        freqs = [0; freqs];
        tesla = [0; tesla];

        % make a grid of interpolation points
        fi = linspace(min(freqs), max(freqs), 15); 
        Bi = linspace(min(tesla), max(tesla), 15);
        [fq,Bq] = meshgrid(fi,Bi);

        losspervollist = reshape(losspervolnonans, [], 1);

        Blist = repmat(tesla, numel(freqs), 1);

        freqlist = reshape(repmat(freqs(:)', [], numel(tesla)), [], 1);

        Pq = griddata(freqlist,Blist,losspervollist,fq,Bq);

        Pq(1,:) = 0;
        Pq(:,1) = 0;
    
    else
        
        [fq,Bq,Pq] = table2vectors(freqs,tesla,losspervol,true);
        
    end

end

