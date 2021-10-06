function [fq, Bq, Pq] = m47assheared24gagecorelossdata(dointerp)
% returns table of losses per kg in M-47 As Sheared 24 Gage (0.025") Steel
% laminations at a number of frequencies and field strengths in Tesla
%
% Using data on Specific Core Loss (W/lb) of As-Sheared 24 Gage M47 Fully
% Processed CRNO at Various Frequencies from AK Steel. Data is converted to
% W/kg using a factor of 2.2.
%
% Syntax
%
% [fq, Bq, Pq] = m47assheared26gagecorelossdata()
% [fq, Bq, Pq] = m47assheared26gagecorelossdata(dointerp)
%
% 

    if nargin < 1
        dointerp = true;
    end
    
    m36density = 7700; % kg / m^3
    freqs = [10,20,30,50,60,100,150,200,300,400,600,1000,1500,2000;];
    tesla = [0.1000;0.2000;0.4000;0.7000;1;1.200;1.300;1.400;1.500;1.550;1.600;1.650;1.700;];
    % loss data is provided in W/lb, convert to W/kg using conversion
    % factor of 2.2 lb/kg
    lossperkg = 2.2 * [ 0.001760,0.003790,0.006070,0.01150,0.01460,0.02920,0.05270,0.08130,0.1540,0.2440,0.4720;
                        0.007520,0.01630,0.02620,0.04900,0.06190,0.1220,0.2160,0.3270,0.6050,0.9340,1.720;
                        0.02540,0.05600,0.09100,0.1730,0.2200,0.4360,0.7720,1.160,2.140,3.300,6.020;
                        0.06200,0.1390,0.2290,0.4440,0.5690,1.160,2.100,3.220,6.110,9.580,18.20;
                        0.1100,0.2470,0.4120,0.8150,1.050,2.210,4.120,6.480,12.70,20.60,40.80;
                        0.1500,0.3390,0.5650,1.130,1.460,3.140,5.950,9.520,19.20,31.40,63.70;
                        0.1730,0.3920,0.6550,1.310,1.700,3.670,7.020,11.30,23.10,38.20,77.50;
                        0.2010,0.4550,0.7590,1.520,1.970,4.260,8.200,13.30,27.40,45.90,92.60;
                        0.2350,0.5310,0.8860,1.770,2.280,4.930,9.510,15.40,31.50,52.80,109;
                        0.2540,0.5740,0.9560,1.910,2.470,5.320,10.20,16.60,33.60,56.40,118;
                        0.2730,0.6160,1.030,2.040,2.650,5.710,11,17.80,35.50,NaN,NaN;
                        0.2910,0.6570,1.100,2.180,2.830,6.110,11.70,19,NaN,NaN,NaN;
                        0.3070,0.6950,1.160,2.310,3,6.480,12.40,20.20,NaN,NaN,NaN;];
                        
    losspervol = lossperkg .* m36density; % (P / kg) * (kg / m^3) = P / m^3
    
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

