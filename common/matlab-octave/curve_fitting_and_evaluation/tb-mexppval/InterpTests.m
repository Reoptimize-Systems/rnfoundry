%function [Tics1D1,Tics1D,Tics2D,Tics3D,Tics4D] = InterpTests
function  Y = InterpTests(nSites,nLoops)

if nargin == 0
    nSites = 1;
end
if nargin <= 1
    nLoops = 1;
end


% univariate test case
X = 0:10;
Y = sqrt(X);
spline = fn2fm(spapi({aptknt(X,3)},{X},Y),'pp');

Tics1Dm = zeros(1,nLoops);

for j = 1:nLoops
      
    X = rand(1,nSites);
    
    tic;yy = ppmval(X,spline); y1 = toc;
    
    % Built in Matlab version
    tic;yyy = ppual(spline,X); y2 = toc;
    
    norm(yy-yyy)
    
    Tics1Dm(j) = y2./y1;
 
end
        
% Bivariate test
X = 0:10;
Y = 0:10;
[XX,YY] = meshgrid(X,Y);
W = sin(sqrt(XX.^2+YY.^2));
spline = fn2fm(spapi({aptknt(X,3),aptknt(Y,3)},{X,Y},W),'pp');

Tics2D = zeros(1,nLoops);

for j = 1:nLoops
   
   X = rand(2,nSites);
    
    tic;yy = ppmval(X,spline); y1 = toc;
    
    % Built in Matlab version
    tic;yyy = ppual(spline,X); y2 = toc;
    
    norm(yy-yyy)
    
    Tics2D(j) = y2./y1;
    
end

        % Trivariate test case
X = 0:10;
Y = 0:10;
Z = 0:10;
[XX,YY,ZZ] = meshgrid(X,Y,Z);
W = sin(sqrt(XX.^2+YY.^2+ZZ.^2));
spline = fn2fm(spapi({aptknt(X,3),aptknt(Y,3),aptknt(Z,3)},{X,Y,Z},W),'pp');

Tics3D = zeros(1,nLoops);

for j = 1:nLoops
   
   X = rand(3,nSites);
    
    tic;yy = ppmval(X,spline); y1 = toc;
    
    % Built in Matlab version
    tic;yyy = ppual(spline,X); y2 = toc;
    
    norm(yy-yyy)
    
    Tics3D(j) = y2./y1;
    
end

                            % 4-D test case
X = 0:10;
Y = 0:10;
Z = 0:10;
W = 0:10;
[XX,YY,ZZ,WW] = ndgrid(X,Y,Z,W);
Q = sin(sqrt(XX.^2+YY.^2+ZZ.^2+WW.^2));
spline = fn2fm(spapi({aptknt(X,3),aptknt(Y,3),aptknt(Z,3),aptknt(W,4)},{X,Y,Z,W},Q),'pp');


Tics4D =  zeros(1,nLoops);
for j = 1:nLoops
   
   X = rand(4,nSites);
    
    tic;yy = ppmval(X,spline); y1 = toc;
    
    % Built in Matlab version
    tic;yyy = ppual(spline,X); y2 = toc;
    
    norm(yy-yyy)
    
    Tics4D(j) = y2./y1;
    
end

            % Finally test the univariate case
X = 0:10;
Y = sqrt(X);
spline = interp1(X,Y,'cubic','pp');

Tics1D =  zeros(1,nLoops);
for j = 1:nLoops
   
   X = rand(1,nSites);
    
    tic;yy = ppuval(X,spline); y1 = toc;
    
    % Built in Matlab version
    tic;yyy = ppual(spline,X); y2 = toc;
    
    norm(yy-yyy)
    
    Tics1D(j) = y2./y1;
    
end


Y = [Tics1D(:),Tics1Dm(:),Tics2D(:),Tics3D(:),Tics4D(:)];









