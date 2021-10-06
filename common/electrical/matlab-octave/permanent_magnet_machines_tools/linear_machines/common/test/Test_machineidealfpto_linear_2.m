% Test_machineidealfpto_linear_2

% load actm_design_and_simoptions.mat

load generator_design_and_simoptions.mat

%%

Fpto = 100000;
xBh = 0;
xBs = 0;
vBh = 1;
vBs = 0;

tstep = 0.001;

ii = 1;

time = 0;

Icoils(1,1:3) = [0,0,0];
[EMF(ii, 1:3), Icoils(ii,1:3), limitIcoils(ii,1:3), Force(ii, 1)] = machineidealfpto_linear(design, simoptions, Fpto, xBh, xBs, vBh, vBs, Icoils(ii,1:3)');
xBh = xBh + (vBh*tstep);
ii = ii + 1;
time(ii) = time(ii-1) + tstep;

ii = 2;
while xBh < 0.5
    
    
    [EMF(ii, 1:3), Icoils(ii,1:3), limitIcoils(ii,1:3), Force(ii, 1)] = machineidealfpto_linear(design, simoptions, Fpto, xBh, xBs, vBh, vBs, Icoils(ii-1,1:3)');
    
    xBh = xBh + (vBh*tstep);
    
    ii = ii + 1;
    
    time(ii) = time(ii-1) + tstep;
    
end

% plot the three-phase EMFs
plot(time(1:end-1), EMF(:,1), 'r', time(1:end-1), EMF(:,2), 'y', time(1:end-1), EMF(:,3), 'b');
ylabel('EMFs');

% Add the currents to the graph
addaxis(time(1:end-1), Icoils(:,1), 'xr');
addaxislabel(2,'Coil Currents');

addaxisplot(time(1:end-1), Icoils(:,2), 2, 'xy');
addaxisplot(time(1:end-1), Icoils(:,3), 2, 'xb');

% Add the force to the graph
addaxis(time(1:end-1), Force, 'm');
addaxislabel(3,'Force');

%% Test with variable pto force

clear Force Icoils EMF

xBh = design.PoleWidth ./ 3;
xBs = 0;
vBh = 1;
vBs = 0;

ii = 1;

Fpto = -200000:1000:200000;

Icoils(1,1:3) = [0,0,0];
[EMF(ii, 1:3), Icoils(ii,1:3), limitIcoils(ii,1:3), Force(ii, 1)] = machineidealfpto_linear(design, simoptions, Fpto(ii), xBh, xBs, vBh, vBs, Icoils(1,1:3)');
    
% ii = ii + 1;
    
for ii = 2:numel(Fpto)
    
    [EMF(ii, 1:3), Icoils(ii,1:3), limitIcoils(ii,1:3), Force(ii, 1)] = machineidealfpto_linear(design, simoptions, Fpto(ii), xBh, xBs, vBh, vBs, Icoils(ii-1,1:3)');
    
%     ii = ii + 1;
    
end

figure; plot(Fpto, Force); 

%% Test limited currents version

clear Force Icoils EMF

xBh = design.PoleWidth ./ 3;
xBs = 0;
vBh = 1;
vBs = 0;

ii = 1;

Fpto = -200000:100:200000;

Icoils(1,1:3) = [0,0,0];
[EMF(ii, 1:3), Icoils(ii,1:3), limitIcoils(ii,1:3), Force(ii, 1)] = machineidealfpto_linear(design, simoptions, Fpto(ii), xBh, xBs, vBh, vBs, Icoils(1,1:3)');
    
% ii = ii + 1;
    
for ii = 2:numel(Fpto)
    
    [EMF(ii, 1:3), Icoils(ii,1:3), limitIcoils(ii,1:3), Force(ii, 1)] = machineidealfpto_linear(design, simoptions, Fpto(ii), xBh, xBs, vBh, vBs, limitIcoils(ii-1,1:3)');
    
%     ii = ii + 1;
    
end

figure; plot(Fpto, Force); 
    


%% Demonstrate how force limit changes with position

Fpto = 100000;
vBh = 0;
xBh = 0:design.Wp/50:3*design.Wp;
xBs = 0;
vBs = 0;

for ii = 1:numel(xBh)
    
    [EMF(ii, 1:3), Icoils(ii,1:3), limitIcoils(ii,1:3), Force(ii, 1)] = machineidealfpto_linear(design, simoptions, Fpto, xBh(ii), xBs, vBh, vBs);
    
end

titlefsz = 16;
labelfsz = 14;

subplot(3,1,1) 
plot(xBh./ design.Wp, repmat(Fpto, size(Force)) ./ 1e3, ':r', xBh ./ design.Wp, Force ./ 1e3, 'b');
%title('Force Requested and Supplied', 'FontSize', titlefsz);
% xlabel('Position Relative to Pole Width (Unitless)', 'FontSize', labelfsz)
ylabel('Force (kN)', 'FontSize', labelfsz)
legend('Requested Force', 'Supplied Force')
legend('boxoff')
ylimits = get(gca, 'YLim');
set(gca, 'YLim', [ylimits(1), ylimits(2) .* 1.1], 'FontSize', 12)

subplot(3,1,2)
plot(xBh./design.Wp, Icoils(:,1), ':r', xBh./design.Wp, Icoils(:,2), 'g', xBh./design.Wp, Icoils(:,3), '--b');
% title('Current Required to Supply Force', 'FontSize', titlefsz);
% xlabel('Position Relative to Pole Width (Unitless)', 'FontSize', labelfsz)
ylabel('Current (A)', 'FontSize', labelfsz)
legend('Phase 1', 'Phase 2', 'Phase 3')
legend('boxoff')
set(gca, 'YLim', [min(min(Icoils)), max(max(Icoils))] .* 1.2, 'FontSize', 12)

subplot(3,1,3)
plot(xBh./design.Wp, limitIcoils(:,1), ':r', xBh./design.Wp, limitIcoils(:,2), 'g', xBh./design.Wp, limitIcoils(:,3), '--b');
% title('Limited Current', 'FontSize', titlefsz);
xlabel('Position Relative to Pole Width (Unitless)', 'FontSize', labelfsz)
ylabel('Current (A)', 'FontSize', labelfsz)
legend('Phase 1', 'Phase 2', 'Phase 3')
legend('boxoff')
ylimits = get(gca, 'YLim');
set(gca, 'YLim', ylimits .* 1.2, 'FontSize', 12)

%%

Fpto = 100000;
xBh = 0;
xBs = 0;
vBh = 1;
vBs = 0;

tstep = 0.001;

ii = 1;

time = 0;

while xBh < 3*design.Wp
    
    
    [EMF(ii, 1:3), Icoils(ii,1:3), limitIcoils(ii,1:3), Force(ii, 1)] = machineidealfpto_linear(design, simoptions, Fpto, xBh, xBs, vBh, vBs);
    
    xBh = xBh + (vBh*tstep);
    
    ii = ii + 1;
    
    time(ii) = time(ii-1) + tstep;
    
end

% plot the three-phase EMFs
% plot(time(1:end-1), EMF(:,1), 'r', time(1:end-1), EMF(:,2), 'y', time(1:end-1), EMF(:,3), 'b');
% ylabel('EMF (V)');

% Add the currents to the graph
figure; plot(time(1:end-1), limitIcoils(:,1), 'r', time(1:end-1), limitIcoils(:,2), '--g', time(1:end-1), limitIcoils(:,3), ':b');
ylabel(2, 'Current (A)');

% Add the force to the graph
addaxis(time(1:end-1), Force, 'k');
addaxislabel(3,'Force');
