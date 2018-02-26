% Example usage of the logger class.

% ==========================================================================
% Author: 	 Richard Crozier
% Organization:  Massachusetts Institute of Technology
% Contact:       <pavan_m@mit.edu>
% Created:       Oct 01, 2011 
% ==========================================================================

n_max = 100;


lgr = wsim.logger ();

lgr.addVariable ('time', [1,1], 'PreallocateStorage', n_max );

lgr.addVariable ( 'weight', [1,1], ...
                  'Description', 'Weight of Subjects', ...
                  'PreallocateStorage', n_max, ...
                  'Indep', 'time' );
              
lgr.addVariable ( 'height', [1,1], ...
                  'Description', 'Height of Subjects', ...
                  'PreallocateStorage', n_max, ...
                  'Indep', 'time' );
              
lgr.addVariable ( 'inertia', [3,3], ...
                  'Description', 'Inertia of Subjects', ...
                  'PreallocateStorage', n_max, ...
                  'Indep', 'time' );

% add a variable, but we don't bother preallocating this one, it will grow
% as necessary
lgr.addVariable ('sum', [1,1]);

tic;
sumval = 0;
for i = 1:n_max
    
	pause(0.001);
    
	weight = 10*rand;
    
	height = 1.5*weight + 5*rand;


	% Log my_output_1 variable as 'weight'
    lgr.logVal('time',i);
	lgr.logVal('weight', weight);
    lgr.logVal('height', height);
    lgr.logVal('inertia', weight * eye (3));

	sumval = sumval + weight + height;

% 	s.r1 = weight;
% 	s.r2 = height;
% 	s.r3 = 'Some Results';

	if weight > 6
		lgr.logWarn(['Height exceeded 6 at ' num2str(i)], true);
	end


	% Store every 10th
	if mod(i,10) == 0
		lgr.logMesg([ num2str(i) ' entries logged with sum = ' num2str(sumval)], 1);
		
		lgr.logVal('sum',sumval);
% 		lgr.logObj('result',s);
	end

end


lgr.setDefaultDesc('Subject ID');

info = lgr.getInfo ({'weight', 'in'});

try
    info = lgr.getInfo ('fjdsaf');
catch
    fprintf ('bad variable name caught by getInfo\n');
end

%%
lgr.setSeries ('inertia', lgr.data.('inertia') .* 10);

%%
try
    lgr.setSeries ('inertia', lgr.data.('weight') .* 10);
catch me
    fprintf ('bad variable size caught by setSeries, error message: \n%s\n', me.message);
end

%%

lgr.addVariable ( 'test1', [1,2], ...
                  'Description', 'test for setSeries without indep var', ...
                  'PreallocateStorage', n_max );

lgr.setSeries ('test1', repmat ([1,1], [20, 1]));

try
    lgr.setSeries ('test1', repmat ([1,1,1], [20, 1]));
catch me
    fprintf ('bad variable size caught by setSeries, error message: \n%s\n', me.message);
end

try
    lgr.setSeries ('test1', repmat ([1,1,1; 1,1,1], [1, 1, 20]));
catch me
    fprintf ('bad variable size caught by setSeries, error message: \n%s\n', me.message);
end

% 
% % Plot several variables from the logger object.
% figure; lgr.plot2Vars('weight','height','LineWidth', 2, 'Color','r'); 
% figure; lgr.plotVar('weight','LineWidth', 2, 'Color','r'); 
% figure; lgr.plotVar('height','LineWidth', 2, 'Color','r'); 
% figure; lgr.plotFofVar('weight',@sin, 'LineWidth', 2, 'Color','r'); 
% figure; lgr.plotFofVar('height',@log, 'LineWidth', 2, 'Color','r'); 
% figure; lgr.plotFofVar('time',@cumsum, 'LineWidth', 2, 'Color','r'); 
