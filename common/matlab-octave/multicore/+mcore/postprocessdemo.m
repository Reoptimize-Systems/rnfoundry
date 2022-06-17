function postprocessdemo(postProcStruct)
%POSTPROCESSDEMO  Example for a user-defined postprocessing function.
%   POSTPROCESSDEMO(POSTPROCSTRUCT) accepts a structure passed by function
%   STARTMULTICOREMASTER and makes some command-line outputs. The functions
%   shall demonstrate the basic use of the postprocessing function.
%
%   Markus Buehren
%   Last modified 19.06.2009
%
%   See also STARTMULTICOREMASTER.

persistent lastDisplayTime
if isempty(lastDisplayTime)
  lastDisplayTime = mcore.mbtime;
end

if strcmp(postProcStruct.state, 'initialization')
  fprintf('Multicore demo started at %s\n', postProcStruct.userData);
end  

if mbtime - lastDisplayTime > 2.0
  lastDisplayTime = mcore.mbtime;

  fprintf('Current time: %s.\n', datestr(now));
  fprintf('Jobs done by master: %2d\n',     postProcStruct.nrOfFilesMaster);
  fprintf('Jobs done by slave:  %2d\n',     postProcStruct.nrOfFilesSlaves);

end