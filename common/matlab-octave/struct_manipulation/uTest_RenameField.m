function uTest_RenameField(doSpeed)
% Automatic test: RenameField
% This is a routine for automatic testing. It is not needed for processing and
% can be deleted or moved to a folder, where it does not bother.
%
% uTest_RenameField(doSpeed)
% INPUT:
%   doSpeed: Optional logical flag to trigger time consuming speed tests.
%            Default: TRUE. If no speed test is defined, this is ignored.
% OUTPUT:
%   On failure the test stops with an error.
%
% Tested: Matlab 6.5, 7.7, 7.8, WinXP, 32 bit
% Author: Jan Simon, Heidelberg, (C) 2009-2011 matlab.THISYEAR(a)nMINUSsimon.de

% $JRev: R-p V:015 Sum:UsDS07hXEq26 Date:12-Feb-2011 00:46:28 $
% $License: BSD (see Docs\BSD_License.txt) $
% $File: Tools\UnitTests_\uTest_RenameField.m $
% History:
% 015: 11-Feb-2011 11:54, Check cell string input.

% Initialize: ==================================================================
ErrID = ['JSim:', mfilename];

if nargin == 0
   doSpeed = true;
end

whichFun  = which('RenameField');
C10 = char(10);

disp(['==== Test RenameField:  ', datestr(now, 0), C10, ...
      'Versions: ', whichFun, C10]);

disp('== Known answer tests:');
disp('  String method');
S = [];
R = RenameField(S, 'old', 'new');
if ~(isempty(R) && isa(R, 'double'))
   error([ErrID, 'KATFailed'], 'Bad reply for S=[]');
end

S = struct([]);
R = RenameField(S, 'old', 'new');
if ~(isempty(R) && isa(R, 'struct') && isempty(fieldnames(R)))
   error([ErrID, 'KATFailed'], 'Bad reply for S=struct([])');
end

S = struct('a', {});
R = RenameField(S, 'old', 'new');
if ~(isempty(R) && isa(R, 'struct') && isequal(fieldnames(R), {'a'}))
   error([ErrID, 'KATFailed'], 'Bad reply for empty S.a');
end

S = struct('a', 1);
R = RenameField(S, 'old', 'new');
if ~(isa(R, 'struct') && isequal(fieldnames(R), {'a'}))
   error([ErrID, 'KATFailed'], 'Bad reply for S.a');
end

S = struct('a', {1, 2});
R = RenameField(S, 'old', 'new');
if ~(isa(R, 'struct') && isequal(fieldnames(R), {'a'}) && ....
      isequal(size(R), [1, 2]))
   error([ErrID, 'KATFailed'], 'Bad reply for empty S(1:2).a');
end

S = struct('a', {}, 'b', {});
R = RenameField(S, 'old', 'new');
if ~(isempty(R) && isa(R, 'struct') && isequal(fieldnames(R), {'a'; 'b'}))
   error([ErrID, 'KATFailed'], 'Bad reply for empty S.a, S.b');
end

S = struct('a', 1, 'b', 2);
R = RenameField(S, 'old', 'new');
if ~isequal(R, S)
   error([ErrID, 'KATFailed'], 'Bad reply for scalar S.a, S.b');
end

S = struct('a', {1, 2}, 'b', {3, 4});
R = RenameField(S, 'old', 'new');
if ~isequal(R, S)
   error([ErrID, 'KATFailed'], 'Bad reply for S(1:2).a, S(1:2).b');
end

% ----------------------
S = struct('old', {});
R = RenameField(S, 'old', 'new');
if ~(isempty(R) && isa(R, 'struct') && isequal(fieldnames(R), {'new'}))
   error([ErrID, 'KATFailed'], 'Bad reply for empty S.old');
end

S = struct('old', 1);
R = RenameField(S, 'old', 'new');
if ~(isa(R, 'struct') && isequal(fieldnames(R), {'new'}))
   error([ErrID, 'KATFailed'], 'Bad reply for S.old');
end

S = struct('old', {1, 2});
R = RenameField(S, 'old', 'new');
if ~(isa(R, 'struct') && isequal(fieldnames(R), {'new'}) && ....
      isequal(size(R), [1, 2]))
   error([ErrID, 'KATFailed'], 'Bad reply for empty S(1:2).old');
end

S = struct('old', {}, 'b', {});
R = RenameField(S, 'old', 'new');
if ~(isempty(R) && isa(R, 'struct') && isequal(fieldnames(R), {'new'; 'b'}))
   error([ErrID, 'KATFailed'], 'Bad reply for empty S.old, S.b');
end

S = struct('old', 1, 'b', 2);
R = RenameField(S, 'old', 'new');
if ~isequal(R, struct('new', 1, 'b', 2))
   error([ErrID, 'KATFailed'], 'Bad reply for scalar S.old, S.b');
end

S = struct('old', {1, 2}, 'b', {3, 4});
R = RenameField(S, 'old', 'new');
if ~isequal(R, struct('new', {1, 2}, 'b', {3, 4}))
   error([ErrID, 'KATFailed'], 'Bad reply for S(1:2).old, S(1:2).b');
end

% ---------------------------
S = struct('a', {}, 'old', {});
R = RenameField(S, 'old', 'new');
if ~(isempty(R) && isa(R, 'struct') && isequal(fieldnames(R), {'a'; 'new'}))
   error([ErrID, 'KATFailed'], 'Bad reply for empty S.a, S.old');
end

S = struct('a', 1, 'old', 2);
R = RenameField(S, 'old', 'new');
if ~isequal(R, struct('a', 1, 'new', 2))
   error([ErrID, 'KATFailed'], 'Bad reply for scalar S.a, S.old');
end

S = struct('a', {1, 2}, 'old', {3, 4});
R = RenameField(S, 'old', 'new');
if ~isequal(R, struct('a', {1, 2}, 'new', {3, 4}))
   error([ErrID, 'KATFailed'], 'Bad reply for S(1:2).a, S(1:2).old');
end

S = struct('a', {1, 2}, 'new', {3, 4});
R = RenameField(S, 'new', 'new');
if ~isequal(R, S)
   error([ErrID, 'KATFailed'], 'Bad reply for S(1:2).a, S(1:2).new');
end

disp('  ok');

% ---------------------------
disp('  Cell string method');
S = [];
R = RenameField(S, {'old'}, {'new'});
if ~(isempty(R) && isa(R, 'double'))
   error([ErrID, 'KATFailed'], 'Bad reply for S=[]');
end

R = RenameField(S, {}, {});
if ~isequal(R, S)
   error([ErrID, 'KATFailed'], 'Bad reply for S=[], New={}');
end

S = struct([]);
R = RenameField(S, {'old'}, {'new'});
if ~(isempty(R) && isa(R, 'struct') && isempty(fieldnames(R)))
   error([ErrID, 'KATFailed'], 'Bad reply for S=struct([])');
end

S = struct('a', {});
R = RenameField(S, {'old'}, {'new'});
if ~(isempty(R) && isa(R, 'struct') && isequal(fieldnames(R), {'a'}))
   error([ErrID, 'KATFailed'], 'Bad reply for empty S.a');
end

R = RenameField(S, {}, {});
if ~isequal(R, S)
   error([ErrID, 'KATFailed'], 'Bad reply for empty S.a and New={}');
end

S = struct('a', 1);
R = RenameField(S, {'old'}, {'new'});
if ~(isa(R, 'struct') && isequal(fieldnames(R), {'a'}))
   error([ErrID, 'KATFailed'], 'Bad reply for S.a');
end

R = RenameField(S, {}, {});
if ~(isa(R, 'struct') && isequal(fieldnames(R), {'a'}))
   error([ErrID, 'KATFailed'], 'Bad reply for S.a and New={}');
end

S = struct('a', {1, 2});
R = RenameField(S, {'old'}, {'new'});
if ~(isa(R, 'struct') && isequal(fieldnames(R), {'a'}) && ....
      isequal(size(R), [1, 2]))
   error([ErrID, 'KATFailed'], 'Bad reply for empty S(1:2).a');
end

S = struct('a', {}, 'b', {});
R = RenameField(S, {'old'}, {'new'});
if ~(isempty(R) && isa(R, 'struct') && isequal(fieldnames(R), {'a'; 'b'}))
   error([ErrID, 'KATFailed'], 'Bad reply for empty S.a, S.b');
end

S = struct('a', 1, 'b', 2);
R = RenameField(S, {'old'}, {'new'});
if ~isequal(R, S)
   error([ErrID, 'KATFailed'], 'Bad reply for scalar S.a, S.b');
end

S = struct('a', {1, 2}, 'b', {3, 4});
R = RenameField(S, {'old'}, {'new'});
if ~isequal(R, S)
   error([ErrID, 'KATFailed'], 'Bad reply for S(1:2).a, S(1:2).b');
end

R = RenameField(S, {}, {});
if ~isequal(R, S)
   error([ErrID, 'KATFailed'], 'Bad reply for S(1:2).a, S(1:2).b and New={}');
end

% ----------------------
S = struct('old', {});
R = RenameField(S, {'old'}, {'new'});
if ~(isempty(R) && isa(R, 'struct') && isequal(fieldnames(R), {'new'}))
   error([ErrID, 'KATFailed'], 'Bad reply for empty S.old');
end

S = struct('old', 1);
R = RenameField(S, {'old'}, {'new'});
if ~(isa(R, 'struct') && isequal(fieldnames(R), {'new'}))
   error([ErrID, 'KATFailed'], 'Bad reply for S.old');
end

S = struct('old', {1, 2});
R = RenameField(S, {'old'}, {'new'});
if ~(isa(R, 'struct') && isequal(fieldnames(R), {'new'}) && ....
      isequal(size(R), [1, 2]))
   error([ErrID, 'KATFailed'], 'Bad reply for empty S(1:2).old');
end

S = struct('old', {}, 'b', {});
R = RenameField(S, {'old'}, {'new'});
if ~(isempty(R) && isa(R, 'struct') && isequal(fieldnames(R), {'new'; 'b'}))
   error([ErrID, 'KATFailed'], 'Bad reply for empty S.old, S.b');
end

S = struct('old', 1, 'b', 2);
R = RenameField(S, {'old'}, {'new'});
if ~isequal(R, struct('new', 1, 'b', 2))
   error([ErrID, 'KATFailed'], 'Bad reply for scalar S.old, S.b');
end

S = struct('old', {1, 2}, 'b', {3, 4});
R = RenameField(S, {'old'}, {'new'});
if ~isequal(R, struct('new', {1, 2}, 'b', {3, 4}))
   error([ErrID, 'KATFailed'], 'Bad reply for S(1:2).old, S(1:2).b');
end

% ---------------------------
S = struct('a', {}, 'old', {});
R = RenameField(S, {'old'}, {'new'});
if ~(isempty(R) && isa(R, 'struct') && isequal(fieldnames(R), {'a'; 'new'}))
   error([ErrID, 'KATFailed'], 'Bad reply for empty S.a, S.old');
end

S = struct('a', 1, 'old', 2);
R = RenameField(S, {'old'}, {'new'});
if ~isequal(R, struct('a', 1, 'new', 2))
   error([ErrID, 'KATFailed'], 'Bad reply for scalar S.a, S.old');
end

S = struct('a', {1, 2}, 'old', {3, 4});
R = RenameField(S, {'old'}, {'new'});
if ~isequal(R, struct('a', {1, 2}, 'new', {3, 4}))
   error([ErrID, 'KATFailed'], 'Bad reply for S(1:2).a, S(1:2).old');
end

S = struct('a', {1, 2}, 'new', {3, 4});
R = RenameField(S, {'new'}, {'new'});
if ~isequal(R, S)
   error([ErrID, 'KATFailed'], 'Bad reply for S(1:2).a, S(1:2).new');
end

S = struct('a', {1, 2}, 'b', {3, 4});
R = RenameField(S, {'a', 'b'}, {'a2', 'b2'});
Want = struct('a2', {1, 2}, 'b2', {3, 4});
if ~isequal(R, Want)
   error([ErrID, 'KATFailed'], ...
      'Bad reply for S(1:2).a, S(1:2).old, {a,b}->{a2,b2}');
end

R = RenameField(S, {'b', 'a'}, {'b2', 'a2'});
if ~isequal(R, Want)
   error([ErrID, 'KATFailed'], ...
      'Bad reply for S(1:2).a, S(1:2).old, {b,a}->{b2,a2}');
end

R = RenameField(S, {'b', 'a', 'notExisting'}, {'b2', 'a2', 'notExisting2'});
if ~isequal(R, Want)
   error([ErrID, 'KATFailed'], ...
      'Bad reply for S(1:2).a, S(1:2).old, {b,a,miss}->{b2,a2,miss2}');
end

R = RenameField(S, {'b', 'notExisting', 'a'}, {'b2', 'notExisting2', 'a2'});
if ~isequal(R, Want)
   error([ErrID, 'KATFailed'], ...
      'Bad reply for S(1:2).a, S(1:2).old, {b,miss,a}->{b2,miss2,a2}');
end

R = RenameField(S, {'notExisting', 'b', 'a'}, {'notExisting2', 'b2', 'a2'});
if ~isequal(R, Want)
   error([ErrID, 'KATFailed'], ...
      'Bad reply for S(1:2).a, S(1:2).old, {miss,b,a}->{miss2,b2,a2}');
end

disp('  ok');

% ---------------------------
disp('== Provoke errors:');
% Length of the name is checked even with disabled validity checks, because
% Matlab's memory management gets corrupted for a too long name!
S = struct('a', 1, 'b', 2, 'c', 3);
tooLazy = false;
try
   R       = RenameField(S, 'a', repmat('a', 1, 64));
   tooLazy = true;
catch
   % Nothing to do
end
if tooLazy
   error([ErrID, 'InvalidNameAccepted'], '64 characters accepted!');
end

S = struct('a', 1, 'b', 2, 'c', 3);
try
   R = RenameField(S, 'a', '_');
   disp('  dangerous: No check for valid field name!');
   CheckEnabled = false;
catch
   disp('  ok: Check for valid field name is enabled!');
   CheckEnabled = true;
end

if CheckEnabled
   try
      R       = RenameField(S, 'a', '0');
      tooLazy = false;
   catch
      % Nothing to do
   end
   if tooLazy
      error([ErrID, 'InvalidNameAccepted'], 'Duplicate name accepted!');
   end
   
   try
      R       = RenameField(S, 'a', '');
      tooLazy = false;
   catch
      % Nothing to do
   end
   if tooLazy
      error([ErrID, 'InvalidNameAccepted'], 'Empty name accepted!');
   end
   
   try
      R       = RenameField(S, {'a'}, {'_a'});
      tooLazy = false;
   catch
      % Nothing to do
   end
   if tooLazy
      error([ErrID, 'InvalidNameAccepted'], 'Leading underscore accepted!');
   end

   try
      R       = RenameField(S, 'a', '0a');
      tooLazy = false;
   catch
      % Nothing to do
   end
   if tooLazy
      error([ErrID, 'InvalidNameAccepted'], 'Leading number accepted!');
   end
   
   try
      R       = RenameField(S, 'a', [98, 99]);
      tooLazy = false;
   catch
      % Nothing to do
   end
   if tooLazy
      error([ErrID, 'InvalidNameAccepted'], 'DOUBLE accepted!');
   end
   
   try
      R       = RenameField(S, 'a', 'b_', 'c_');
      tooLazy = false;
   catch
      % Nothing to do
   end
   if tooLazy
      error([ErrID, 'InvalidNameAccepted'], '4 inputs accepted!');
   end
   
   try
      R       = RenameField(S, 'a', {'b'});
      tooLazy = false;
   catch
      % Nothing to do
   end
   if tooLazy
      error([ErrID, 'InvalidNameAccepted'], 'String+Cell input accepted!');
   end
   
   try
      R       = RenameField(S, {'a'}, 'b');
      tooLazy = false;
   catch
      % Nothing to do
   end
   if tooLazy
      error([ErrID, 'InvalidNameAccepted'], 'Cell+String input accepted!');
   end
   
   try
      R       = RenameField(S, {'a'}, {'b', 'c'});
      tooLazy = false;
   catch
      % Nothing to do
   end
   if tooLazy
      error([ErrID, 'InvalidNameAccepted'], ...
         'Cell input with different sizes accepted!');
   end

   try
      R       = RenameField(S, cell(1, 2), {'b', 'c'});
      tooLazy = false;
   catch
      % Nothing to do
   end
   if tooLazy
      error([ErrID, 'InvalidNameAccepted'], ...
         'Unitialized cell elements in [Old] accepted!');
   end

   try
      R       = RenameField(S, {'a', 'b'}, cell(1, 2));
      tooLazy = false;
   catch
      % Nothing to do
   end
   if tooLazy
      error([ErrID, 'InvalidNameAccepted'], ...
         'Unitialized cell elements in [New] accepted!');
   end
end

disp('  ok');

% Speed: -----------------------------------------------------------------------
if doSpeed
   TestTime = 1.0;  % sec
   fprintf('== Speed test (test time: %g sec):\n', TestTime);
else
   TestTime = 0.1;
   fprintf('\n== Speed test (test time: %g sec - may be inaccurate):\n', ...
      TestTime);
end

% Determine number of loops:
C = cell(1000, 1);
D = cell(1000, 1);
E = cell(1000, 1);
for iC = 1:1000
   C{iC} = sprintf('field%d', iC);
   D{iC} = rand(8, 8);
   E{iC} = sprintf('new%d', iC);
end

S1000 = cell2struct(D, C);
iTime = cputime;
iLoop = 0;
while cputime - iTime < TestTime
   v = RenameField(S1000, 'field500', 'newname');  %#ok<*NASGU>
   clear('v');   % Suppress JIT acceleration for realistic times
   iLoop = iLoop + 1;
end
nDigit = max(1, floor(log10(max(1, iLoop))) - 1);
nLoop  = max(2, round(iLoop / 10 ^ nDigit) * 10 ^ nDigit);

fprintf('STRING - rename 1 field of a struct:  (%d loops)\n', nLoop);
for nField = [1, 10, 100, 1000]
   fprintf('  Struct has %d fields:\n', ...
      nField);
   pause(0.02);  % Flush Java event queue to show the message
   
   S = cell2struct(D(1:nField), C(1:nField));
   tic;
   for iLoop = 1:nLoop
      T = RenameField_local(S, C{rem(iLoop, nField) + 1}, 'newname');
      clear('T');
   end
   M_time = toc + eps;
   
   tic;
   for iLoop = 1:nLoop
      T = RenameField(S, C{rem(iLoop, nField) + 1}, 'newname');
      clear('T');
   end
   Mex_time = toc;
   
   fprintf('    M:  %7.3f sec   Mex:%7.3f sec  ==> %4.1f%% of M-version\n', ...
      M_time, Mex_time, 100 * Mex_time / M_time);
end

nLoop = nLoop * 10;
fprintf('\nCELL STRING - rename all fields of a struct:  (%d loops)\n', nLoop);
for nField = [1, 10, 100, 1000]
   fprintf('  Struct has %d fields:\n', nField);
   pause(0.02);  % Flush Java event queue to show the message
   
   C_ = C(1:nField);
   E_ = E(1:nField);
   
   S = cell2struct(D(1:nField), C_);
   tic;
   for iLoop = 1:nLoop / nField
      T = RenameFieldC_local(S, C_, E_);
      clear('T');
   end
   M_time = toc + eps;
   
   tic;
   for iLoop = 1:nLoop / nField
      T = RenameField(S, C_, E_);
      clear('T');
   end
   Mex_time = toc;
   fprintf('    M:  %7.3f sec   Mex:%7.3f sec  ==> %4.1f%% of M-version\n', ...
      M_time, Mex_time, 100 * Mex_time / M_time);
end

fprintf('\n== RenameField passed the tests.\n');

return;

% ==============================================================================
function S = RenameField_local(S, Old, New)
% A simple and fairly efficient M-version:

Data  = struct2cell(S);
Field = fieldnames(S);
Field{strcmp(Field, Old)} = New;
S     = cell2struct(Data, Field);

return;

% ==============================================================================
function S = RenameFieldC_local(S, Old, New)
% A simple and fairly efficient M-version:

Data  = struct2cell(S);
Field = fieldnames(S);
for i = 1:numel(Old)
   Field{strcmp(Field, Old{i})} = New{i};
end
S     = cell2struct(Data, Field);

return;
