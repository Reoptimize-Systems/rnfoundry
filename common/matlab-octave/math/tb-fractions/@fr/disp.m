function disp(FRAC)
% fr/disp: Disp a fraction object
%
% arguments:
%  FRAC     - a fraction object
%
% If the fraction parts are recognized types, it displays as K + N / D
% otherwise displays the parts separately.
% MxN arrays are displayed as single columns.
%
%
%  See also: display

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release (using vpi/disp as a template)
%   28/7/09 - removed display of zero parts
%   20/8/09 - fixed bug when large vpi numbers have multiline output


nf = numel(FRAC);
sf = size(FRAC);

if nf == 0
  % an empty array
  disp('   []');
elseif nf > 1
  % an array with 2 or more elements: dump everything into rows

  fprintf(' %d',sf(1));
  fprintf('x%d',sf(2:end));
  fprintf(' fraction array with elements (reading down columns)\n');
  for i=1:nf
    disp(FRAC(i));
  end;

else
  % A scalar fraction
  
  N=FRAC.numer;
  D=FRAC.denom;
  K=FRAC.whole;
  
  cn=class(N);
  cd=class(D);
  ck=class(K);

  classok={'vpi','double','single','int8','uint8','int16', ...
    'uint16','int32','uint32','int64','uint64','logical'};
  if any(~ismember({cn,cd,ck},classok))
    % display using disp of the unknown data type
    disp(K);
    fprintf(' +\n');
    disp(N);
    fprintf(' /\n');
    disp(D);
    fprintf('\n');
  else
    % can write the fraction more compactly
    if isfinite(FRAC)
      % write as K or K + N/D or N/D, depending whether K or N are zero
      Knz=K~=0;
      Nnz=N~=0;
      if Knz || ~Nnz
        % write K if it is nonzero or if N=0
        if isequal(ck,'vpi'),
          s=['   ',disp1(K)];
        else
          s=['   ',int2str(K)];
        end;
        if Nnz
          % follow by plus sign for N/D
          s=[s,' + '];
        end;
      else
        % for single N/D term
        s = '   ';
      end;
      if Nnz
        % include N/D term
        if isequal(cn,'vpi'),
          s=[s,disp1(N)];
        else
          s=[s,int2str(N)];
        end;
        if isequal(cd,'vpi'),
          s=[s,' / ',disp1(D)];
        else
          s=[s,' / ',int2str(D)];
        end;
      end;
    else
      % inf, nan
      if isnan(FRAC)
        s='   NaN';
      elseif FRAC>0
        % +inf
        s='   Inf';
      else
        % -inf
        s='  -Inf';
      end;
    end;
    fprintf('%s\n',s);
  end;
end

function s=disp1(n)
% large vpi numbers are displayed in a char array
% - displays vpi number in 1 line and trims leading/trailing blanks
s=strtrim(reshape(disp(n)',1,[]));



