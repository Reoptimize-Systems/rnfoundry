function [rules,vars] = MMakefile ()

%     thisfilepath = mfemmdeps.getmfilepath (mfilename);
% 
%     if ispc
%         trilibraryflag = '-DCPU86';
%     else
%         trilibraryflag = '-DLINUX';
%     end
% 
%     % flags that will be passed direct to mex
%     vars.MEXFLAGS = ['${MEXFLAGS} -I"../cfemm/fmesher" -I"../cfemm/libfemm" -I"../cfemm/libfemm/liblua" ', trilibraryflag];
%     
%     cfemmpath = fullfile (thisfilepath, '..', 'cfemm');
%     libfemmpath = fullfile (cfemmpath, 'libfemm');
%     libluapath = fullfile (libfemmpath, 'liblua');

    vars.OBJS = { ...
      'wire_size.${OBJ_EXT}', ... 
      'wire.${OBJ_EXT}', ...
      'starofslot.${OBJ_EXT}', ... 
      'coil.${OBJ_EXT}', ... 
      'winding.${OBJ_EXT}', ... 
      'm_phase_winding.${OBJ_EXT}', ...
      'mexmPhaseWL.cpp' };

    % mexfmesher.${MEX_EXT}: ${OBJS}
    %     mex $^ -output $@
    rules(1).target = {'mexmPhaseWL.${MEX_EXT}'};
    rules(1).deps = vars.OBJS;
    rules(1).commands = 'mex ${MEXFLAGS} $^ -output $@';
    
    % created the following using:
    % clc
    % for i = 1:numel (vars.OBJS)
    %     fprintf ('rules(end+1).target = ''%s.${OBJ_EXT}'';\nrules(end).deps = ''%s.h'';\n\n', vars.OBJS{i}(1:end-11), vars.OBJS{i}(1:end-11));
    % end

    rules(end+1).target = 'wire_size.${OBJ_EXT}';
    rules(end).deps = {'wire_size.h', 'constant.h'};

    rules(end+1).target = 'wire.${OBJ_EXT}';
    rules(end).deps = 'wire.h';

    rules(end+1).target = 'starofslot.${OBJ_EXT}';
    rules(end).deps = 'starofslot.h';
    
    rules(end+1).target = 'coil.${OBJ_EXT}';
    rules(end).deps = 'coil.h';
    
    rules(end+1).target = 'winding.${OBJ_EXT}';
    rules(end).deps = 'winding.h';

    rules(end+1).target = 'm_phase_winding.${OBJ_EXT}';
    rules(end).deps = {'m_phase_winding.h', 'winding.h'};


    rules(3).target = 'tidy';
    rules(3).commands = {'try; delete(''*.${OBJ_EXT}''); catch; end;'};
    
    rules(4).target = 'clean';
    rules(4).commands = [ rules(3).commands, ...
                          {'try; delete(''*.${MEX_EXT}''); catch; end;'}];

end
