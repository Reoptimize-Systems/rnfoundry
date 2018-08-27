function [rules,vars] = MMakefile (varargin)

    options.DoCrossBuildWin64 = false;
    options.W64CrossBuildMexLibsDir = '';
    options.Verbose = false;
%     options.Debug = false;
    
    options = parse_pv_pairs (options, varargin);
    
    vars.OBJS = { ...
      'wire_size.${OBJ_EXT}', ... 
      'wire.${OBJ_EXT}', ...
      'starofslot.${OBJ_EXT}', ... 
      'coil.${OBJ_EXT}', ... 
      'winding.${OBJ_EXT}', ... 
      'm_phase_winding.${OBJ_EXT}', ...
      'mexmPhaseWL.cpp' };
  
    extra_mex_args = '';
  
    if options.DoCrossBuildWin64 
        
        assert ( ~isempty (options.W64CrossBuildMexLibsDir), ...
                 TextWrapper.wraplines ( ['W64CrossBuild is true, but W64CrossBuildMexLibsDir is empty. ', ...
                                          'You must supply the location of the windows mingw64 windows libraries'] ) ...
               );
           
        assert ( all ( [ exist(fullfile (options.W64CrossBuildMexLibsDir, 'libmex.a'), 'file'), ...
                         exist(fullfile (options.W64CrossBuildMexLibsDir, 'libmx.a'), 'file'), ...
                         exist(fullfile (options.W64CrossBuildMexLibsDir, 'libmat.a'), 'file') ] ), ...
                 'One of libmex.a, libmx.a or libmat.a was not found in %s', ...
                 options.W64CrossBuildMexLibsDir );

        cross_full_path = mmake.cross.cross_prefix_full_path ('x86_64-w64-mingw32.static');
        
        vars.COMPILER = ['"', cross_full_path, '-gcc"'];
        
        extra_mex_args = ['-L"', options.W64CrossBuildMexLibsDir, '"'];
    end
    
    if options.Verbose
        extra_mex_args = [extra_mex_args, ' -v'];
    end

    % mexfmesher.${MEX_EXT}: ${OBJS}
    %     mex $^ -output $@
    if options.DoCrossBuildWin64
        rules(1).target = {'mexmPhaseWL.mexw64'};
        rules(1).commands = sprintf ('mex ${MEXFLAGS} ${COMPILERKEY}="${COMPILER}" ${OPTIMFLAGSKEY}="${OPTIMFLAGS}" ${CXXFLAGSKEY}="${CXXFLAGS}" ${LDFLAGSKEY}="${LDFLAGS}" %s $^ EXE="mexmPhaseWL.mexw64"', extra_mex_args);
    else
        rules(1).target = {'mexmPhaseWL.${MEX_EXT}'};
        if isoctave ()
            rules(1).commands = sprintf ('mex ${MEXFLAGS} %s $^ --output $@', extra_mex_args);
        else
            rules(1).commands = sprintf ('mex ${MEXFLAGS} ${COMPILERKEY}="${COMPILER}" ${OPTIMFLAGSKEY}="${OPTIMFLAGS}" ${CXXFLAGSKEY}="${CXXFLAGS}" ${LDFLAGSKEY}="${LDFLAGS}" %s $^ -output $@', extra_mex_args);
        end
    end
    rules(1).deps = vars.OBJS;
    
    
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
    if options.DoCrossBuildWin64
        rules(4).commands = [ rules(3).commands, ...
                              {'try; delete(''*.mexw64''); catch; end;'}];        
    else
        rules(4).commands = [ rules(3).commands, ...
                              {'try; delete(''*.${MEX_EXT}''); catch; end;'}];
    end

end
