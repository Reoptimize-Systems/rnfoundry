function str = latexsipre (val,sgf,trz)

    switch nargin
        
        case 1
            str = sipre (val);
            
        case 2
            str = sipre (val,sgf);
            
        case 3
            str = sipre (val,sgf,false,trz);
            
    end
    
    switch str (end)
        
        case 'y'
            
            str = [str(1:end-2), '~y'];
        
        case 'z'
            str = [str(1:end-2), '~z'];
        
        case 'a'
            str = [str(1:end-2), '~a'];
            
        case 'f'
            str = [str(1:end-2), '~f'];
        
        case 'p'
            str = [str(1:end-2), '~p'];
            
        case 'n'
            str = [str(1:end-2), '~n'];
            
        case 'u'
            str = [str(1:end-2), '~$\mu$'];
        
        case 'm'
            str = [str(1:end-2), '~m'];
        
        case 'k'
            str = [str(1:end-2), '~k'];
        
        case 'M'
            str = [str(1:end-2), '~M'];
        
        case 'G'
            str = [str(1:end-2), '~G'];
            
        case 'T'
            str = [str(1:end-2), '~T'];
        
        case 'P'
            str = [str(1:end-2), '~P'];
            
        case 'E'
            str = [str(1:end-2), '~E'];
            
        case 'Z'
            str = [str(1:end-2), '~Z'];
        
        case 'Y'
            str = [str(1:end-2), '~Y'];
            
            
    end
end