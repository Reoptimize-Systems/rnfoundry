function evaloptions = designandevaloptions_TM(evaloptions)

    if nargin == 0
            evaloptions = [];
    end
    
    evaloptions = designandevaloptions_linear(evaloptions);
    
end