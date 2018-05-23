%%


check_wrap = @(result, expect) assert (all(strcmp (result, expect)), 'test failed');


%% simple test

text = 'Hello there, how are you this fine day?  I''m glad to hear it!';

expect = {'Hello there,', ...
          'how are you', ...
          'this fine', ...
          'day?  I''m', ...
          'glad to hear', ...
          'it!'};
      
result = TextWrapper.wraplines (text, 'width', 12);

check_wrap (result, expect);

%% 

text = 'Hello there, how are you this fine day?  I''m glad to hear it!';

expect = {'Hello there, how are you this fine day?', ...
          'I''m glad to hear it!'};
                     
result = TextWrapper.wraplines (text, 'width', 42);

check_wrap (result, expect);      
               
%%

text = 'Hello there, how are you this fine day?  I''m glad to hear it!';

result = TextWrapper.wraplines (text, 'width', 80);

check_wrap (result, {text});


%%

%Check that wrapping the empty string returns an empty list.
result = TextWrapper.wraplines ('', 'Width', 6);

check_wrap(result, {''});

result = TextWrapper.wraplines ('', 'Width', 6, 'DropWhitespace', false);
check_wrap(result, {''});

%% function tests

test_empty_string_with_initial_indent(check_wrap)
% test_whitespace(check_wrap)
% test_fix_sentence_endings(check_wrap)
test_wrap_short(check_wrap)
test_wrap_short_1line(check_wrap)
% test_hyphenated(check_wrap)
test_hyphenated_numbers(check_wrap)

%% function test definitions

function  test_empty_string_with_initial_indent(check_wrap)
        % Check that the empty string is not indented.
        
	result = TextWrapper.wraplines ('', 'Width', 6, 'InitialIndent', '++');
    
    check_wrap (result, {''});
    
    result = TextWrapper.wraplines ('', 'InitialIndent', '++', 'DropWhiteSpaceChars', false);
    
    check_wrap (result, {''});

end

% function  test_whitespace(check_wrap)
%         % Whitespace munging and end-of-sentence detection
% 
%         text = [ ...
% 'This is a paragraph that already has\n', ... 
% 'line breaks.  But some of its lines are much longer than the others,\n', ...
% 'so it needs to be wrapped.\n', ...
% 'Some lines are \ttabbed too.\n', ...
% 'What a mess!\n', ...
% ];
% 
%         expect = ['This is a paragraph that already has line', ...
%                   'breaks.  But some of its lines are much', ...
%                   'longer than the others, so it needs to be', ...
%                   'wrapped.  Some lines are  tabbed too.  What a', ...
%                   'mess!'];
% 
%         wrapper = TextWrapper(45, 'FixSentenceEndings', true);
%         result = wrapper.wrap(text);
%         self.check(result, expect);
% 
%         result = wrapper.fill(text);
%         self.check(result, '\n'.join(expect));
% 
%         text = '\tTest\tdefault\t\ttabsize.';
%         expect = ['        Test    default         tabsize.'];
%         check_wrap (text, 80, expect);
% 
%         text = '\tTest\tcustom\t\ttabsize.';
%         expect = ['    Test    custom      tabsize.'];
%         check_wrap (text, 80, expect, tabsize=4);
% end

% function  test_fix_sentence_endings(check_wrap)
% 
%         wrapper = TextWrapper(60, 'FixSentenceEndings', true);
% 
%         % SF %847346: ensure that fix_sentence_endings', true does the
%         % right thing even on input short enough that it doesn't need to
%         % be wrapped.
%         text = 'A short line. Note the single space.';
%         expect = {'A short line.  Note the single space.'};
%         self.check(wrapper.wrap(text), expect);
% 
%         % Test some of the hairy end cases that _fix_sentence_endings()
%         % is supposed to handle (the easy stuff is tested in
%         % test_whitespace() above).
%         text = 'Well, Doctor? What do you think?';
%         expect = {'Well, Doctor?  What do you think?'};
%         self.check(wrapper.wrap(text), expect);
% 
%         text = 'Well, Doctor?\nWhat do you think?';
%         self.check(wrapper.wrap(text), expect);
% 
%         text = 'I say, chaps! Anyone for ''tennis?''\nHmmph!';
%         expect = {'I say, chaps!  Anyone for ''tennis?''  Hmmph!'};
%         self.check(wrapper.wrap(text), expect);
% 
%         wrapper.width = 20;
%         expect = {'I say, chaps!', 'Anyone for ''tennis?''', 'Hmmph!'};
%         self.check(wrapper.wrap(text), expect);
% 
%         text = 'And she said, ''Go to hell!''\nCan you believe that?';
%         expect = { 'And she said, ''Go to', ...
%                    'hell!''  Can you', ...
%                    'believe that?' };
%         self.check(wrapper.wrap(text), expect);
% 
%         wrapper.width = 60;
%         expect = { 'And she said, ''Go to hell!''  Can you believe that?'};
%         self.check(wrapper.wrap(text), expect);
% 
%         text = 'File stdio.h is nice.';
%         expect = {'File stdio.h is nice.'};
%         self.check(wrapper.wrap(text), expect);
%         
% end

function  test_wrap_short(check_wrap)
        % Wrapping to make short lines longer

        text = sprintf('This is a\nshort paragraph.');

        result = TextWrapper.wraplines (text, 'Width', 20);
        
        check_wrap(result, {'This is a short', ...
                                   'paragraph.'});
                               
        result = TextWrapper.wraplines (text, 'Width', 40);

        check_wrap(result, {'This is a short paragraph.'});

end

function  test_wrap_short_1line(check_wrap)
        % Test endcases

        text = 'This is a short line.';

        result = TextWrapper.wraplines (text, 'Width', 30);
        
        check_wrap(result, {'This is a short line.'});
        
        result = TextWrapper.wraplines (text, 'Width', 30, 'InitialIndent', '(1) ');
        
        check_wrap(result, {'(1) This is a short line.'});

end

% function  test_hyphenated(check_wrap)
%         % Test breaking hyphenated words
% 
%         text = ['this-is-a-useful-feature-for-', ...
%                 'reformatting-posts-from-tim-peters''ly'];
% 
%         check_wrap(text, 40, ...
%                         {'this-is-a-useful-feature-for-', ...
%                          'reformatting-posts-from-tim-peters''ly'});
%         check_wrap(text, 41, ...
%                         {'this-is-a-useful-feature-for-', ...
%                          'reformatting-posts-from-tim-peters''ly'});
%         check_wrap(text, 42, ...
%                         {'this-is-a-useful-feature-for-reformatting-', ...
%                          'posts-from-tim-peters''ly'});
%         % The test tests current behavior but is not testing parts of the API.
%         expect = strsplit (['this-|is-|a-|useful-|feature-|for-|', ...
%                             'reformatting-|posts-|from-|tim-|peters''ly'], '|');
%         check_wrap(text, 1, expect, 'BreakLongWords', false)
%         self.check_split(text, expect);
% 
%         self.check_split('e-mail', {'e-mail'});
%         self.check_split('Jelly-O', {'Jelly-O'});
%         % The test tests current behavior but is not testing parts of the API.
%         self.check_split('half-a-crown', strsplit ('half-|a-|crown', '|'));
% end

function  test_hyphenated_numbers(check_wrap)
    % Test that hyphenated numbers (eg. dates) are not broken like words.
    text = sprintf (['Python 1.0.0 was released on 1994-01-26.  Python 1.0.1 was\n', ...
                     'released on 1994-02-15.']);

    result = TextWrapper.wraplines (text, 'Width', 30);

    check_wrap(result, {'Python 1.0.0 was released on', ...
                        '1994-01-26.  Python 1.0.1 was', ...
                        'released on 1994-02-15.'});

    result = TextWrapper.wraplines (text, 'Width', 30, 'BreakLongWords', false);

    check_wrap(result, {'Python 1.0.0 was released on 1994-01-26.', ...
                        'Python 1.0.1 was released on 1994-02-15.'});
                    
    check_wrap(result, strsplit (text));

    text = 'I do all my shopping at 7-11.';

    result = TextWrapper.wraplines (text, 'Width', 25);

    check_wrap(result, {'I do all my shopping at', ...
                          '7-11.'});
                      
    result = TextWrapper.wraplines (text, 'Width', 27);

    check_wrap(result, {'I do all my shopping at', ...
                          '7-11.'});
                      
    result = TextWrapper.wraplines (text, 'Width', 29);

    check_wrap(result, {'I do all my shopping at 7-11.'});
    
    result = TextWrapper.wraplines (text, 'Width', 1, 'BreakLongWords', false);
    
    check_wrap(result, strsplit (text));
    
end
