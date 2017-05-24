%%


check_wrap = @(result, expect) all(strcmp (result, expect));


%% simple test

text = 'Hello there, how are you this fine day?  I''m glad to hear it!';

expect = {'Hello there,', ...
          'how are you', ...
          'this fine', ...
          'day?  I''m', ...
          'glad to hear', ...
          'it!'};
      
result = TextWrapper.wraplines (text, 'width', 12);

assert (check_wrap (result, expect));

%% 

text = 'Hello there, how are you this fine day?  I''m glad to hear it!';

expect = {'Hello there, how are you this fine day?', ...
          'I''m glad to hear it!'};
                     
result = TextWrapper.wraplines (text, 'width', 42);

assert (check_wrap (result, expect));      
               
%%

text = 'Hello there, how are you this fine day?  I''m glad to hear it!';

result = TextWrapper.wraplines (text, 'width', 80);

assert (check_wrap (result, {text}));


%%

%Check that wrapping the empty string returns an empty list.
result = TextWrapper.wraplines ('', 'width', 6);

        self.check_wrap('', 6, [])
        
        
        self.check_wrap('', 6, [], drop_whitespace=False)

    def test_empty_string_with_initial_indent(self):
        # Check that the empty string is not indented.
        self.check_wrap('', 6, [], initial_indent='++')
        self.check_wrap('', 6, [], initial_indent='++', drop_whitespace=False)

    def test_whitespace(self):
        # Whitespace munging and end-of-sentence detection

        text = '''\
This is a paragraph that already has
line breaks.  But some of its lines are much longer than the others,
so it needs to be wrapped.
Some lines are \ttabbed too.
What a mess!
'''

        expect = ['This is a paragraph that already has line',
                  'breaks.  But some of its lines are much',
                  'longer than the others, so it needs to be',
                  'wrapped.  Some lines are  tabbed too.  What a',
                  'mess!']

        wrapper = TextWrapper(45, fix_sentence_endings=True)
        result = wrapper.wrap(text)
        self.check(result, expect)

        result = wrapper.fill(text)
        self.check(result, '\n'.join(expect))

        text = '\tTest\tdefault\t\ttabsize.'
        expect = ['        Test    default         tabsize.']
        self.check_wrap(text, 80, expect)

        text = '\tTest\tcustom\t\ttabsize.'
        expect = ['    Test    custom      tabsize.']
        self.check_wrap(text, 80, expect, tabsize=4)

    def test_fix_sentence_endings(self):
        wrapper = TextWrapper(60, fix_sentence_endings=True)

        # SF #847346: ensure that fix_sentence_endings=True does the
        # right thing even on input short enough that it doesn't need to
        # be wrapped.
        text = 'A short line. Note the single space.'
        expect = ['A short line.  Note the single space.']
        self.check(wrapper.wrap(text), expect)

        # Test some of the hairy end cases that _fix_sentence_endings()
        # is supposed to handle (the easy stuff is tested in
        # test_whitespace() above).
        text = 'Well, Doctor? What do you think?'
        expect = ['Well, Doctor?  What do you think?']
        self.check(wrapper.wrap(text), expect)

        text = 'Well, Doctor?\nWhat do you think?'
        self.check(wrapper.wrap(text), expect)

        text = 'I say, chaps! Anyone for 'tennis?'\nHmmph!'
        expect = ['I say, chaps!  Anyone for 'tennis?'  Hmmph!']
        self.check(wrapper.wrap(text), expect)

        wrapper.width = 20
        expect = ['I say, chaps!', 'Anyone for 'tennis?'', 'Hmmph!']
        self.check(wrapper.wrap(text), expect)

        text = 'And she said, 'Go to hell!'\nCan you believe that?'
        expect = ['And she said, 'Go to',
                  'hell!'  Can you',
                  'believe that?']
        self.check(wrapper.wrap(text), expect)

        wrapper.width = 60
        expect = ['And she said, 'Go to hell!'  Can you believe that?']
        self.check(wrapper.wrap(text), expect)

        text = 'File stdio.h is nice.'
        expect = ['File stdio.h is nice.']
        self.check(wrapper.wrap(text), expect)

    def test_wrap_short(self):
        # Wrapping to make short lines longer

        text = 'This is a\nshort paragraph.'

        self.check_wrap(text, 20, ['This is a short',
                                   'paragraph.'])
        self.check_wrap(text, 40, ['This is a short paragraph.'])


    def test_wrap_short_1line(self):
        # Test endcases

        text = 'This is a short line.'

        self.check_wrap(text, 30, ['This is a short line.'])
        self.check_wrap(text, 30, ['(1) This is a short line.'],
                        initial_indent='(1) ')


    def test_hyphenated(self):
        # Test breaking hyphenated words

        text = ('this-is-a-useful-feature-for-'
                'reformatting-posts-from-tim-peters'ly')

        self.check_wrap(text, 40,
                        ['this-is-a-useful-feature-for-',
                         'reformatting-posts-from-tim-peters'ly'])
        self.check_wrap(text, 41,
                        ['this-is-a-useful-feature-for-',
                         'reformatting-posts-from-tim-peters'ly'])
        self.check_wrap(text, 42,
                        ['this-is-a-useful-feature-for-reformatting-',
                         'posts-from-tim-peters'ly'])
        # The test tests current behavior but is not testing parts of the API.
        expect = ('this-|is-|a-|useful-|feature-|for-|'
                  'reformatting-|posts-|from-|tim-|peters'ly').split('|')
        self.check_wrap(text, 1, expect, break_long_words=False)
        self.check_split(text, expect)

        self.check_split('e-mail', ['e-mail'])
        self.check_split('Jelly-O', ['Jelly-O'])
        # The test tests current behavior but is not testing parts of the API.
        self.check_split('half-a-crown', 'half-|a-|crown'.split('|'))

    def test_hyphenated_numbers(self):
        # Test that hyphenated numbers (eg. dates) are not broken like words.
        text = ('Python 1.0.0 was released on 1994-01-26.  Python 1.0.1 was\n'
                'released on 1994-02-15.')

        self.check_wrap(text, 30, ['Python 1.0.0 was released on',
                                   '1994-01-26.  Python 1.0.1 was',
                                   'released on 1994-02-15.'])
        self.check_wrap(text, 40, ['Python 1.0.0 was released on 1994-01-26.',
                                   'Python 1.0.1 was released on 1994-02-15.'])
        self.check_wrap(text, 1, text.split(), break_long_words=False)

        text = 'I do all my shopping at 7-11.'
        self.check_wrap(text, 25, ['I do all my shopping at',
                                   '7-11.'])
        self.check_wrap(text, 27, ['I do all my shopping at',
                                   '7-11.'])
        self.check_wrap(text, 29, ['I do all my shopping at 7-11.'])
        self.check_wrap(text, 1, text.split(), break_long_words=False)