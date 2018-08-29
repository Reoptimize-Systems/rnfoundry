classdef TextWrapper
% Port of the python TextWrapper class for wrapping long lines of text
%
% Syntax
%
%
% twobj = TextWrapper ()
% twobj = TextWrapper ('Parameter', Value)
%
% Description
%
% TexWrapper is a class for wrapping/filling text. Text can be wrapped to
% given line width. See the help for the TexWrapper constructor for details
% of it's use. The static methods wraplines and wraptext are provided as
% convenince functions and can be used like normal Matlab funcitons (they
% set up a TexWrapper class internally).
%
% TextWrapper Methods:
%
%   TextWrapper - TextWrapper constructor
%   fill - wrap paragraph to desired width and return as a character vector
%   wrap - wrap paragraph to desired width and return lines in cell array
%   wraplines - wrap paragraph to desired width and returns lines in cell array
%   wraptext - wrap paragraph to desired width and returns as a character vector
%
%

%  
% Copyright (C) 1999-2001 Gregory P. Ward.
% Copyright (C) 2002, 2003 Python Software Foundation.
% Written by Greg Ward <gward@python.net>
% Ported to m-code by Richard Crozier
%

    properties
        
        % Hardcode the recognized whitespace characters to the US-ASCII
        % whitespace characters.  The main reason for doing this is that in
        % ISO-8859-1, 0xa0 is non-breaking whitespace, so in certain locales
        % that character winds up in string.whitespace.  Respecting
        % string.whitespace in those cases would 1) make textwrap treat 0xa0 the
        % same as any other whitespace char, which is clearly wrong (it's a
        % *non-breaking* space), 2) possibly cause problems with Unicode,
        % since 0xa0 is not in range(128).
        whitespacechars = { sprintf('\t'), sprintf('\n'), sprintf('\x0b'), sprintf('\x0c'), sprintf('\r') };
        
        wordsep_re;
        wordsep_simple_re;
        
        width;
        initialIndent;
        subsequentIndent;
        expandTabs;
        ReplaceWhiteSpaceChars;
        fixSentenceEndings;
        breakLongWords;
        dropWhiteSpaceChars;
        breakOnHyphens;
        sentence_end_re;
        
    end
    
% beginning to implement dedent funciton (see bottom of this file)    
%     properties (Constant)
%         whitespacechars_only_re = re.compile('^[ \t]+$', re.MULTILINE)
%         leadingwhitespacechars_re = re.compile('(^[ \t]*)(?:[^ \t\n])', re.MULTILINE)    
%     end
   

	methods 

        function self = TextWrapper (varargin)
            % TextWrapper constructor
            %
            % Syntax
            %
            % twobj = TextWrapper ()
            % twobj = TextWrapper ('Parameter', Value)
            %
            % Description
            %
            % TexWrapper is a class for wrapping/filling text. Text can be
            % wrapped to given line width.
            %
            % Input
            %
            % Arguments are supplied as parameter-value pairs. The
            % available options are:
            %
            %  'Width' - the maximum width of wrapped lines (unless
            %    breakLongWords is false). Default is 70.
            %
            %  'InitialIndent' - character vector that will be prepended to
            %    the first line of wrapped output; Counts towards each
            %    line's width. Default is ''.
            %
            %  'SubsequentIndent' - character vector that will be prepended to
            %    the all lines of wrapped output except the first line.
            %    Counts towards each line's width. Default is ''.
            %
            %  'ExpandTabs' - true/false flag indicating whether to expand 
            %    tabs in input text to spaces before further processing.
            %    Each tab will become 1 .. 8 spaces, depending on its
            %    position in its line. If false, each tab is treated as a
            %    single character. Default is true.
            %
            %  'ReplaceWhiteSpaceChars' - true/false flag indicating
            %    whether to Replace all whitespace characters in the input
            %    text by spaces after tab expansion. Note that if
            %    expandTabs is false and ReplaceWhiteSpaceChars is true,
            %    every tab will be converted to a single space! Default is
            %    true.
            %
            %  'FixSentenceEndings' - true/false flag indicating whether to
            %    ensure that sentence-ending punctuation is always followed
            %    by two spaces. Off by default because the algorithm is
            %    (unavoidably) imperfect. Default is false.
            %
            %  'BreakLongWords' - true/false flag indicating whether to
            %    break words longer than 'width'. If false, those words
            %    will not be broken, and some lines might be longer than
            %    'width'. Default is true.
            %
            %  'DropWhiteSpaceChars' - true/false flag indicating whether
            %    to drop leading and trailing whitespace from lines.
            %    Default is true.
            %
            %  'BreakOnHyphens' - true/false flag indicating whether to
            %    allow breaking hyphenated words. If true, wrapping will
            %    occur preferably on whitespaces and right after hyphens
            %    part of compound words. Default is true.
            %
            % Output
            %
            %  twobj - TextWrapper object
            %
            %
            % See Also: TextWrapper.wrap, TextWrapper.fill, 
            %           TextWrapper.wraplines, TextWrapper.wraptext 
            %
            
            % The public interface of TextWrapper consists of the wrap()
            % and fill() methods; the other methods are just there for
            % subclasses to override in order to tweak the default
            % behaviour. If you want to completely replace the main
            % wrapping algorithm, you'll probably have to override
            % _wrap_chunks().
                     
            options.Width = 70;
            options.InitialIndent = '';
            options.SubsequentIndent = '';
            options.ExpandTabs = true;
            options.ReplaceWhiteSpaceChars = true;
            options.FixSentenceEndings = false;
            options.BreakLongWords = true;
            options.DropWhiteSpaceChars = true;
            options.BreakOnHyphens = true;
            
            options = parse_pv_pairs (options, varargin);
            
            assert (isscalar (options.Width) && (options.Width > 0) && (options.Width >= 1), ...
                'Width must be a scalar value >= 1' );
            
            assert (ischar (options.InitialIndent), ...
                'InitialIndent must be character vector' );
            assert (ischar (options.SubsequentIndent), ...
                'SubsequentIndent must be character vector' );
            
            assert (isscalar (options.ExpandTabs) && islogical (options.ExpandTabs), ...
                'ExpandTabs must be a scalar logical value (true/false)' );
            assert (isscalar (options.ReplaceWhiteSpaceChars) && islogical (options.ReplaceWhiteSpaceChars), ...
                'ReplaceWhiteSpaceChars must be a scalar logical value (true/false)' );
            assert (isscalar (options.FixSentenceEndings) && islogical (options.FixSentenceEndings), ...
                'FixSentenceEndings must be a scalar logical value (true/false)' );
            assert (isscalar (options.BreakLongWords) && islogical (options.BreakLongWords), ...
                'BreakLongWords must be a scalar logical value (true/false)' );
            assert (isscalar (options.DropWhiteSpaceChars) && islogical (options.DropWhiteSpaceChars), ...
                'DropWhiteSpaceChars must be a scalar logical value (true/false)' );
            assert (isscalar (options.BreakOnHyphens) && islogical (options.BreakOnHyphens), ...
                'BreakOnHyphens must be a scalar logical value (true/false)' );
            
%             
%             whitespace_trans = string.maketrans(whitespacechars, ' ' * len(whitespacechars))
% 
%             unicodewhitespacechars_trans = {}
%             uspace = ord(u' ')
%             for x in map(ord, whitespacechars):
%                 unicodewhitespacechars_trans[x] = uspace

            % This funky little regex is just the trick for splitting
            % text up into word-wrappable chunks.  E.g.
            %   "Hello there -- you goof-ball, use the -b option!"
            % splits into
            %   Hello/ /there/ /--/ /you/ /goof-/ball,/ /use/ /the/ /-b/ /option!
            % (after stripping out empty strings).
            self.wordsep_re = [ ...
                '(\s+|', ...                                   % any whitespace
                '[^\s\w]*\w+[^0-9\W]-(?=\w+[^0-9\W])|', ...    % hyphenated words
                '(?<=[\w\!\"\''\&\.\,\?])-{2,}(?=\w))'];       % em-dash

            % This less funky little regex just split on recognized spaces. E.g.
            %   "Hello there -- you goof-ball, use the -b option!"
            % splits into
            %   Hello/ /there/ /--/ /you/ /goof-ball,/ /use/ /the/ /-b/ /option!/
            self.wordsep_simple_re = '(\s+)';

            % XXX this is not locale- or charset-aware -- string.lowercase
            % is US-ASCII only (and therefore English-only)
            self.sentence_end_re = [ '[%s]', ...        % lowercase letter
                                     '[\.\!\?]', ...    % sentence-ending punct.
                                     '[\"\'']?', ...    % optional end-of-quote
                                     '\Z'  ];           % end of chunk
                                         % string.lowercase)
            
                     
            self.width = round (options.Width);
            self.initialIndent = options.InitialIndent;
            self.subsequentIndent = options.SubsequentIndent;
            self.expandTabs = options.ExpandTabs;
            self.ReplaceWhiteSpaceChars = options.ReplaceWhiteSpaceChars;
            self.fixSentenceEndings = options.FixSentenceEndings;
            self.breakLongWords = options.BreakLongWords;
            self.dropWhiteSpaceChars = options.DropWhiteSpaceChars;
            self.breakOnHyphens = options.BreakOnHyphens;

%             % recompile the regexes for Unicode mode -- done in this clumsy way for
%             % backwards compatibility because it's rather common to monkey-patch
%             % the TextWrapper class' wordsep_re attribute.
%             self.wordsep_re_uni = re.compile(self.wordsep_re.pattern, re.U)
%             self.wordsep_simple_re_uni = re.compile(
%                 self.wordsep_simple_re.pattern, re.U)

        end
        
        
        function lines = wrap (self, text)
            % wrap paragraph to desired width and return lines in cell array
            %
            % Syntax
            %
            % lines = wrap(twobj, text)
            %
            % Description
            %
            % Reformat the single paragraph in 'text' so it fits in lines
            % of no more than the number of columns specified in the
            % 'width' property, and return a cell array of wrapped lines.
            % By default Tabs in 'text' are expanded, and all other
            % whitespace characters (including newline) are converted to
            % space.
            %
            % Input
            %
            %  twobj - TextWrapper object
            %
            %  text - character vector containing the paragraph to be
            %    wrapped
            %
            % Output
            %
            %  lines - cell array containg each line of the wrapped text
            %
            %
            % See also: TexWrapper.fill, TextWrapper.wraplines
            %

            if numel (text) < self.width
                
                lines = {text};
                
            else
                text = self.mungewhitespacechars (text);

                chunks = self.split (text);

                if self.fixSentenceEndings
                    chunks = self.adjust_sentence_endings (chunks);
                end

                lines = self.wrap_chunks (chunks);
            end

        end
        

        function text = fill (self, text)
            % wrap paragraph to desired width and return as a character vector
            %
            % Syntax
            %
            % lines = fill(twobj, text)
            %
            % Description
            %
            % Reformat the single paragraph in 'text' so it fits in lines
            % of no more than the number of columns specified in the
            % 'width' property, and return a cell array of wrapped lines.
            % By default Tabs in 'text' are expanded, and all other
            % whitespace characters (including newline) are converted to
            % space.
            %
            % Input
            %
            %  twobj -  TextWrapper object
            %
            %  text - character vector containing the paragraph to be
            %    wrapped
            %
            % Output
            %
            %  text - character vector containing the wrapped paragraph
            %
            %
            % See also: TexWrapper.wrap, TextWrapper.wraptext
            %

            text = strjoin (self.wrap(text), sprintf ('\n')); %#ok<*SPRINTFN>

        end
            
    end

    methods (Access = protected)
        % (possibly useful for subclasses to override)

        function text = mungewhitespacechars(self, text)
            % _mungewhitespacechars(text : string) -> string
            % 
            % Munge whitespace in text: expand tabs and convert all other
            % whitespace characters to spaces.  Eg. " foo\tbar\n\nbaz"
            % becomes " foo    bar  baz".

            if self.expandTabs
                text = strrep (text, sprintf ('\t'), repmat (' ' , 1, 8));
            end
            
            if self.ReplaceWhiteSpaceChars
                for ind = 1:numel (self.whitespacechars)
                    text = strrep (text, self.whitespacechars{ind}, ' ');
                end
            end
            
        end


        function chunks = split(self, text)
            % _split(text : string) -> [string]
            % 
            % Split the text to wrap into indivisible chunks.  Chunks are
            % not quite the same as words; see _wrap_chunks() for full
            % details.  As an example, the text
            %   Look, goof-ball -- use the -b option!
            % breaks into the following chunks:
            %   'Look,', ' ', 'goof-', 'ball', ' ', '--', ' ',
            %   'use', ' ', 'the', ' ', '-b', ' ', 'option!'
            % if breakOnHyphens is True, or in:
            %   'Look,', ' ', 'goof-ball', ' ', '--', ' ',
            %   'use', ' ', 'the', ' ', '-b', ' ', option!'
            % otherwise.

            if self.breakOnHyphens
                pat = self.wordsep_re;
            else
                pat = self.wordsep_simple_re;
            end
            [split, match, tokeninds] = regexp (text, pat, 'split', 'match', 'tokenExtents');
            
            if ~isempty (tokeninds)
                
                if tokeninds{1}(1) == 1
                    if numel (match) == 1 ...
                            && all(cellfun (@isempty, split, 'UniformOutput', true))
                        % edge case handles text = '   ', i./e. just space
                        chunks = {};
                    else
                        if isempty (split{1})
                            split(1) = [];
                        end
                        chunks = cell (1, numel(split) + numel(match));
                        chunks(1:2:end-1) = match;
                        chunks(2:2:end) = split;
                    end
                else
                    chunks = cell (1, numel(split) + numel(match));
                    chunks(1:2:end) = split;
                    chunks(2:2:end) = match;
                end

                % remove empty chunks
                rmind = [];
                for ind = 1:numel (chunks)
                    if isempty (chunks{ind})
                        rmind = [rmind,ind];
                    end
                end
                chunks(rmind) = [];
            else
                chunks = {};
            end
        end

        function chunks = adjust_sentence_endings(self, chunks)
            % _fixSentenceEndings(chunks : [string])
            % 
            % Correct for sentence endings buried in 'chunks'.  Eg. when the
            % original text contains "... foo.\nBar ...", mungewhitespacechars()
            % and split() will convert that to [..., "foo.", " ", "Bar", ...]
            % which has one too few spaces; this method simply changes the one
            % space to two.

            i = 1;
            patsearch = self.sentence_end_re;
            while i < numel(chunks)
                if strcmp (chunks{i+1}, ' ') ...
                        && ~isempty(regexp(chunks{i}, patsearch, 'once'))
                    
                    chunks{i+1} = '  ';
                    i = i + 2;
                    
                else
                    
                    i = i + 1;
                    
                end
            end
            
        end

        function [reversed_chunks, cur_line] = handle_long_word(self, reversed_chunks, cur_line, cur_len, width)
            % _handle_long_word(chunks : [string],
            %                      cur_line : [string],
            %                      cur_len : int, width : int)
            % 
            % Handle a chunk of text (most likely a word, not whitespace) that
            % is too long to fit in any line.

            % Figure out when indent is larger than the specified width, and make
            % sure at least one character is stripped off on every pass
            if width < 1
                space_left = 1;
            else
                space_left = width - cur_len;
            end

            if self.breakLongWords
                % If we're allowed to break long words, then do so: put as much
                % of the next chunk onto the current line as will fit.
            
%                 cur_line.append(reversed_chunks{end}(1:space_left))
                cur_line = [cur_line, {reversed_chunks{end}(1:space_left)}];

                reversed_chunks{end} = reversed_chunks{end}(space_left:end);
            
            elseif isempty(cur_line)
                % Otherwise, we have to preserve the long word intact.  Only add
                % it to the current line if there's nothing already there --
                % that minimizes how much we violate the width constraint.
%                 cur_line.append(reversed_chunks.pop())
                cur_line = [cur_line, reversed_chunks(end)];
                reversed_chunks(end) = [];
                
            end

            % If we're not allowed to break long words, and there's already
            % text on the current line, do nothing.  Next time through the
            % main loop of _wrap_chunks(), we'll wind up here again, but
            % cur_len will be zero, so the next line will be entirely
            % devoted to the long word that we can't handle right now.
            
        end

        function lines = wrap_chunks(self, chunks)
            % _wrap_chunks(chunks : [string]) -> [string]
            % 
            % Wrap a sequence of text chunks and return a list of lines of
            % length 'self.width' or less.  (If 'breakLongWords' is false,
            % some lines may be longer than this.)  Chunks correspond roughly
            % to words and the whitespace between them: each chunk is
            % indivisible (modulo 'breakLongWords'), but a line break can
            % come between any two chunks.  Chunks should not have internal
            % whitespace; ie. a chunk is either all whitespace or a "word".
            % Whitespace chunks will be removed from the beginning and end of
            % lines, but apart from that whitespace is preserved.

            lines = {};
            if self.width <= 0
                error ('invalid width %d (must be > 0)', self.width);
            end
            % Arrange in reverse order so items can be efficiently popped
            % from a stack of chucks.
            chunks = fliplr (chunks);

            while ~isempty(chunks)

                % Start the list of chunks that will make up the current line.
                % cur_len is just the length of all the chunks in cur_line.
                cur_line = {};
                cur_len = 0;

                % Figure out which static string will prefix this line.
                if ~isempty(lines)
                    indent = self.subsequentIndent;
                else
                    indent = self.initialIndent;
                end
                % Maximum width for this line.
                thislinewidth = self.width - numel(indent);

                % First chunk on line is whitespace -- drop it, unless this
                % is the very beginning of the text (ie. no lines started yet).
                if self.dropWhiteSpaceChars ...
                        && isempty(strtrim (chunks{end})) ...
                        && ~isempty(lines)
                    
                    chunks(end) = [];
                    
                end
                
                while ~isempty (chunks)
                    
                    l = numel(chunks{end});

                    if cur_len + l <= thislinewidth
                        % Can at least squeeze this chunk onto the current line.
                        cur_line(end+1) = chunks(end);
                        chunks(end) = [];
                        
                        cur_len = cur_len + l;
                    
                    else
                        % Nope, this line is full.
                        break;
                    end
                end
                
                % The current line is full, and the next chunk is too big to
                % fit on *any* line (not just this one).
                if ~isempty(chunks) && (numel(chunks{end}) > thislinewidth)
                    
                    [chunks, cur_line] = self.handle_long_word(chunks, cur_line, cur_len, thislinewidth);
                    
                end
                
                % If the last chunk on this line is all whitespace, drop it.
                if self.dropWhiteSpaceChars ...
                        && ~isempty(cur_line) ...
                        && isempty(strtrim(cur_line{end}))
                    
                    cur_line(end) = [];
                    
                end
                
                % Convert current line back to a string and store it in list
                % of all lines (return value).
                if ~isempty(cur_line)
                    lines = [lines, [indent, strjoin(cur_line, '')]];
                end
                
            end
            
        end
        
    end
    
    methods (Static)
        
        function wrappedlines = wraplines (text, varargin)
            % wrap paragraph to desired width and returns lines in cell array
            %
            % Syntax
            %
            % wrappedlines = TextWrapper.wraplines (text)
            % wrappedlines = TextWrapper.wraplines (..., 'Parameter', Value)
            %
            % Description
            %
            % Reformat the single paragraph in 'text' so it fits in lines
            % of no more than the number of columns specified in the
            % 'width' property, and return a cell array of wrapped lines.
            % By default Tabs in 'text' are expanded, and all other
            % whitespace characters (including newline) are converted to
            % space.
            %
            % This is a static method of the TextWrapper class provided as
            % a convenience function.
            %
            % Input
            %
            %  text - character vector containing the paragraph to be
            %    wrapped
            %
            % Additional arguments may be supplied as parameter-value
            % pairs. The available options are:
            %
            %  'Width' - the maximum width of wrapped lines (unless
            %    breakLongWords is false). Default is 70.
            %
            %  'InitialIndent' - character vector that will be prepended to
            %    the first line of wrapped output; Counts towards each
            %    line's width. Default is ''.
            %
            %  'SubsequentIndent' - character vector that will be prepended to
            %    the all lines of wrapped output except the first line.
            %    Counts towards each line's width. Default is ''.
            %
            %  'ExpandTabs' - true/false flag indicating whether to expand 
            %    tabs in input text to spaces before further processing.
            %    Each tab will become 1 .. 8 spaces, depending on its
            %    position in its line. If false, each tab is treated as a
            %    single character. Default is true.
            %
            %  'ReplaceWhiteSpaceChars' - true/false flag indicating
            %    whether to Replace all whitespace characters in the input
            %    text by spaces after tab expansion. Note that if
            %    expandTabs is false and ReplaceWhiteSpaceChars is true,
            %    every tab will be converted to a single space! Default is
            %    true.
            %
            %  'FixSentenceEndings' - true/false flag indicating whether to
            %    ensure that sentence-ending punctuation is always followed
            %    by two spaces. Off by default because the algorithm is
            %    (unavoidably) imperfect. Default is false.
            %
            %  'BreakLongWords' - true/false flag indicating whether to
            %    break words longer than 'width'. If false, those words
            %    will not be broken, and some lines might be longer than
            %    'width'. Default is true.
            %
            %  'DropWhiteSpaceChars' - true/false flag indicating whether
            %    to drop leading and trailing whitespace from lines.
            %    Default is true.
            %
            %  'BreakOnHyphens' - true/false flag indicating whether to
            %    allow breaking hyphenated words. If true, wrapping will
            %    occur preferably on whitespaces and right after hyphens
            %    part of compound words. Default is true.
            %
            % Output
            %
            %  lines - cell array containg each line of the wrapped text
            %
            %
            % See also: TexWrapper.wraptext, TextWrapper
            %
            
            w = TextWrapper(varargin{:});
            
            wrappedlines = w.wrap (text);
            
        end
        
        function text = wraptext (text, varargin)
            % wrap paragraph to desired width and returns as a character vector
            %
            % Syntax
            %
            % text = TextWrapper.wraptext (text)
            % text = TextWrapper.wraptext (..., 'Parameter', Value)
            %
            % Description
            %
            % Reformat the single paragraph in 'text' so it fits in lines
            % of no more than the number of columns specified in the
            % 'width' property, and return a character vector containing
            % the wrapped text. By default Tabs in 'text' are expanded, and
            % all other whitespace characters (including newline) are
            % converted to space.
            %
            % This is a static method of the TextWrapper class provided as
            % a convenience function.
            %
            % Input
            %
            %  text - character vector containing the paragraph to be
            %    wrapped
            %
            % Additional arguments may be supplied as parameter-value
            % pairs. The available options are:
            %
            %  'Width' - the maximum width of wrapped lines (unless
            %    breakLongWords is false). Default is 70.
            %
            %  'InitialIndent' - character vector that will be prepended to
            %    the first line of wrapped output; Counts towards each
            %    line's width. Default is ''.
            %
            %  'SubsequentIndent' - character vector that will be prepended to
            %    the all lines of wrapped output except the first line.
            %    Counts towards each line's width. Default is ''.
            %
            %  'ExpandTabs' - true/false flag indicating whether to expand 
            %    tabs in input text to spaces before further processing.
            %    Each tab will become 1 .. 8 spaces, depending on its
            %    position in its line. If false, each tab is treated as a
            %    single character. Default is true.
            %
            %  'ReplaceWhiteSpaceChars' - true/false flag indicating
            %    whether to Replace all whitespace characters in the input
            %    text by spaces after tab expansion. Note that if
            %    expandTabs is false and ReplaceWhiteSpaceChars is true,
            %    every tab will be converted to a single space! Default is
            %    true.
            %
            %  'FixSentenceEndings' - true/false flag indicating whether to
            %    ensure that sentence-ending punctuation is always followed
            %    by two spaces. Off by default because the algorithm is
            %    (unavoidably) imperfect. Default is false.
            %
            %  'BreakLongWords' - true/false flag indicating whether to
            %    break words longer than 'width'. If false, those words
            %    will not be broken, and some lines might be longer than
            %    'width'. Default is true.
            %
            %  'DropWhiteSpaceChars' - true/false flag indicating whether
            %    to drop leading and trailing whitespace from lines.
            %    Default is true.
            %
            %  'BreakOnHyphens' - true/false flag indicating whether to
            %    allow breaking hyphenated words. If true, wrapping will
            %    occur preferably on whitespaces and right after hyphens
            %    part of compound words. Default is true.
            %
            % Output
            %
            %  lines - cell array containg each line of the wrapped text
            %
            %
            % See also: TexWrapper.wraplines, TextWrapper
            %
            
            w = TextWrapper(varargin{:});
            
            text = w.fill (text);
            
        end
        
    end
    
end

% % -- Convenience interface ---------------------------------------------
% 
% def wrap(text, width=70, **kwargs):
%     """Wrap a single paragraph of text, returning a list of wrapped lines.
% 
%     Reformat the single paragraph in 'text' so it fits in lines of no
%     more than 'width' columns, and return a list of wrapped lines.  By
%     default, tabs in 'text' are expanded with string.expandtabs(), and
%     all other whitespace characters (including newline) are converted to
%     space.  See TextWrapper class for available keyword args to customize
%     wrapping behaviour.
%     """
%     w = TextWrapper(width=width, **kwargs)
%     return w.wrap(text)
% 
% def fill(text, width=70, **kwargs):
%     """Fill a single paragraph of text, returning a new string.
% 
%     Reformat the single paragraph in 'text' to fit in lines of no more
%     than 'width' columns, and return a new string containing the entire
%     wrapped paragraph.  As with wrap(), tabs are expanded and other
%     whitespace characters converted to space.  See TextWrapper class for
%     available keyword args to customize wrapping behaviour.
%     """
%     w = TextWrapper(width=width, **kwargs)
%     return w.fill(text)
% 
% 
% % -- Loosely related functionality -------------------------------------
% 
% whitespacechars_only_re = re.compile('^[ \t]+$', re.MULTILINE)
% _leadingwhitespacechars_re = re.compile('(^[ \t]*)(?:[^ \t\n])', re.MULTILINE)
% 
% function dedent(text):
%     % Remove any common leading whitespace from every line in `text`.
%     % 
%     % This can be used to make triple-quoted strings line up with the left
%     % edge of the display, while still presenting them in the source code
%     % in indented form.
%     % 
%     % Note that tabs and spaces are both treated as whitespace, but they
%     % are not equal: the lines "  hello" and "\thello" are
%     % considered to have no common leading whitespace.  (This behaviour is
%     % new in Python 2.5; older versions of this module incorrectly
%     % expanded tabs before searching for common leading whitespace.)
%     
%     % Look for the longest leading string of spaces and tabs common to
%     % all lines.
%     margin = None
%     text = whitespacechars_only_re.sub('', text)
%     indents = _leadingwhitespacechars_re.findall(text)
%     for indent in indents:
%         if margin is None:
%             margin = indent
% 
%         % Current line more deeply indented than previous winner:
%         % no change (previous winner is still on top).
%         elif indent.startswith(margin):
%             pass
% 
%         % Current line consistent with and no deeper than previous winner:
%         % it's the new winner.
%         elif margin.startswith(indent):
%             margin = indent
% 
%         % Current line and previous winner have no common whitespace:
%         % there is no margin.
%         else:
%             margin = ""
%             break
% 
%     % sanity check (testing/debugging only)
%     if 0 and margin:
%         for line in text.split("\n"):
%             assert not line or line.startswith(margin), \
%                    "line = %r, margin = %r" % (line, margin)
% 
%     if margin:
%         text = re.sub(r'(?m)^' + margin, '', text)
%     return text
