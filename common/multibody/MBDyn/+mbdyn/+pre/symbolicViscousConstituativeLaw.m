classdef symbolicViscousConstituativeLaw < mbdyn.pre.constituativeLaw
% mbdyn.pre.symbolicViscousConstituativeLaw class
%
% Syntax
%
% symvlaw = symbolicViscousConstituativeLaw (scalar_function)
%
% Description
%
% Applies damping according to a symbolic expression defined by the user,
% with the values of the dimensions of the law as variables in the
% expression.
%
% NB: Requires that MBDyn has been compiled with the GiNAC library.
%
% Example 1
%
% % For a 1D constituative law
% 
% law = mbdyn.pre.symbolicViscousConstituativeLaw ('eps', '1000.*eps + 5.*eps^3');
% 
% Example 2
% 
% % For a 3D constituative law
%
% law = mbdyn.pre.symbolicViscousConstituativeLaw ( {'eps1', 'eps3', 'eps3'}, ...
%                                                   {'1000.*eps1 + 5.*eps1^3 - 10.*eps2*eps3', ...
%                                                    '1000.*eps2 + 5.*eps2^3 - 10.*eps3*eps1', ...
%                                                    '1000.*eps3 + 5.*eps3^3 - 10.*eps1*eps2'} );
%       
% mbdyn.pre.symbolicViscousConstituativeLaw Methods:
%
%   symbolicViscousConstituativeLaw - mbdyn.pre.symbolicViscousConstituativeLaw 
%     constructor
%   generateMBDynInputString - generates MBDyn input string for the 
%     symbolicViscousConstituativeLaw
%

    properties
        
        epsilon_prime;
        expression;
        
    end
    
    methods
        
        function self = symbolicViscousConstituativeLaw (epsilon_prime, expression)
            % mbdyn.pre.symbolicViscousConstituativeLaw constructor
            %
            % Syntax
            %
            % symvlaw = symbolicViscousConstituativeLaw (scalar_function)
            %
            % Description
            %
            % Applies damping according to a symbolic expression defined by
            % the user, with the values of the dimensions of the law as
            % variables in the expression.
            %
            % NB: Requires that MBDyn has been compiled with the GiNAC
            % library.
            %
            % Input
            %
            %  epsilon_prime - either a string, a character vector, a cell
            %   string array or array of strings. This contains the names
            %   of the variables in the expression which correspond to the
            %   differrent dimensions of the constituitive law, and which
            %   are used in the expression.
            %
            %  expression - either a string, a character vector, a cell
            %   string array or array of strings. This contains the
            %   symbolic expression, which uses the variables in
            %   epsilon_prime, which correspond to the value of each
            %   dimension of the constituitive law.
            %
            % Output
            %
            %  symvlaw - mbdyn.pre.symbolicViscousConstituativeLaw
            %   object
            %
            % Example 1
            %
            % % For a 1D constituative law
            % 
            % law = mbdyn.pre.symbolicViscousConstituativeLaw ('eps', '1000.*eps + 5.*eps^3');
            % 
            % Example 2
            % 
            % % For a 3D constituative law
            %
            % law = mbdyn.pre.symbolicViscousConstituativeLaw ( {'eps1', 'eps3', 'eps3'}, ...
            %                                                   {'1000.*eps1 + 5.*eps1^3 - 10.*eps2*eps3', ...
            %                                                    '1000.*eps2 + 5.*eps2^3 - 10.*eps3*eps1', ...
            %                                                    '1000.*eps3 + 5.*eps3^3 - 10.*eps1*eps2'} );
            %
            % See Also:
            %
            
            assert (ischar (expression) || iscellstr (epsilon_prime) || isstring (epsilon_prime), ...
                    'epsilon_prime must be a character vector, a cell array of character vectors or a string array' );
            assert (ischar (expression) || iscellstr (expression) || isstring (expression), ...
                    'expression must be a character vector, a cell array of character vectors or a string array' );
                
            % convert whatever was supplied into a cell array of strings
            epsilon_prime = cellstr (epsilon_prime);
            expression = cellstr (expression);
            
            assert (numel (epsilon_prime) == numel (expression), ...
                    'the number of elements in epsilon_prime must be the same as the number of elements in expression' );
            
            
            
            self.type = 'symbolic viscous';
            self.epsilon_prime = epsilon_prime;
            self.expression = expression;
            
        end
        
        function str = generateMBDynInputString (self)
            % generates MBDyn input string for the symbolicViscousConstituativeLaw
            % 
            % Syntax
            %  
            % str = generateMBDynInputString (symvlaw)
            %  
            % Description
            %  
            % generateMBDynInputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %  
            % Input
            %  
            %  symvlaw - mbdyn.pre.symbolicViscousConstituativeLaw object
            %  
            % Output
            %  
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %
            
            str = sprintf ('%s,', self.type);
            
            epsilon_prime_str = '';
            for ind = 1:numel (self.epsilon_prime)
                epsilon_prime_str = [epsilon_prime_str, '"', self.epsilon_prime{ind}, '", '];
            end
            
            str = self.addOutputLine (str, self.commaSepList ('epsilon prime', epsilon_prime_str), 1, false);
            
            str = self.addOutputLine (str, 'expression', 1, true);
            
            for ind = 1:numel (self.expression)
                str = self.addOutputLine (str, ['"', self.expression{ind}, '"'], 2, ind ~= numel (self.expression));
            end

        end
        
    end
    
end