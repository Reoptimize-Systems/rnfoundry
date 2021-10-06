function idot = multiphasecurrentderiv(I, E, R, M)
% calculates the derivatives of the currents in a multiphase RL circuit
% with self and mutual inductances
%
% Syntax
% 
% idot = multiphasecurrentderiv(I, E, R, M)
%
% Description
%
% Solves the right hand side of the differential equation describing the
% circuit below:
%
%                    Vr        L dI/dt
%                  <------      <-----
%                   _____                                                   
%            ______|     |______      _____                                 
%           |      |_____|      UUUUUU     |                                
%           |                              |                                
%           |         R           L        |                                
%           |                              |                                
%         /   \                            |                                
%     E  |  ~  |                           |                                
%         \   /                           \|/ I                               
%           |                              |                                
%           |                              |                                
%           |                              |                                
%           |______________________________|                                
%                             |                                             
%                            _|_                                            
%                           \ \ \                                           
%
% However the circuit may be single or multiphase, with mutual inductances
% between Phases. 
%
% Input
%
% I is a vector of values of the currentin the the machine Phases at the
%   current time.
%
% E is a vector of values of the emfs in each phase at the current time
%
% R is an (n x n) matrix with diagonals the phase resistances of the
%   machine, i.e.
%
%   R = [ Ra  0   0;
%         0   Rb  0;
%         0   0   Rc ];
%
% M is an inductance matrix with diagonal terms the main inductances of the
%   Phases, and off-diagonal terms the mutual inductances between Phases
%
%   M = [ La   Lba  Lca;
%         Lba  Lb   Lcb;
%         Lba  Lcb  Lc  ];
% 
%   Often it is assumed that: La = Lb = Lc = Lp, and Lba = Lca = Lcb = Mpp
%   such that
%
%   M = [ Lp   Mpp  Mpp;
%         Mpp  Lp   Mpp;
%         Mpp  Mpp  Lp ];
%
% For an explaination, see Pillay, P., and R. Krishnan. 1988. 'Modeling of
% permanent magnet motor drives.', IEEE Transactions on Industrial
% Electronics 35(4): 537-541.
% http://ieeexplore.ieee.org/xpl/articleDetails.jsp?tp=&arnumber=9176&
% contentType=Journals+%26+Magazines&searchField%3DSearch_All%26queryText%
% 3DModeling+of+Permanent+Magnet+Motor+Drives
% (Accessed May 3, 2012).

    % a circuit consisting of a resistance and series inductance can be
    % described using the following equation:
    %                                                              
    %
    % E = Vr + L dI/dt 
    % E = RI + L dI/dt
    % E - RI = L dI/dt
    % (E - RI) / L = dI/dt
    %
    % dI/dt = (E - RI) / L
    % 
    % As we have several circuits we are solving at once this is actually
    % division of a matrix by a matrix, or equivalent to 
    %
    % dI/dt = inv(L) * (E - RI)
    %
    % 
    % Calculate di/dt
    idot = ((E(:) - R * I(:))' / M)';
    
%     % Fancier version to take position dependant inductance
%     ydot = ((E - design.R * (y))' / feval(simparams.Lfcn, x))';


end