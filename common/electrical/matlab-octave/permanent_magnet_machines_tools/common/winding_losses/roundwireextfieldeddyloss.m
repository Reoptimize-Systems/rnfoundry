function P = roundwireextfieldeddyloss(dc, length, rho, dBdt)
% calculates the instantaneous power loss due to eddy currents induced in a
% winding made up of cylindrical conductors in a magnetic field which is
% changing in time but uniform spacially over the conductor cross-section
%
% Syntax
%
% P = roundwireextfieldeddyloss(dc, length, rho, dBdt)
%
% Input
%
%   dc - diameter of the wire cross-section
%
%   length - length of the wire in the magnetic field
%
%   rho - resistivity of the wire
%
%   dBdt - derivative w.r.t. time of the time-changing magnetic field,
%     asumed to be constnt throughout the wire cross-section
%
% Description
%
% 
%                         |                                        
%                         |                                        
%                         |                                        
%                         |   dc                                     
%             <----------------------->                                
%                         |                                        
%                         |                                        
%                     ... | ...                                    
%                  ...    |    ...                                  
%                ..       |       ..                                
%              ..         |         ..                                
%             ..          |          ..                               
%             .           |           .                               
% ------------.-----------------------.-----------                    
%             .           |           .                               
%             .           |           .                               
%              .          |          .                                
%               ..        |        ..                                
%                 ...     |     ...                                  
%                     ....|....                                    
%                         |                                        
%                         |                                        
%                ^ ^ ^ ^ ^|^ ^ ^ ^ ^ ^                             
%                ¦ ¦ ¦ ¦ ¦|¦ ¦ ¦ ¦ ¦ ¦                             
%                ¦ ¦ ¦ ¦ ¦|¦ ¦ ¦ ¦ ¦ ¦                             
%                ¦ ¦ ¦ ¦ ¦|¦ ¦ ¦ ¦ ¦ ¦     dB / dt                        
%                ¦ ¦ ¦ ¦ ¦|¦ ¦ ¦ ¦ ¦ ¦                             
%                ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦                             
%                ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦    
% 
% 
%
% Uses the Squared Field Derivative (SVD) formula taken from:
%
% C. R. Sullivan, "Computationally Efficient Winding Loss Calculation with
% Multiple Windings, Arbitrary Waveforms and Two- or Three-Dimensional
% Field Geometry", IEEE Transactions on Power Electronics, vol. 16, no. 1,
% pp. 142– 150.
%
% This formula is valid for round wire provided dc is not >> skin depth at the
% frequency of applied field.
%

    P = pi .* length .* dc.^4 .* (dBdt).^2 ./ (64 .* rho);

end