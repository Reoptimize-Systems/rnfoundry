classdef hydroBodyMexForces < wsim.hydroBody & cppinterface
% represents one hydrodynamically interacting body in a hydrodynamic system
%
% Syntax
%
% hb = hydroBody (filename)
%
% Description
%
% wsim.hydroBody represents one hydrodynamically interacting body in a
% hydrodynamic system. It is intended to be used in conjunction with the
% wsim.hydroSystem class which manages a collection of wsim.hydroBody
% objects and is used to perform the time domain simulation of the bodies.
%
% The body simulation data is contained in a case directory. This case
% directory should contain two subdirectories, 'hydroData' and 'geometry'.
%
% The hydroData subdirectory should contain either one .h5 file containing
% the output of the BEMIO files which process the output of various BEM
% solvers to a format understood by the hydroBody class, or, a collection
% of .mat files containing hydroData structures which can be directly
% loaded by the body. In this case, there will be one mat file for each
% body.
%
% The geometry subdirectory should contain a collection of STL files, one
% for each body. 
%
% wsim.hydroBody Methods:
%
%   hydroBody - constructor for the hydroBody class
%   adjustMassMatrix - Merges diagonal term of added mass matrix to the mass matrix
%   advanceStep - advance to the next time step, accepting the current time
%   bodyGeo - Reads an STL mesh file and calculates areas and centroids
%   checkInputs - Checks the user inputs
%   forceAddedMass - Recomputes the real added mass force time history for the
%   getVelHist - not documented
%   hydroForcePre - performs pre-processing calculations to populate hydroForce structure
%   hydroForces - hydroForces calculates the hydrodynamic forces acting on a
%   hydrostaticForces - calculates the hydrostatic forces acting on the body
%   lagrangeInterp - not documented
%   linearExcitationForces - calculates linear wave excitation forces during transient
%   linearInterp - not documented
%   listInfo - Display some information about the body at the command line
%   loadHydroData - load hydrodynamic data from file or variable
%   makeMBDynComponents - creates mbdyn components for the hydroBody
%   morrisonElementForce - not documented
%   nonlinearExcitationForces - calculates the non-linear excitation forces on the body
%   offsetXYZ - Function to move the position vertices
%   plotStl - Plots the body's mesh and normal vectors
%   radForceODEOutputfcn - OutputFcn to be called after every completed ode time step
%   radForceSSDerivatives - wsim.hydroBody/radForceSSDerivatives is a function.
%   radiationForces - calculates the wave radiation forces
%   readH5File - Reads an HDF5 file containing the hydrodynamic data for the body
%   restoreMassMatrix - Restore the mass and added-mass matrix back to the original value
%   rotateXYZ - Function to rotate a point about an arbitrary axis
%   saveHydroData - saves the body's hydrodata structure to a .mat file
%   setCaseDirectory - set the case directory for the simulation the body is part of
%   setInitDisp - Sets the initial displacement when having initial rotation
%   storeForceAddedMass - Store the modified added mass and total forces history (inputs)
%   timeDomainSimReset - resets the body in readiness for a transient simulation
%   timeDomainSimSetup - sets up the body in preparation for a transient simulation
%   viscousDamping - not documented
%   waveElevation - calculate the wave elevation at centroids of triangulated surface
%   write_paraview_vtp - Writes vtp files for visualization with ParaView
%
%
% See Also: wsim.hydroSystem
%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2014 the National Renewable Energy Laboratory and Sandia Corporation
% Modified 2017 by The University of Edinburgh
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    methods
        
        function obj = hydroBodyMexForces (filename)
            
            obj@wsim.hydroBody (filename);
            
            if isoctave
                mexfcn = str2func ('wsim.mex_hydrobody');
                % work around for Ocave bug #46659
            else
                mexfcn = @wsim.mex_hydrobody;
            end
            
            obj@cppinterface (mexfcn);
            
        end
        
        function timeDomainSimSetup (obj, waves, simu, bodynum)
            % sets up the body in preparation for a transient simulation
            %
            % Syntax
            %
            % timeDomainSimSetup (hb, waves, simu, bodynum)
            %
            % Desciription
            %
            % timeDomainSimSetup initialises various parmaters and settings in
            % preparation for performing a transient simulation based on
            % the ODE solver routines. 
            %
            % Input
            %
            %   hb - hydroBody object
            %
            %   waves - waveClass object with the desired wave parameters
            %     to be used in the simulation
            %
            %   simu - simulationClass Object with the desired simulation
            %     parameters to be used in the simulation
            %
            %   bodynum - the number associated with this body. Typically
            %     this is generated by a parent hydrosys object
            %
            % 
            
            timeDomainSimSetup@wsim.hydroBody (obj, waves, simu, bodynum);
            
            if isempty (obj.radForce_IRKB_interp)
                obj.radForce_IRKB_interp = zeros ([12,601,6]);
                 obj.CIdt = 0;
            end
            
                               
            check.multicheck ( @(x) ~isempty(x), ...
                               'input cannot be empty', ...
                               obj.bodyNumber, ...
                               obj.doViscousDamping, ...
                               obj.doLinearDamping, ...
                               obj.doNonLinearFKExcitation, ...
                               obj.doMorrisonElementViscousDrag, ...
                               obj.excitationMethodNum, ...
                               obj.freeSurfaceMethodNum, ...
                               obj.radiationMethodNum, ...
                               obj.hydroForce.fExt.re, ...
                               obj.hydroForce.fExt.im, ...
                               obj.hydroForce.visDrag, ...
                               obj.hydroForce.linearHydroRestCoef, ...
                               obj.hydroForce.fDamping, ...
                               obj.hydroForce.fAddedMass, ...
                               obj.waves.A, ...
                               obj.waves.w, ...
                               obj.waves.phase, ...
                               obj.waves.dw, ...
                               obj.waves.k, ...
                               ...obj.waves.S, ...
                               obj.linearDamping, ...
                               obj.cg, ...
                               obj.cb, ...
                               obj.CIdt, ...
                               obj.radForce_IRKB_interp(:,:,1), ...
                               obj.radForce_IRKB_interp(:,:,2), ...
                               obj.radForce_IRKB_interp(:,:,3), ...
                               obj.radForce_IRKB_interp(:,:,4), ...
                               obj.radForce_IRKB_interp(:,:,5), ...
                               obj.radForce_IRKB_interp(:,:,6), ...
                               obj.simu.g, ...
                               obj.simu.rho, ...
                               obj.simu.dtFeNonlin, ...
                               obj.simu.dt, ...
                               obj.simu.startTime, ...
                               obj.simu.b2b, ...
                               obj.bodyTotal, ...
                               obj.simu.CTTime, ...
                               obj.dof, ...
                               obj.simu.rampT, ...
                               obj.hydroForce.storage.mass, ...
                               obj.dispVol ...
                             );
                   
            if isvector (obj.linearDamping)
                linearDamping = diag (obj.linearDamping);
            else
                linearDamping = obj.linearDamping;
            end
            
            assert ( isnumeric (linearDamping) && isreal (linearDamping) && size (linearDamping, 1) == 6 && size (linearDamping, 2) == 6, ...
                     'linearDamping should be a 6x6 real numeric matrix' );
            
                 
            fprintf (1, 'bodyNumber:\n'); disp(obj.bodyNumber);
            fprintf (1, 'doViscousDamping:\n'); disp(obj.doViscousDamping);
            fprintf (1, 'doLinearDamping:\n'); disp(obj.doLinearDamping);
            fprintf (1, 'doNonLinearFKExcitation:\n'); disp(obj.doNonLinearFKExcitation);
            fprintf (1, 'doMorrisonElementViscousDrag:\n'); disp(obj.doMorrisonElementViscousDrag);
            fprintf (1, 'excitationMethodNum:\n'); disp(obj.excitationMethodNum);
            fprintf (1, 'freeSurfaceMethodNum:\n'); disp(obj.freeSurfaceMethodNum);
            fprintf (1, 'radiationMethodNum:\n'); disp(obj.radiationMethodNum);
            fprintf (1, 'obj.hydroForce.fExt.re:\n'); disp(obj.hydroForce.fExt.re);
            fprintf (1, 'obj.hydroForce.fExt.im:\n'); disp(obj.hydroForce.fExt.im);
            fprintf (1, 'hydroForce.visDrag:\n'); disp(obj.hydroForce.visDrag);
            fprintf (1, 'hydroForce.linearHydroRestCoef:\n'); disp(obj.hydroForce.linearHydroRestCoef);
            fprintf (1, 'hydroForce.fDamping:\n'); disp(obj.hydroForce.fDamping);
            fprintf (1, 'hydroForce.fAddedMass:\n'); disp(obj.hydroForce.fAddedMass);
            fprintf (1, 'waves.A:\n'); disp(obj.waves.A);
            fprintf (1, 'waves.w:\n'); disp(obj.waves.w);
            fprintf (1, 'waves.phase:\n'); disp(obj.waves.phase);
            fprintf (1, 'waves.dw:\n'); disp(obj.waves.dw);
            fprintf (1, 'waves.k:\n'); disp(obj.waves.k);
            fprintf (1, 'waves.S:\n'); disp(obj.waves.S);
%             fprintf (1, 'linearDamping:\n'); disp(obj.linearDamping);
            fprintf (1, 'linearDamping:\n'); disp(linearDamping);
            fprintf (1, 'cg:\n'); disp(obj.cg);
            fprintf (1, 'cb:\n'); disp(obj.cb);
            fprintf (1, 'CIdt:\n'); disp(obj.CIdt);
%             fprintf (1, 'radForce_IRKB_interp(:,:,1):\n'); disp(obj.radForce_IRKB_interp(:,:,1));
%             fprintf (1, 'radForce_IRKB_interp(:,:,2):\n'); disp(obj.radForce_IRKB_interp(:,:,2));
%             fprintf (1, 'radForce_IRKB_interp(:,:,3):\n'); disp(obj.radForce_IRKB_interp(:,:,3));
%             fprintf (1, 'radForce_IRKB_interp(:,:,4):\n'); disp(obj.radForce_IRKB_interp(:,:,4));
%             fprintf (1, 'radForce_IRKB_interp(:,:,5):\n'); disp(obj.radForce_IRKB_interp(:,:,5));
%             fprintf (1, 'radForce_IRKB_interp(:,:,6):\n'); disp(obj.radForce_IRKB_interp(:,:,6));
            fprintf (1, 'simu.g:\n'); disp(obj.simu.g);
            fprintf (1, 'simu.rho:\n'); disp(obj.simu.rho);
            fprintf (1, 'simu.dtFeNonlin:\n'); disp(obj.simu.dtFeNonlin);
            fprintf (1, 'simu.dt:\n'); disp(obj.simu.dt);
            fprintf (1, 'simu.startTime:\n'); disp(obj.simu.startTime);
            fprintf (1, 'simu.b2b:\n'); disp(obj.simu.b2b);
            fprintf (1, 'bodyTotal:\n'); disp(obj.bodyTotal);
%             fprintf (1, 'simu.CTTime:\n'); disp(obj.simu.CTTime);
            fprintf (1, 'dof:\n'); disp(obj.dof);
            fprintf (1, 'simu.rampT:\n'); disp(obj.simu.rampT);
            fprintf (1, 'hydroForce.storage.mass:\n'); disp(obj.hydroForce.storage.mass);
            fprintf (1, 'dispVol:\n'); disp(obj.dispVol);
            
            
            obj.cppcall ( 'Initialize', ...
                          obj.bodyNumber, ...
                          obj.doViscousDamping, ...
                          obj.doLinearDamping, ...
                          obj.doNonLinearFKExcitation, ...
                          obj.doMorrisonElementViscousDrag, ...
                          obj.excitationMethodNum, ...
                          obj.freeSurfaceMethodNum, ...
                          obj.radiationMethodNum, ...
                          obj.hydroForce.fExt.re, ...
                          obj.hydroForce.fExt.im, ...
                          obj.hydroForce.visDrag, ...
                          obj.hydroForce.linearHydroRestCoef, ...
                          obj.hydroForce.fDamping, ...
                          obj.hydroForce.fAddedMass, ...
                          obj.waves.A, ...
                          obj.waves.w, ...
                          obj.waves.phase, ...
                          obj.waves.dw, ...
                          obj.waves.k, ...
                          obj.waves.S, ...
                          linearDamping, ...
                          obj.cg, ...
                          obj.cb, ...
                          obj.CIdt, ...
                          obj.radForce_IRKB_interp(:,:,1), ...
                          obj.radForce_IRKB_interp(:,:,2), ...
                          obj.radForce_IRKB_interp(:,:,3), ...
                          obj.radForce_IRKB_interp(:,:,4), ...
                          obj.radForce_IRKB_interp(:,:,5), ...
                          obj.radForce_IRKB_interp(:,:,6), ...
                          obj.simu.g, ...
                          obj.simu.rho, ...
                          obj.simu.dtFeNonlin, ...
                          obj.simu.dt, ...
                          obj.simu.startTime, ...
                          obj.simu.b2b, ...
                          obj.bodyTotal, ...
                          obj.simu.CTTime, ...
                          obj.dof, ...
                          obj.simu.rampT, ...
                          obj.hydroForce.storage.mass, ...
                          obj.dispVol );

        end

        
        function advanceStep (obj, t, vel, accel)
            % advance to the next time step, accepting the current time
            % step and data into stored solution histories
            %
            % Syntax
            %
            % advanceStep (hb, t, vel, accel)
            %
            % Description
            %
            % advanceStep must be called at the end of each integration
            % time step to update the history of the solution. This is
            % required for the radiation forces when using either the
            % convolution integral method, or the state-space integration
            % using the default internal integration provided by this class
            % (invoked with simu.ssCalc = 1).
            %
            % Input
            %
            %  hb - hydroBody object
            %
            %  t - the last computed time step (the step at which the vel
            %    and accel values are being provided).
            %
            %  vel - vector of velocities and angular velocities, of size
            %    (6 x 1) if there is no body to body interaction (just the
            %    velocities of this body), or size (6*nbodies x 1) if there
            %    is body to body interaction (the velocities of all bodies
            %    in the system).
            %
            %  accel - vector of accelerations and angular accelerations,
            %    of size (6 x 1) if there is no body to body interaction
            %    (just the velocities of this body), or size (6*nbodies x
            %    1) if there is body to body interaction (the accelerations
            %    of all bodies in the system).
            %
            % 
            
%             fprintf (1, 't:\n'); disp(t);
%             fprintf (1, 'vel:\n'); disp(vel);
%             fprintf (1, 'accel:\n'); disp(accel);
            
            advanceStep@wsim.hydroBody (obj, t, vel, accel);
            
            obj.cppcall ('advanceStep', t, vel, accel);
            
        end
        
        function [forces, breakdown] = hydroForces (obj, t, pos, vel, accel)
            
            [forces, breakdown] = obj.cppcall ('hydroForces', t, pos, vel, accel);
            
        end
        
    end
end
