function Tmat = tThermal_PMSM(varargin)
% tThermal_PMSM: a function for calculating the steady state temperature of
% the linear PM synchronous machine based on the rotary version in Grauers
% Thesis, Appendix B

    % Thermal constants from Grauers Thesis, we will package these in a
    % structure for convenience
    Tamb = 30;          %Ambient temperature in [oC]
    tConstants.q_vcair = 2.8;      % Volumetric flow rate of cooling air
    tConstants.rho_cair = 1.1;     % Density of cooling air [kg / m^3]
    tConstants.k_thcair = 1010;    % Specific heat capacitivity of cooling air [J / (kg.K)]

    tConstants.alpha(1) = 60;      % heat transfer coeff at stator yoke back [W / (K.m^2)]
    tConstants.alpha(2) = 40;      % heat transfer coeff in air gap [W / (K.m^2)]
    tConstants.alpha(3) = 25;      % heat transfer coeff at end shields [W / (K.m^2)]
    tConstants.alpha(4) = 25;      % heat transfer coeff at end windings [W / (K.m^2)]
    tConstants.alpha(5) = 25;      % heat transfer coeff at rotor yoke back [W / (K.m^2)]

    tConstants.lambda(1) = 38;     % lamda_Fe, thermal conductivity of iron [W / (K.m)]
    tConstants.lambda(2) = 1.8;    % lamda_Coil, thermal conductivity through the coil [W / (K.m)]
    tConstants.lambda(3) = 400;    % lamda_Cu, thermal conductivity of copper, along the coil [W / (K.m)]
    tConstants.lambda(4) = 0.2;    % lamda_i, thermal conductivity of insulation [W / (K.m)]
    tConstants.lambda(5) = 9;      % lamda_m, thermal conductivity of NdFeB [W / (K.m)]
    tConstants.lambda(6) = 0.7;    % lamda_glue, thermal conductivity of magnet glue [W / (K.m)]
    tConstants.lambda(7) = 0.2;    % lamda_GRP, thermal conductivity of GRP magnet protection [W / (K.m)]

    % Constants for Stator yoke hysteresis and eddy current losses
    lossConstants.k_Hysy = 2;         % Empirical hysteresis loss factor for the stator yoke
    lossConstants.k_Ecsy = 1.8;       % Empirical eddy current loss factor for the stator yoke (kftys in Grauer)
    lossConstants.p_Hy = 2.04;        % Hysteresis loss density at 50 Hz and 1.5 T
    lossConstants.p_Ec = 0.76;        % Eddy current loss density at 50 Hz and 1.5 T
    m_Fesy = 467;                     % Stator yoke weight [kg]

    % Constants for Stator teeth hysteresis and eddy current losses
    lossConstants.k_Hyt = 1.2;        % Empirical hysteresis loss factor for the teeth
    lossConstants.k_Ect = 2.5;        % Empirical eddy current loss factor for the teeth (kftd in Grauer)
    m_Fet = 888;                      % Translator yoke weight [kg]

    % Constants fo calculating copper loss 
    lossConstants.alpha_Cu = 3.9e-3;  % Coefficient resistivity for copper, i.e. increase in resistivity per Kelvin [K^-1]
    lossConstants.rho_Cu_25 = 1.7e-8; % Resistivity of copper at 25 degrees celcius [ohm.m]

    % Simulation parameters
    Q = 1;              % Number of parallel models, use to simulate entire pole, or machine
    pp = 99;            % Number of Pole-Pairs
    Nphases = 3;        % Number of Phases
    Vpeak = 2.2;        % Peak machine velocity

    % Physical parameters and dimensions (all in [m] except k_Cu which is unitless)
    k_Cu = 0.8;         % Copper fill factor in the stranded coils excluding coil insulation
    hm0 = 0.1e-3;       % Thickness of magnet glue 
    hm1 = 0.5e-3;       % Thickness of magnet reinforcement 
    hi = 1e-3;          % Coil insulation thickness 
    bt = 11.1e-3;       % Tooth width (bd in Grauer)
    bs = 11.7e-3;       % Slot height
    hsy = 15.9e-3;      % Thickness of stator (toothed part) yoke ((middle bit)
    hry = 15.4e-3;      % Thickness of translator yoke (back iron) 
    ls = 0.55;          % Stator length (Active length of the generator)
    hs = 88e-3;         % Slot height
    hs1 = 1e-3;         % Distance from base of tooth to the point where tooth begins to slope inwards
    hs2 = 4e-3;         % Distance from where tooth begins to slope inwards to point where sloping stops
    Taup = 68.3e-3;     % Pole pitch/width
    bs1 = 3e-3;         % Distance between base of teeth
    bp = 16.2e-3;       % Magnet width
    Bd0 = 1.5;          % Peak teeth flux density
    Bsy = 1.5;          % Peak flux density in the stator yoke     
    hm = 6.3e-3;        % magnet height (lm in my simulations)
    dsi = 2.15;         % Internal diameter of stator 

    % Dynamic parameters
    I = 50;              % Peak current in the winding

    % RMS coil current
    Is = I/sqrt(2); 
    % Copper width in slot (slot width minus insulation thickness on each side)
    bCu = bs - (2*hi);               
    % Total length of end winding
    lb = 2 * Taup * pp * Q;                   
    % RMS fequency at peak velocity
    frms = (Vpeak/(2*sqrt(2)*Taup)); 
    % Slot pitch/width
    Taus = Taup / Nphases;
    % Height of conductor region (two conductor regions per slot)
    hCu = (hs / 2) - 4 * hi;
    % Area of conducting material in the coils
    A_Cus = bCu * hCu * k_Cu;
    % Total winding length for one winding
    l_Cus = 2*(ls + lb); 
    % Width of base of tooth, in the case of the linear PMSM the teeth are the
    % same width at all points, i.e. there is no 'foot' at the end of the teeth
    % so we set this equal to bt 
    Tautb  =  bt;                    

    % Calculate the thermal resistances
    R = tThermalResistances(tConstants.q_vcair, tConstants.rho_cair, tConstants.k_thcair, k_Cu, bt, bs, hsy, hry, ls, dsi, hm0, hm1, hm, hi, hs, hs1, hs2, Taup, bs1, bp, Nphases, pp, Tautb, Q, tConstants.alpha, tConstants.lambda);

    % Hysteresis loss in stator yoke, dep on operational speed; not on T
    P_Hysy = lossConstants.k_Hysy * m_Fesy * lossConstants.p_Hy * (frms / 50) * (Bsy / 1.5)^2; 
    % Eddy current loss in stator yoke, dep on operational speed; not on T
    P_Ecsy = lossConstants.k_Ecsy * m_Fesy * lossConstants.p_Ec * (frms / 50)^2 * (Bsy / 1.5)^2; 
    % Hysteresis loss in stator teeth, dep on operational speed; not on T
    P_Hyt = lossConstants.k_Hyt * m_Fet * lossConstants.p_Hy * (frms / 50) * (Bd0 / 1.5)^2; 
    % Eddy current loss in stator teeth, dep on operational speed; not on T
    P_Ect = lossConstants.k_Ect * m_Fet * lossConstants.p_Ec * (frms / 50)^2 * (Bd0 / 1.5)^2; 

    % Magnet eddy current loss: Dep on size of magnets - how many, split up to
    % stop eddy current paths. Not dep on T. Grauers treats this as a fixed loss,
    % Joekel formula f(v^2)
    p_Ecm = 50; %[W / m^2]
    P_Ecm = 2 * pp * bp * ls * p_Ecm; %Magnet eddy current loss 
    %Could scale with n or om_e

    %Stray losses: At rated load, stray losses assumed to be 20% of core losses
    %at no load, located at stator tooth tip. scale with I^2.
    P_ad_nom = 0.2 * (P_Hysy + P_Ecsy + P_Hyt + P_Ect);
    I_nom = 1; % Nominal current for scaling
    P_ad = P_ad_nom * (I / I_nom)^2; %stray losses scaled

    % Construct conductance matrix, [Gmat], 12by12
    Gmat = tConductanceMatrix_PMSM(R(52:66));

    % Iteratively (i.e. taking account of resistivity temperature dependence)
    % find the steady-state temperature
    [Tmat, Rs] = tSteadyStateTemp_PMSM(Gmat, Tamb, Is, R, lossConstants.rho_Cu_25, lossConstants.alpha_Cu, P_Hyt, P_Hysy, P_Ect, P_Ecsy, P_Ecm, P_ad, ls, lb, bt, bs, Taus, l_Cus, A_Cus, Nphases);

end
