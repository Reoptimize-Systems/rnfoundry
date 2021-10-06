function DelA = GenericAngularRotAtARightEndFixed(CHM, CVM, CMM, HA, VA, MA, LFM)
% GenericAngularRotAtARightEndFixed: calculates the angular rotation at point
% A for a frame with it's right side fixed at the base for the cases 5 to
% 12 in Table 4 of Roark's Formulas for Stress and Strain
%

    DelA = CHM.*HA + CVM.*VA + CMM.*MA - LFM;

end