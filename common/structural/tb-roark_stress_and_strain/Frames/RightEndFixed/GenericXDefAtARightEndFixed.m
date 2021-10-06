function DelA = GenericXDefAtARightEndFixed(CHH, CHV, CHM, HA, VA, MA, LFH)
% GenericXDefAtARightEndFixed: calculates the horizontal deflection at point
% A for a frame with it's right side fixed at the base for the cases 5 to
% 12 in Table 4 of Roark's Formulas for Stress and Strain
%

    DelA = CHH.*HA + CHV.*VA + CHM.*MA - LFH;

end