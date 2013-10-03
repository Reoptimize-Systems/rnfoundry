function DelA = GenericYDefAtARightEndFixed(CHV, CVV, CVM, HA, VA, MA, LFV)
% GenericYDefAtARightEndFixed: calculates the vertical deflection at point
% A for a frame with it's right side fixed at the base for the cases 5 to
% 12 in Table 4 of Roark's Formulas for Stress and Strain
%

    DelA = CHV.*HA + CVV.*VA + CVM.*MA - LFV;

end