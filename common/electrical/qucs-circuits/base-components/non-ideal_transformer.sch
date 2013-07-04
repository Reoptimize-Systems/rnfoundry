<Qucs Schematic 0.0.17>
<Properties>
  <View=0,0,974,800,1,0,0>
  <Grid=5,5,1>
  <DataSet=non-ideal_transformer.dat>
  <DataDisplay=non-ideal_transformer.dpl>
  <OpenDisplay=1>
  <Script=non-ideal_transformer.m>
  <RunScript=0>
  <showFrame=0>
  <FrameText0=Title>
  <FrameText1=Drawn By:>
  <FrameText2=Date:>
  <FrameText3=Revision:>
</Properties>
<Symbol>
  <.ID -120 -66 NTR "1=R1_val=1=Primary Winding Resistiance=Resistance" "1=R2_val=1=Secondary Winding Resistance=Resistance" "1=Lm_val=0.1=Primary Leakage Inductance=Inductance" "1=V_ratio=1=Primary to Secondary Voltage Ratio=Unitless" "1=k_val=0.98=Coupling coefficient=">
  <Line -120 -70 80 0 #000080 2 1>
  <Line -120 -150 0 80 #000080 2 1>
  <.PortSym -130 -140 1 0>
  <Line -130 -140 10 0 #000080 2 1>
  <Line -130 -80 10 0 #000080 2 1>
  <.PortSym -130 -80 3 0>
  <EArc -105 -130 0 0 0 5760 #000000 0 1>
  <EArc -75 -120 20 20 1494 2902 #000000 0 1>
  <Line -85 -140 0 60 #000000 0 1>
  <Line -80 -80 0 -60 #000000 0 1>
  <EArc -75 -140 20 20 1494 2902 #000000 0 1>
  <EArc -75 -100 20 20 1494 2902 #000000 0 1>
  <EArc -110 -100 20 20 4374 2902 #000000 0 1>
  <EArc -110 -120 20 20 4374 2902 #000000 0 1>
  <EArc -110 -140 20 20 4374 2902 #000000 0 1>
  <Ellipse -60 -140 10 10 #000000 0 1 #000000 1 1>
  <Ellipse -115 -140 10 10 #000000 0 1 #000000 1 1>
  <.PortSym -30 -140 2 180>
  <.PortSym -30 -80 4 180>
  <Line -40 -140 10 0 #000080 2 1>
  <Line -40 -80 10 0 #000080 2 1>
  <Line -40 -150 0 80 #000080 2 1>
  <Line -120 -150 80 0 #000080 2 1>
</Symbol>
<Components>
  <Port P1 1 70 435 -23 -56 1 0 "1" 1 "analog" 0>
  <Port P2 1 70 595 -23 12 0 0 "3" 1 "analog" 0>
  <R R1 1 185 435 -26 15 0 0 "R1_val" 1 "26.85" 0 "0.0" 0 "0.0" 0 "26.85" 0 "european" 0>
  <Port S2 1 605 595 4 12 1 2 "4" 1 "analog" 0>
  <Port S1 1 605 435 4 -56 0 2 "2" 1 "analog" 0>
  <R R2 1 465 435 -26 15 0 0 "R2_val" 1 "26.85" 0 "0.0" 0 "0.0" 0 "26.85" 0 "european" 0>
  <MUT Tr1 1 340 510 -29 38 0 0 "L1_val" 0 "L2_val" 0 "k_val" 0>
  <Eqn Eqn1 1 720 435 -24 15 0 0 "L1_val=Lm_val" 1 "L2_val=V_ratio^2 * Lm_val" 1 "no" 0>
</Components>
<Wires>
  <70 435 155 435 "" 0 0 0 "">
  <70 595 310 595 "" 0 0 0 "">
  <310 540 310 595 "" 0 0 0 "">
  <370 595 605 595 "" 0 0 0 "">
  <370 540 370 595 "" 0 0 0 "">
  <215 435 310 435 "" 0 0 0 "">
  <310 435 310 480 "" 0 0 0 "">
  <495 435 605 435 "" 0 0 0 "">
  <370 435 370 480 "" 0 0 0 "">
  <370 435 435 435 "" 0 0 0 "">
</Wires>
<Diagrams>
</Diagrams>
<Paintings>
  <Text 40 20 11 #000000 0 "Non-Ideal Transformer Circuit\n--------------------------------------------------------\nA linear two-winding transformer can be represented by two mutual inductance coupled \ncircuit loops linking the transformer's five impedance constants.\n\nM is mutual inductance\nL1 & L1 are primary and secondary winding self-inductances\nR1 & R2 are primary and secondary winding resistances\n\nConstants M, L1, L2, R1 & R2 are measurable at the transformer's terminals\n\nThe coupling coefficient must be specified, but is k can be calculated from:\n\n            k = M / sqrt{L_P * L_S}, with 0 < k < 1\n\nTypical values for k are in the region 0.99 -- 0.98. Here you specify the parameters Lm, the\nprimary winding inductance and V_ratio, the output voltage ratio. The secondary winding \nratio is determined from the equation:\n\n            L2 = L1 * V_ratio">
</Paintings>
