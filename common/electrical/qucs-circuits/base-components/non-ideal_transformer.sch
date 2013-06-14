<Qucs Schematic 0.0.16>
<Properties>
  <View=0,0,974,800,1,0,0>
  <Grid=10,10,1>
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
  <.ID -120 -66 NTR "1=R1_val=1=Primary Winding Resistiance=Resistance" "1=R2_val=1=Secondary Winding Resistance=Resistance" "1=L1_val=0.1=Primary Leakage Inductance=Inductance" "1=L2_val=0.1=Secondary Leakage Inductance=Inductance" "1=Lm_val=0.01=Magnetising Inductance=Inductance" "1=Rcl_val=0.01=Core Loss Resistance=Resistance" "1=V_ratio=1=Primary to Secondary Voltage Ratio=Unitless">
  <Line -120 -70 80 0 #000080 2 1>
  <Line -120 -150 0 80 #000080 2 1>
  <Line -120 -150 80 0 #000080 2 1>
  <Line -40 -150 0 80 #000080 2 1>
  <Line -40 -140 10 0 #000080 2 1>
  <.PortSym -30 -140 2 180>
  <Line -40 -80 10 0 #000080 2 1>
  <.PortSym -30 -80 4 180>
  <.PortSym -130 -140 1 0>
  <Line -130 -140 10 0 #000080 2 1>
  <Line -130 -80 10 0 #000080 2 1>
  <.PortSym -130 -80 3 0>
  <Ellipse -60 -140 10 10 #000000 0 1 #000000 1 1>
  <Ellipse -110 -140 10 10 #000000 0 1 #000000 1 1>
</Symbol>
<Components>
  <Port P1 1 70 260 -23 -56 1 0 "1" 1 "analog" 0>
  <Port S2 1 880 420 4 12 1 2 "4" 1 "analog" 0>
  <Port S1 1 880 260 4 -56 0 2 "2" 1 "analog" 0>
  <Port P2 1 70 420 -23 12 0 0 "3" 1 "analog" 0>
  <R R1 1 160 260 -26 15 0 0 "R1_val" 1 "26.85" 0 "0.0" 0 "0.0" 0 "26.85" 0 "european" 0>
  <L L1 1 250 260 -26 10 0 0 "L1_val" 1 "" 0>
  <R Rcl 1 310 330 15 -26 0 1 "Rcl_val" 1 "26.85" 0 "0.0" 0 "0.0" 0 "26.85" 0 "european" 0>
  <L Lm 1 430 330 10 -26 0 1 "Lm_val" 1 "" 0>
  <R R2 1 670 260 -26 15 0 0 "R2_val" 1 "26.85" 0 "0.0" 0 "0.0" 0 "26.85" 0 "european" 0>
  <L L2 1 780 260 -26 10 0 0 "L2_val" 1 "" 0>
  <Tr Tr1 1 560 330 -29 38 0 0 "V_ratio" 1>
</Components>
<Wires>
  <700 260 750 260 "" 0 0 0 "">
  <810 260 880 260 "" 0 0 0 "">
  <190 260 220 260 "" 0 0 0 "">
  <70 260 130 260 "" 0 0 0 "">
  <280 260 310 260 "" 0 0 0 "">
  <310 260 430 260 "" 0 0 0 "">
  <310 260 310 300 "" 0 0 0 "">
  <70 420 310 420 "" 0 0 0 "">
  <310 420 430 420 "" 0 0 0 "">
  <310 360 310 420 "" 0 0 0 "">
  <430 260 430 300 "" 0 0 0 "">
  <430 360 430 420 "" 0 0 0 "">
  <430 420 530 420 "" 0 0 0 "">
  <530 360 530 420 "" 0 0 0 "">
  <590 420 880 420 "" 0 0 0 "">
  <590 360 590 420 "" 0 0 0 "">
  <590 260 640 260 "" 0 0 0 "">
  <590 260 590 300 "" 0 0 0 "">
  <430 260 530 260 "" 0 0 0 "">
  <530 260 530 300 "" 0 0 0 "">
</Wires>
<Diagrams>
</Diagrams>
<Paintings>
  <Text 40 30 11 #000000 0 "Non-Ideal Transformer SubCircuit\n--------------------------------------------------------\nThis circuit is a parameterised non-ideal transformer equivalent circuit. R1 and R2 are the main\nwinding resistances for the primary and secondary windings respectively. L1 and L2 are the  \nleakage inductances for the windings. Rcl is used to represent the losses in the magnetic core \nand Lm the magnetising reactance. \n\nSee Electric Machinery and Transformers by Guru and Hizeroglu or similar textbooks.">
</Paintings>
