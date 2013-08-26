<Qucs Schematic 0.0.17>
<Properties>
  <View=0,0,954,800,1,0,0>
  <Grid=5,5,1>
  <DataSet=odd_harmonic_voltage_source.dat>
  <DataDisplay=odd_harmonic_voltage_source.dpl>
  <OpenDisplay=1>
  <Script=odd_harmonic_voltage_source.m>
  <RunScript=0>
  <showFrame=0>
  <FrameText0=Title>
  <FrameText1=Drawn By:>
  <FrameText2=Date:>
  <FrameText3=Revision:>
</Properties>
<Symbol>
  <Ellipse -30 -150 40 40 #000000 0 1 #c0c0c0 1 0>
  <.PortSym -40 -130 1 0>
  <Line -40 -130 10 0 #000080 2 1>
  <.ID -30 -106 HV "1=f=50=Base Frequency (primary harmonic)=Frequency" "1=V0=0=DC voltage component amplitude=Voltage" "1=V1=1=1st Harmonic Amplitude=Voltage" "1=V3=0.5=3rd Harmonic Amplitude=Voltage" "1=V5=0=5th Harmonic Amplitude=Voltage" "0=V7=0=7th Harmonic Amplitude=Voltage" "0=V9=0=9th Harmonic Amplitude=Voltage" "0=V11=0=11th Harmonic Amplitude=Voltage" "0=V13=0=13th Harmonic Amplitude=Voltage">
  <Ellipse -5 -145 30 30 #000000 0 1 #c0c0c0 1 0>
  <.PortSym 50 -130 2 180>
  <Line 40 -130 10 0 #000080 2 1>
  <Ellipse 20 -140 20 20 #000000 0 1 #c0c0c0 1 0>
  <Line -35 -145 0 10 #aa0000 2 1>
  <Line -40 -140 10 0 #aa0000 2 1>
  <Line 40 -140 10 0 #aa0000 2 1>
</Symbol>
<Components>
  <Port P1 1 90 280 -54 4 1 1 "1" 1 "analog" 0>
  <Port P2 1 840 280 12 4 0 1 "2" 1 "analog" 0>
  <Vdc VH0 1 160 260 -26 -50 0 2 "V0" 1>
  <Eqn FrequenciesEqn 1 110 355 -24 15 0 0 "f3=f * 3" 1 "f5=f * 5" 1 "f7=f * 7" 1 "f9=f * 9" 1 "f11=f * 11" 1 "f13=f * 13" 1 "no" 0>
  <Vac VH9 1 560 260 -26 -66 0 2 "V9" 1 "f9" 1 "0" 0 "0" 0>
  <Vac VH5 1 460 260 -26 -66 0 2 "V5" 1 "f5" 1 "0" 0 "0" 0>
  <Vac VH3 1 360 260 -26 -66 0 2 "V3" 1 "f3" 1 "0" 0 "0" 0>
  <Vac VH1 1 260 260 -26 -66 0 2 "V1" 1 "f" 1 "0" 0 "0" 0>
  <Vac VH11 1 660 260 -26 -66 0 2 "V11" 1 "f11" 1 "0" 0 "0" 0>
  <Vac VH13 1 760 260 -26 -66 0 2 "V13" 1 "f13" 1 "0" 0 "0" 0>
</Components>
<Wires>
  <290 260 330 260 "" 0 0 0 "">
  <390 260 430 260 "" 0 0 0 "">
  <490 260 530 260 "" 0 0 0 "">
  <590 260 630 260 "" 0 0 0 "">
  <690 260 730 260 "" 0 0 0 "">
  <190 260 230 260 "" 0 0 0 "">
  <90 260 130 260 "" 0 0 0 "">
  <90 260 90 280 "" 0 0 0 "">
  <790 260 840 260 "" 0 0 0 "">
  <840 260 840 280 "" 0 0 0 "">
</Wires>
<Diagrams>
</Diagrams>
<Paintings>
  <Text 60 40 12 #000000 0 "Odd Harmonic Voltage Source\n-------------------------------------------\n\nVoltage source subcircuit which produces the first 6 odd harmonics as might\nbe obtained from a fourier analysis.">
</Paintings>
