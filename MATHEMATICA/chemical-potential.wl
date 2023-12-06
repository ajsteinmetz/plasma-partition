(* ::Package:: *)

(* ::Section:: *)
(*Chemical Potential Setup*)


(* ::Text:: *)
(*Here we define the physical and cosmological constants requires to evaluate the early universe electron-positron gas.*)


g = 2; (*Gyro-magnetic factor (g-factor) of electrons*)
me = 511; (*Mass of electrons in [keV]*)
loT = 10; hiT = 300; (*Cosmic temperature range in [keV]*)
Xp = 0.878; (*Baryon fraction within protons (excluding neutrons)*)
nBs = 0.865 10^(-10); (*Comoving baryon density to entropy density ratio*)
gs = 3.91; (*Relativistic fermion and boson degrees of freedom*)


(* ::Text:: *)
(*We define the functional dependencies which are fed into the chemical potential.*)


(* ::Item:: *)
(*The x function is the dimensionless Boltzmann factor "E/T". *)


(* ::Item:: *)
(*The \[Omega] function is the Bessel function which defines the dominant contribution to the chemical potential for the electron-positron gas.*)


(* ::Item:: *)
(*The "LeftHandSide" defines the trig. dependance on chemical potential given  by the charge neutrality of the universe.*)


x[T_, b_, s_, g_] := Sqrt[m^(2)/T^(2) + b(1 - s g/2)];
\[Omega][T_, b_, s_, g_]:=x[T, b, s, g]^(2)BesselK[2,x[T, b, s, g]]+\
(b/2)x[T, b, s, g]BesselK[1,x[T, b, s, g]]+\
(b^(2)/12)BesselK[0,x[T, b, s, g]];
LeftHandSide[T_,s_,M_,\[Eta]_]:=Sinh[M-s*\[Eta]/T];


(* ::Section:: *)
(*Chemical Potential Plotting*)


(* ::Text:: *)
(*Define some housekeeping functions for the axis.*)


reversal[T_]:=-T+ loT + hiT;
Tlist[n_]:=Table[x,{x,n,9n,n}];


(* ::Text:: *)
(*Make a table of contours to fit. E.g. Find where LeftHandSide equals RightHandSide.*)


(* ::Item:: *)
(*Table1: The chemical potential function for various cosmic field strengths*)


(* ::Subitem:: *)
(*Style1: Blue dotted lines*)


(* ::Item:: *)
(*Table2: The chemical potential function for various spin polarization values*)


(* ::Subitem:: *)
(*Style2: Black dashed lines*)


(* ::Item:: *)
(*Table3: The free particle chemical potential function with only thermal dependency*)


(* ::Subitem:: *)
(*Style3: Black solid line*)


(*Define the curves to be printed*)
table1 = Table[LeftHandSide[T,1,M,0]\[Omega][T,b,1,g] + LeftHandSide[T,-1,M,0]\[Omega][T,b,-1,g],{b,{25,50,100,300}}]/.m->me;
table2 = Table[LeftHandSide[T,1,M,\[Eta]]\[Omega][T,0,1,g] + LeftHandSide[T,-1,M,\[Eta]]\[Omega][T,0,-1,g],{\[Eta],{100,200,300,400,500}}]/.m->me;
table3 = Table[LeftHandSide[T,1,M,0]\[Omega][T,b,1,g] + LeftHandSide[T,-1,M,0]\[Omega][T,b,-1,g],{b,{0}}]/.m->me;
(*Define line styles for each curve*)
style1 = ContourStyle -> {Directive[Blue,Dotted]};
style2 = ContourStyle -> {Directive[Black,Dashed]};
style3 = ContourStyle -> {Directive[Black]};


(* ::Text:: *)
(*This is the shared plotting code used by all contours.*)


potentialplot[table_,style_] := ContourPlot[
 (Pi^(2) Xp nBs (2 Pi^(2) gs)/45) == table,
 {T, loT, hiT},{M,10^(-12),10^(2)},
 PlotRange -> {{10,300},{10^(3),10^(-12)}},Frame -> True,
 AspectRatio -> 3/2,
 FrameLabel -> {"T [keV]", "\[Mu]/T"},
 LabelStyle -> Directive[Black,14,FontFamily -> "Times"],
 style,
 PlotPoints->30,
 FrameStyle -> Directive[Black,20],
 Background -> White,
 ScalingFunctions -> {
       {-Log[#]&,InverseFunction[-Log[#]&]},
       {Log,InverseFunction[Log]}},
 FrameTicks -> {{{#, Superscript[10, Log10@#]} & /@ ({10^2, 10^-2, 10^-6, 10^-10}), None},
       {{#, Superscript[10, Log10@#]} & /@ ({10^3, 10^2, 10^1}), None}},
 GridLines -> {Drop[Flatten[Table[Tlist[n],{n,{10,100,1000}}]],-8],Table[10^(n),{n,1,-33,-1}]},
 GridLinesStyle -> Directive[Line, Lighter[Gray,.65]]];


(* ::Text:: *)
(*Print the final plot.*)


Show[potentialplot[table1,style1],potentialplot[table2,style2],potentialplot[table3,style3],PlotRange->All]
