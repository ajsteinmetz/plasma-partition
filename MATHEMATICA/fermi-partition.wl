(* ::Package:: *)

(* ::Chapter:: *)
(*Fermi Partition Function*)


(* ::Text:: *)
(*This is a notebook which defines the partition function for the magnetized Fermi gas in the Boltzmann approximation and calculates various thermodynamic quantities of that gas. The planned application for this is to describe the electron-positron plasma of the early universe between the temperatures of 20 keV and 2,000 keV.*)


(* ::Section:: *)
(*Initial Setup*)


(* ::Subsection:: *)
(*Definitions*)


(* ::Text:: *)
(*Define the Bessel parameter X.*)


x[T_, b_, s_, g_] := Sqrt[m^(2)/T^(2) + b (1 - s g/2)];


x[T,b,-1,2]//TraditionalForm;
x[T,b,1,2]//TraditionalForm;


(* ::Text:: *)
(*Define the magnetic field strength*)


B[T_, b_] := (T/511)^2 b


(* ::Text:: *)
(*Define the spin potential and spin fugacity.*)


\[Xi][T_,s_,\[Eta]_]:=Exp[s \[Eta]/T];


(* ::Text:: *)
(*Define the partition function.*)


lnZ[V_, u_, T_, b_, s_, g_,\[Eta]_] := (T^(3) V/(2 Pi)^(2)) (2 Cosh[u/T])\[Xi][T,s,\[Eta]](x[T, b, s, g]^(2) BesselK[2, x[T, b, s, g]] + \
              b x[T, b, s, g] BesselK[1, x[T, b, s, g]]/2 + b^(2) BesselK[0, x[T, b, s, g]]/12);


lnZ[V,u,T,b,-1,2,0]+lnZ[V,u,T,b,1,2,0]//TraditionalForm;
lnZ[V,u,T,b,-1,2,0]+lnZ[V,u,T,b,1,2,0]//FullSimplify//TraditionalForm;


(* ::Subsection::Closed:: *)
(*Zero B-field limit test*)


(T/V)(q/T^(2))D[lnZ[V,u,T,b,1,2,0],b]/.b-> 0/.Sqrt[m^2/T^2]->m/T//FullSimplify//TraditionalForm;
(T/V)(q/T^(2))D[lnZ[V,u,T,b,-1,2,0],b]/.b-> 0/.Sqrt[m^2/T^2]->m/T//FullSimplify//TraditionalForm;


(* ::Text:: *)
(*Check that b=0 limit magnetization is zero.*)


Clear[\[Eta]];
(q^(2)T^(2)/(2 Pi^(2)m^(2)))^(-1)(q/m^(2))(T/V)(q/T^(2))(D[lnZ[V,u,T,b,1,2,\[Eta]],b])/.Sqrt[m^2/T^2]->m/T/.b->0//FullSimplify//TraditionalForm
(q^(2)T^(2)/(2 Pi^(2)m^(2)))^(-1)(q/m^(2))(T/V)(q/T^(2))(D[lnZ[V,u,T,b,-1,2,\[Eta]],b])/.Sqrt[m^2/T^2]->m/T/.b->0//FullSimplify//TraditionalForm
(q^(2)T^(2)/(2 Pi^(2)m^(2)))^(-1)(q/m^(2))(T/V)(q/T^(2))(D[lnZ[V,u,T,b,1,2,\[Eta]],b]+D[lnZ[V,u,T,b,-1,2,\[Eta]],b])/.b->0/.Sqrt[m^2/T^2]->m/T//FullSimplify//TraditionalForm


(* ::Text:: *)
(*Introduce spin potential as an exponential weight. Doing this, I expect that the spin potential should manifest as a Sinh[us/T] term in the partition function. Here is my thinking: M[=]B/Subscript[\[Mu], 0] and we have b=eB/T^2. So we then need for mean field self magnetization M->Subscript[M, m]=Subscript[B, m]/Subscript[\[Mu], 0]=bT^2/Subscript[e\[Mu], 0].*)


(m q T BesselK[1, m/T] Cosh[u/T])/(4 \[Pi]^2) Exp[us/T] - (m q T BesselK[1, m/T] Cosh[u/T])/(4 \[Pi]^2) Exp[-us/T] // FullSimplify


(* ::Section:: *)
(*Magnetization*)


(* ::Subsection:: *)
(*Evaluating magnetization*)


(* ::Text:: *)
(*Perform the magnetization calculation and sum the two polarization*)


(T/V)(q/T^(2))(D[lnZ[V,u,T,b,1,2,0],b]+\
D[lnZ[V,u,T,b,-1,2,0],b])/.Sqrt[m^2/T^2]->m/T//FullSimplify//TraditionalForm


(* ::Text:: *)
(*Let' s determine the positive polarization magnetization (with g = 2). We've divided out the coefficient (q^2/2 \[Pi]^2)(T^2/m^2) shared between the two expressions. We've also introduced the spin fugacity \[Xi].*)


(q^(2)T^(2)/(2 Pi^(2)m^(2)))^(-1)(q/m^(2))(T/V)(q/T^(2))(D[lnZ[V,u,T,b,1,2,0],b])/.Sqrt[m^2/T^2]->m/T//FullSimplify//TraditionalForm
(q^(2)T^(2)/(2 Pi^(2)m^(2)))^(-1)(q/m^(2))(T/V)(q/T^(2))(D[lnZ[V,u,T,b,-1,2,0],b])/.Sqrt[m^2/T^2]->m/T//FullSimplify//TraditionalForm
(q^(2)T^(2)/(2 Pi^(2)m^(2)))^(-1)(q/m^(2))(T/V)(q/T^(2))(D[lnZ[V,u,T,b,1,2,\[Eta]],b]+D[lnZ[V,u,T,b,-1,2,\[Eta]],b])/.Sqrt[m^2/T^2]->m/T//FullSimplify//TraditionalForm


(* ::Subsection:: *)
(*Plotting the Magnetization*)


(* ::Text:: *)
(*Take the output magnetization, divide by Cosh[u/T], and plot over temperature range. Replace q->4\[Pi] \[Alpha]/(Subscript[m, e])^2 to obtain unit-less magnetization using M/Subscript[H, C]. Sub in electron mass in [keV] and cosmic magnetic scale. Critical magnetization defined as Subscript[H, C]=Subscript[B, C]/Subscript[\[Mu], 0] where we use the critical field limit Subscript[B, C]=(Subscript[m, e] c^2)^2/(2\[HBar] c^2). Technically, we're not plotting magnetization proper because we've divided out the chemical potential portion Cosh[\[Mu]/T]. Full plot would need to put the chemical potential part back in.*)


g = {1.18, 2, 5};
g = 2;
\[Eta] = 0;
f[T_, b_] = (T/V) (q/T^(2)) (D[lnZ[V, u, T, b, 1, g, \[Eta]], b] + D[lnZ[V, u, T, b, -1, g, \[Eta]], b])/Cosh[u/T] /. q -> (4 Pi/137)/m^(2) /. m -> 511;
bValues = {10^(-3), 10^(-11)};
LogLogPlot[Evaluate[f[T, #] & /@ bValues], {T, 20, 2000},
 	PlotRange -> All,
 	AxesLabel -> {"T [keV]", "f[T]/\!\(\*SubscriptBox[\(H\), \(c\)]\)"},
 	PlotStyle -> {Blue, Red},
 	PlotLegends -> Placed[LineLegend[{"b = \!\(\*SuperscriptBox[\(10\), \(-3\)]\)", "b = \!\(\*SuperscriptBox[\(10\), \(-11\)]\)"}], {0.8, 0.4}],
 	PlotLabel -> "Plot of f[T] from T = 20 to 2,000 keV",
 	GridLines -> Automatic, GridLinesStyle -> LightGray]


(* ::Text:: *)
(*Let' s also plot the external fields over the same temperature range. The expressions for the external fields in dimensional units are as follows:*)


LogLogPlot[{B[T, 10^(-3)], B[T, 10^(-11)]}, {T, 20, 2000},
  PlotStyle -> {Directive[Dashed, Blue], Directive[Dashed, Red]},
  PlotLegends -> Placed[LineLegend[{"b = \!\(\*SuperscriptBox[\(10\), \(-3\)]\)", "b = \!\(\*SuperscriptBox[\(10\), \(-11\)]\)"}], {0.8, 0.6}],
  AxesLabel -> {"T [keV]", "M[T]/\!\(\*SubscriptBox[\(H\), \(c\)]\)"},
  GridLines -> Automatic, GridLinesStyle -> LightGray,
  PlotRange -> All]


g = 2;
Bc = 4.41 10^(13);
\[Eta] = 0;
f[T_, b_] = (T/V) (q/T^(2)) (D[lnZ[V, u, T, b, 1, g, \[Eta]], b] + D[lnZ[V, u, T, b, -1, g, \[Eta]], b])/Cosh[u/T] /. q -> (4 Pi/137)/m^(2) /. m -> 511;
bValues = {10^(-3), 10^(-11)};
Plot[Evaluate[{Bc f[T, #],Bc B[T,#]} & /@ bValues], {T, 10, 2000},
 	PlotRange -> {{10,2000},{10^(-29),10^(14)}},
 	Frame -> True,
 	FrameLabel -> {"T [keV]", "B [G]"},
 	PlotStyle -> {Blue,Directive[Dashed,Blue],Red,Directive[Dashed, Red]},
 	PlotLegends -> Placed[LineLegend[{"M(\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-3\)]\))","H(\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-3\)]\))","M(\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-11\)]\))","H(\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-11\)]\))"}], {0.8, 0.3}],
 	Background -> White,
 	ScalingFunctions -> {"Log", "Log"},
 	GridLines -> Automatic, GridLinesStyle -> LightGray]


(* ::Subsection:: *)
(*Wrangling the behavior of temperature*)


(* ::Text:: *)
(*Here we plot the magnetization versus temperature for a fixed magnetic field, rather than fixed scale. Essentially, we're allowing the temperature to vary independently while the magnetic field strength remains constant. Below we see the magnetization at high temperature is suppressed compared to the comoving conserved magnetic flux case which has the magnetic field grow with temperature accordingly.*)


g = 2;
Bc = 4.41 10^(13);
Bc = 1;
\[Eta] = 0;
loT = 10; hiT = 2000;
Tlist[n_]:=Table[x,{x,n,9n,n}];
reversal[T_]:=-T+ loT + hiT;
f[T_, b_] = (T/V) (q/T^(2)) (D[lnZ[V, u, T, b, 1, g, \[Eta]], b] + D[lnZ[V, u, T, b, -1, g, \[Eta]], b])/Cosh[u/T] /. q -> (4 Pi/137)/m^(2) /. m -> 511;
bValues = {10^(-3), 10^(-3)*(1000)^(2)/T^(2)};
Plot[Evaluate[{Bc f[T, #],Bc B[T,#]} & /@ bValues], {T, loT, hiT},
 	PlotRange -> {{10,1000},{10^(-31),10^(1)}},
 	Frame -> True,
 	FrameLabel -> {"T [keV]", "\!\(\*OverscriptBox[\(\[ScriptCapitalM]\), \(_\)]\)"},
 	PlotStyle -> {Blue,Directive[Dashed,Blue],Red,Directive[Dashed, Red]},
 	PlotLegends -> Placed[LineLegend[{"\[ScriptCapitalM](\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-3\)]\))","\[ScriptCapitalB](\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-3\)]\))","\[ScriptCapitalM](\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-3\)]\)(\!\(\*SubscriptBox[\(T\), \(0\)]\)/T\!\(\*SuperscriptBox[\()\), \(2\)]\))","\[ScriptCapitalB](\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-3\)]\)(\!\(\*SubscriptBox[\(T\), \(0\)]\)/T\!\(\*SuperscriptBox[\()\), \(2\)]\))"}],{0.25, 0.25}],
 	LabelStyle -> Directive[Black,14,FontFamily -> "Times"],
 	FrameStyle -> Directive[Black,20],
 	Background -> White,
 	ScalingFunctions -> {
       {-Log[#]&,InverseFunction[-Log[#]&]},
       {Log,InverseFunction[Log]}},
 	GridLines -> {Drop[Flatten[Table[Tlist[n],{n,{10,100,1000}}]],-8],Table[10^(n),{n,1,-33,-1}]},
 	GridLinesStyle -> Directive[Line, Lighter[Gray,.9]]]


g = {2};
\[Eta] = 0;
bValues = {10^(-3), 10^(-3)*(20)^(2)/T^(2)};
f[T_, b_] = (T/V) (q/T^(2)) (D[lnZ[V, u, T, b, 1, g, \[Eta]], b] + D[lnZ[V, u, T, b, -1, g, \[Eta]], b])/Cosh[u/T];
LogLogPlot[Evaluate[f[T, #]/. q -> (4 Pi/137)/m^(2) /. m -> 511 & /@ bValues], {T,20,2000},
 	PlotRange -> All,
 	AxesLabel -> {"T [keV]", "f[T]/\!\(\*SubscriptBox[\(H\), \(c\)]\)"},
 	PlotStyle -> {Blue, Red},
 	PlotLegends -> Placed[LineLegend[{"b = \!\(\*SuperscriptBox[\(10\), \(-3\)]\)", "b = \!\(\*SuperscriptBox[\(10\), \(-3\)]\)(\!\(\*SubscriptBox[\(T\), \(0\)]\)/T\!\(\*SuperscriptBox[\()\), \(2\)]\)"}], {0.8, 0.4}],
 	PlotLabel -> "Plot of f[T] from T = 20 to 2,000 keV",
 	GridLines -> Automatic, GridLinesStyle -> LightGray]


(* ::Text:: *)
(*However, we still have this strange behavior that magnetization increases with temperature, which seems backwards. The magnetic field is no long being "boosted" by the FLRW scale factor a(t) as we travel to higher temperatures. The magnetization should in principle drop as the temperature increases.*)


(* ::Text:: *)
(*In the next input function we determine the magnetization as a function of magnetic field and temperature as independent quantities. The variable "b" is no longer cosmic scale, but magnetic field strength.*)


g = 2;
\[Eta] = 0;
f[T_, b_] = (T/V)(q/T^(2))(D[lnZ[V, u, T, b, 1, g, \[Eta]],b]+D[lnZ[V,u,T,b,-1, g, \[Eta]],b])/Cosh[u/T] /. b->b/T^(2)//FullSimplify//TraditionalForm


(* ::Text:: *)
(*And using this result, we plot the magnetization as a function of magnetic field at several fixed temperatures. Again note that the variable "b" is no longer cosmic scale, but magnetic field strength.*)


g = 2;
\[Eta] = 0;
f[T_, b_] = (T/V)(q/T^(2))(D[lnZ[V, u, T, b, 1, g, \[Eta]],b]+D[lnZ[V,u,T,b,-1, g, \[Eta]],b])/Cosh[u/T] /. q -> (4 Pi/137)/m^(2)/.b->b/T^(2)/.m->511;
TValues = {20, 50,100,200};
LogPlot[Evaluate[f[#, b] & /@ TValues], {b,0,10^(-3)},
 	PlotRange -> All,
 	AxesLabel -> {"eB", "f[eB]/\!\(\*SubscriptBox[\(H\), \(c\)]\)"},
 	PlotStyle -> {Blue,Red,Black,Orange},
 	PlotLegends -> Placed[LineLegend[{"T = 20 [keV]", "T = 50 [keV]","T = 100 [keV]","T = 200 [keV]"}], {0.8, 0.4}],
 	PlotLabel -> "Plot of f[eB] from eB = 0 to \!\(\*SuperscriptBox[\(10\), \(-3\)]\)",
 	GridLines -> Automatic, GridLinesStyle -> LightGray]


(* ::Text:: *)
(*The thing that is unclear to me is why increased temperatures lead to more magnetization rather than less. Is the partition function sick in some manner? Something feels wrong with this behavior. The Bessel K[x] functions mostly decay to zero for large values of x. Therefore, having x=m/T be the controlling variable leads to situations where for low temperature, the Bessel functions are killed, but for high temperature, they become large. To be fair, in the Boltzmann approximation, you do not expect the low temperature behavior to be correct as the Boltzmann distribution is very unlike the Fermi-Dirac distribution. I think the issue must be that the Boltzmann k=1 approximation must be only valid within a particular domain.*)


(* ::Section:: *)
(*Chemical Potential*)


(* ::Subsection:: *)
(*Definitions*)


(* ::Text:: *)
(*Define the chemical potential.*)


Xp = 0.878;
nBs = 0.865 10^(-10);
gs = 3.91;
entropy[T_]:=T^(3)(2 Pi^(2) gs)/45;
np[T_] := Xp nBs entropy[T];
chempot[T_,b_, g_,\[Eta]_] := ArcSinh[Pi^(2) Xp nBs (2 Pi^(2) gs)/45 (
			  \[Xi][T,1,\[Eta]](x[T, b, 1, g]^(2) BesselK[2, x[T, b, 1, g]] + b x[T, b, 1, g] BesselK[1, x[T, b, 1, g]]/2 + b^(2) BesselK[0, x[T, b, 1, g]]/12) + \
              \[Xi][T,-1,\[Eta]](x[T, b, -1, g]^(2) BesselK[2, x[T, b, -1, g]] + b x[T, b, -1, g] BesselK[1, x[T, b, -1, g]]/2 + b^(2) BesselK[0, x[T, b, -1, g]]/12))^(-1)];
chempot[T,0,2,10]


g = 2;
Bc = 4.41 10^(13);
Bc = 1;
\[Eta] = {0,100,200,300,400,500};
me = 511;
b0 = 10^(-3);
bValues = {25,50,100,300};
loT = 10; hiT = 1000;
Tlist[n_]:=Table[x,{x,n,9n,n}];
reversal[T_]:=-T+ loT + hiT;
bzero = {0};
etaValues = {0};
Plot[{Evaluate[chempot[T, #,g,\[Eta]]/.m->me & /@ bzero],Evaluate[chempot[T, #,g,etaValues]/.m->me & /@ bValues]}, {T, loT, hiT},
 	PlotRange -> {{10,200},{10^(2),10^(-12)}},
 	Frame -> True,
 	AspectRatio -> 3/2,
 	FrameLabel -> {"T [keV]", "\[Mu]/T"},
 	PlotStyle -> {Black,Directive[Dashed,Black],Directive[Dashed,Black],Directive[Dashed,Black],Directive[Dashed,Black],Directive[Dashed,Black],Directive[Dotted,Blue],Directive[Dotted,Blue],Directive[Dotted,Blue],Directive[Dotted,Blue],Directive[Dotted,Blue]},
 	LabelStyle -> Directive[Black,14,FontFamily -> "Times"],
 	FrameStyle -> Directive[Black,20],
 	Background -> White,
 	ScalingFunctions -> {
       {-Log[#]&,InverseFunction[-Log[#]&]},
       {Log,InverseFunction[Log]}},
     FrameTicks -> {{{#, Superscript[10, Log10@#]} & /@ ({10^2, 10^-2, 10^-6, 10^-10}), None},
       {{#, Superscript[10, Log10@#]} & /@ ({10^3, 10^2, 10^1}), None}},
 	GridLines -> {Drop[Flatten[Table[Tlist[n],{n,{10,100,1000}}]],-8],Table[10^(n),{n,1,-33,-1}]},
 	GridLinesStyle -> Directive[Line, Lighter[Gray,.8]]]


g = 2;
Bc = 4.41 10^(13);
Bc = 1;
\[Eta] = {0};
me = 511;
b0 = 10^(-3);
loT = 10; hiT = 1000;
Tlist[n_]:=Table[x,{x,n,9n,n}];
reversal[T_]:=-T+ loT + hiT;
bValues = {10,20,50,100};
Plot[Evaluate[chempot[T, #,g,\[Eta]]/.m->me & /@ bValues], {T, loT, hiT},
 	PlotRange -> {{10,200},{10^(2),10^(-12)}},
 	Frame -> True,
 	FrameLabel -> {"T [keV]", "\[Mu]/T"},
 	PlotStyle -> {Directive[Dashed,Blue],Directive[Dashed,Blue],Directive[Dashed,Blue],Directive[Dashed,Blue],Directive[Dashed,Blue]},
 	LabelStyle -> Directive[Black,14,FontFamily -> "Times"],
 	FrameStyle -> Directive[Black,20],
 	Background -> White,
 	ScalingFunctions -> {
       {-Log[#]&,InverseFunction[-Log[#]&]},
       {Log,InverseFunction[Log]}},
     FrameTicks -> {{{#, Superscript[10, Log10@#]} & /@ ({10^2, 10^-2, 10^-6, 10^-10}), None},
       {{#, Superscript[10, Log10@#]} & /@ ({10^3, 10^2, 10^1}), None}},
 	GridLines -> {Drop[Flatten[Table[Tlist[n],{n,{10,100,1000}}]],-8],Table[10^(n),{n,1,-33,-1}]},
 	GridLinesStyle -> Directive[Line, Lighter[Gray,.8]]]


(* ::InheritFromParent:: *)
(**)


g = {2};
Bc = 4.41 10^(13);
Bc = 1;
\[Eta] = 0;
loT = 10; hiT = 2000;
Tlist[n_]:=Table[x,{x,n,9n,n}];
reversal[T_]:=-T+ loT + hiT;
f[T_, b_] = (T/V) (q/T^(2)) (D[lnZ[V, chempot[T,0,1,g]T, T, b, 1, g, \[Eta]], b] + D[lnZ[V, chempot[T,0,-1,g]T, T, b, -1, g, \[Eta]], b])/. q -> (4 Pi/137)/m^(2) /. m -> 511;
bValues = {10^(-3), 10^(-3)*(10)^(2)/T^(2)};
Plot[Evaluate[{Bc f[T, #],Bc B[T,#]} & /@ bValues], {T, loT, hiT},
 	PlotRange -> {{10,1000},{10^(-25),10^(1)}},
 	Frame -> True,
 	FrameLabel -> {"T [keV]", "\!\(\*OverscriptBox[\(\[ScriptCapitalM]\), \(_\)]\)"},
 	PlotStyle -> {Blue,Directive[Dashed,Blue],Red,Directive[Dashed, Red]},
 	PlotLegends -> Placed[LineLegend[{"\[ScriptCapitalM](\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-3\)]\))","\[ScriptCapitalB](\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-3\)]\))","\[ScriptCapitalM](\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-3\)]\)(\!\(\*SubscriptBox[\(T\), \(0\)]\)/T\!\(\*SuperscriptBox[\()\), \(2\)]\))","\[ScriptCapitalB](\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-3\)]\)(\!\(\*SubscriptBox[\(T\), \(0\)]\)/T\!\(\*SuperscriptBox[\()\), \(2\)]\))"}],{0.25, 0.25}],
 	LabelStyle -> Directive[Black,14,FontFamily -> "Times"],
 	FrameStyle -> Directive[Black,20],
 	Background -> White,
 	ScalingFunctions -> {
       {-Log[#]&,InverseFunction[-Log[#]&]},
       {Log,InverseFunction[Log]}},
     FrameTicks -> {{{#, Superscript[10, Log10@#]} & /@ ({10^0, 10^-5, 10^-10, 10^-15, 10^-20, 10^-25}), None},
       {{#, Superscript[10, Log10@#]} & /@ ({10^3, 10^2, 10^1}), None}},
 	GridLines -> {Drop[Flatten[Table[Tlist[n],{n,{10,100,1000}}]],-8],Table[10^(n),{n,1,-33,-1}]},
 	GridLinesStyle -> Directive[Line, Lighter[Gray,.8]]]


ScientificForm[{10.,100.,10.^(-2),1000000.}]
ScientificForm[{123450000.0,0.000012345,123.45}]


g = {2};
Bc = 4.41 10^(13);
Bc = 1;
\[Eta] = 0;
etaValues = {10^-9,10^-6,10^-3,10^-2,10^-1,10^0,10^1};
loT = 10; hiT = 2000;
Tlist[n_]:=Table[x,{x,n,9n,n}];
reversal[T_]:=-T+ loT + hiT;
ff[T_, b_,eta_] = (T/V) (q/T^(2)) (D[lnZ[V, chempot[T,0,g,\[Eta]]T, T, b, 1, g, eta], b] + D[lnZ[V, chempot[T,0,g,\[Eta]]T, T, b, -1, g, eta], b])/. q -> (4 Pi/137)/m^(2) /. m -> 511;
bValues = {10^(-3),10^(-11)};
Plot[Evaluate[{Bc ff[T, bValues, 0],Bc ff[T, 0, #] & /@ etaValues}], {T, loT, hiT},
 	PlotRange -> {{10,1000},{10^(-30),10^(0)}},
 	Frame -> True,
 	AspectRatio -> 3/2,
 	FrameLabel -> {"T [keV]", "M"},
 	PlotStyle -> {Blue,Red,Directive[Dashed,Black],Directive[Dashed,Black],Directive[Dashed,Black],Directive[Dashed,Black],Directive[Dashed,Black],Directive[Dashed,Black],Directive[Dashed,Black]},
 	LabelStyle -> Directive[Black,18,FontFamily -> "Times"],
 	FrameStyle -> Directive[Black,24],
 	Background -> White,
 	ScalingFunctions -> {
       {-Log[#]&,InverseFunction[-Log[#]&]},
       {Log,InverseFunction[Log]}},
     FrameTicks -> {{{#, Superscript[10, Log10@#]} & /@ ({10^0, 10^-5, 10^-10, 10^-15, 10^-20, 10^-25,10^-30}), None},
       {{#, Superscript[10, Log10@#]} & /@ ({10^3, 10^2, 10^1}), None}},
 	GridLines -> {Drop[Flatten[Table[Tlist[n],{n,{10,100,1000}}]],-8],Table[10^(n),{n,1,-33,-1}]},
 	GridLinesStyle -> Directive[Line, Lighter[Gray,.8]]]


g = {2};
Bc = 4.41 10^(13);
Bc = 1;
\[Eta] = 0;
bValues = {10^(-3)};
Tlist[n_]:=Table[x,{x,n,9n,n}];
ff[T_, b_,eta_,gee_] = (T/V) (q/T^(2)) (D[lnZ[V, chempot[T,0,gee,\[Eta]]T, T, b, 1, gee, eta], b] + D[lnZ[V, chempot[T,0,gee,\[Eta]]T, T, b, -1, gee, eta], b])/. q -> (4 Pi/137)/m^(2) /. m -> 511;
Plot[{Bc ff[150, bValues, 0, gee],Bc ff[120, bValues, 0, gee],Bc ff[100, bValues, 0, gee],Bc ff[80, bValues, 0, gee]}, {gee, -4, 4},
 	PlotRange -> {Automatic,{-2 10^(-9),10^(-8)}},
 	Frame -> True,
 	AspectRatio -> 3/2,
 	FrameLabel -> {"g-factor", "\[GothicCapitalM]"},
 	PlotStyle -> {Black,Red,Blue,Green},
 	LabelStyle -> Directive[Black,24,FontFamily -> "Times"],
 	FrameStyle -> Directive[Black,22],
 	Background -> White,
     GridLines->Automatic,
 	GridLinesStyle -> Directive[Line, Lighter[Gray,.8]]]


(* ::Section:: *)
(*Analysis of the pair effect*)


(* ::Subsection:: *)
(*Spin fugacity sensitivity*)


(* ::Text:: *)
(*If there is even a small asymmetry in the spins of the electron-positron gas, the gas is easily magnetized (and even self magnetizing) for small seed primordial fields. The gas with a small spin asymmetry and small PMF self magnetizes in a similar manner to the non-self magnetized spin symmetric but stronger PMF case. The two become indistinguishable as the temperature cools and the plasma vanishes. The results here, and in the below case, all presume that the spin potential is NOT parameterized to temperature. Whether this is a good assumption is not yet determined.*)


g = 2;
Bc = 4.41 10^(13);
f[T_, b_,eta_] = (T/V) (q/T^(2)) (D[lnZ[V, u, T, b, 1, g, eta], b] + D[lnZ[V, u, T, b, -1, g, eta], b])/Cosh[u/T] /. q -> (4 Pi/137)/m^(2) /. m -> 511;
bValues = {10^(-3), 10^(-11)};
Manipulate[Plot[Evaluate[{Bc f[T, #,eta],Bc B[T,#]} & /@ bValues], {T, 10, 2000},
 	PlotRange -> {{10,2000},{10^(-29),10^(14)}},
 	Frame -> True,
 	FrameLabel -> {"T [keV]", "B [G]"},
 	PlotStyle -> {Blue,Directive[Dashed,Blue],Red,Directive[Dashed, Red]},
 	PlotLegends -> Placed[LineLegend[{"M(\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-3\)]\))","H(\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-3\)]\))","M(\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-11\)]\))","H(\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-11\)]\))"}], {0.8, 0.3}],
 	Background -> White,
 	ScalingFunctions -> {"Log", "Log"},
 	GridLines -> Automatic, GridLinesStyle -> LightGray],{eta,0,10^(-3),10^(-4),Appearance->"Labeled"}]


(* ::Text:: *)
(*Conversely, the anti-parallel asymmetric case (with \[Eta]<0) presents very differently between the strong PMF case and the weak PMF case. In the strong case, there exists a critical temperature where, for a given spin potential, the magnetization flips from paramagnetic to diamagnetic. The switch is sudden and dramatic as the whole gas reorients itself to a whole new magnetization. In the weak PMF case, the case is either fully paramagnetic or fully diamagnetic with the spin potential controlling which phase the gas is in over the whole temperature range.*)


g = 2;
Bc = 4.41 10^(13);
f[T_, b_,eta_] = (T/V) (q/T^(2)) (D[lnZ[V, u, T, b, 1, g, eta], b] + D[lnZ[V, u, T, b, -1, g, eta], b])/Cosh[u/T] /. q -> (4 Pi/137)/m^(2) /. m -> 511;
bValues = {10^(-3), 10^(-11)};
Manipulate[Plot[Evaluate[{Bc Abs[f[T, #,eta]],Bc B[T,#]} & /@ bValues], {T, 10, 2000},
 	PlotRange -> {{10,2000},{10^(-29),10^(14)}},
 	Frame -> True,
 	FrameLabel -> {"T [keV]", "B [G]"},
 	PlotStyle -> {Blue,Directive[Dashed,Blue],Red,Directive[Dashed, Red]},
 	PlotLegends -> Placed[LineLegend[{"M(\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-3\)]\))","H(\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-3\)]\))","M(\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-11\)]\))","H(\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-11\)]\))"}], {0.8, 0.3}],
 	Background -> White,
 	ScalingFunctions -> {"Log", "Log"},
 	GridLines -> Automatic, GridLinesStyle -> LightGray],{eta,-10^(-3),10^(-3),10^(-4),Appearance->"Labeled"}]


(* ::Text:: *)
(*In the plot below we explore the possibility for self-magnetization with zero externally applied magnetic flux.*)


g = 2;
Bc = 4.41 10^(13);
Bc = 1;
\[Eta] = 0;
f[T_, b_,eta_] = (T/V) (q/T^(2)) (D[lnZ[V, chempot[T,0,g,\[Eta]]T, T, b, 1, g, eta], b] + D[lnZ[V, chempot[T,0,g,\[Eta]]T, T, b, -1, g, eta], b]) /. q -> (4 Pi/137)/m^(2) /. m -> 511;
bValues = {10^(-3), 0};
Manipulate[Plot[{Bc f[T, 10^(-3),0],Evaluate[{Bc f[T, #,eta],Bc B[T,#]} & /@ bValues]}, {T, 10, 2000},
 	PlotRange -> {{10,1000},{10^(-25),10^(1)}},
 	Frame -> True,
 	FrameLabel -> {"T [keV]", "\!\(\*OverscriptBox[\(\[ScriptCapitalM]\), \(_\)]\)"},
 	PlotStyle -> {Directive[Black,Dashed],Blue,Directive[Blue,Dashed],Red},
 	PlotLegends -> Placed[LineLegend[{"\[ScriptCapitalM](\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-3\)]\); \[Eta]=0)","\[ScriptCapitalM](\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-3\)]\); \[Eta]=\!\(\*SuperscriptBox[\(10\), \(-3\)]\) [keV])","\[ScriptCapitalB](\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-3\)]\))","\[ScriptCapitalM](\!\(\*SubscriptBox[\(b\), \(0\)]\)=0; \[Eta]=\!\(\*SuperscriptBox[\(10\), \(-3\)]\) [keV])"}],{0.25, 0.25}],
 	LabelStyle -> Directive[Black,14,FontFamily -> "Times"],
 	FrameStyle -> Directive[Black,20],
 	Background -> White,
 	ScalingFunctions -> {
       {-Log[#]&,InverseFunction[-Log[#]&]},
       {Log,InverseFunction[Log]}},
     FrameTicks -> {{{#, Superscript[10, Log10@#]} & /@ ({10^0, 10^-5, 10^-10, 10^-15, 10^-20, 10^-25}), None},
       {{#, Superscript[10, Log10@#]} & /@ ({10^3, 10^2, 10^1}), None}},
 	GridLines -> {Drop[Flatten[Table[Tlist[n],{n,{10,100,1000}}]],-8],Table[10^(n),{n,1,-33,-1}]},
 	GridLinesStyle -> Directive[Line, Lighter[Gray,.8]]],{eta,0,10^(-3),10^(-4),Appearance->"Labeled"}]


g = 2;
Bc = 1;
eta1 = {0,10^(-3)};
f[T_, b_,eta_] = (T/V) (q/T^(2)) (D[lnZ[V, chempot[T,0,1,g]T, T, b, 1, g, eta], b] + D[lnZ[V, chempot[T,0,-1,g]T, T, b, -1, g, eta], b]) /. q -> (4 Pi/137)/m^(2) /. m -> 511;
bValues = {10^(-3)};
Plot[Bc f[T, bValues,eta1], {T, 10, 10000},
 	PlotRange -> {{10,10000},{10^(-25),10^(1)}},
 	Frame -> True,
 	FrameLabel -> {"T [keV]", "\!\(\*OverscriptBox[\(\[ScriptCapitalM]\), \(_\)]\)"},
 	PlotStyle -> {Directive[Blue,Dashed],Directive[Black,Dashed],Blue,Red},
 	PlotLegends -> Placed[LineLegend[{"\[ScriptCapitalM](\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-3\)]\); \[Eta]=0)","\[ScriptCapitalM](\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-3\)]\); \[Eta]=\!\(\*SuperscriptBox[\(10\), \(-3\)]\) [keV])","\[ScriptCapitalB](\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-3\)]\))","\[ScriptCapitalM](\!\(\*SubscriptBox[\(b\), \(0\)]\)=0; \[Eta]=\!\(\*SuperscriptBox[\(10\), \(-3\)]\) [keV])"}],{0.25, 0.25}],
 	LabelStyle -> Directive[Black,14,FontFamily -> "Times"],
 	FrameStyle -> Directive[Black,20],
 	Background -> White,
 	ScalingFunctions -> {
       {-Log[#]&,InverseFunction[-Log[#]&]},
       {Log,InverseFunction[Log]}},
     FrameTicks -> {{{#, Superscript[10, Log10@#]} & /@ ({10^0, 10^-5, 10^-10, 10^-15, 10^-20, 10^-25}), None},
       {{#, Superscript[10, Log10@#]} & /@ ({10^3, 10^2, 10^1}), None}},
 	GridLines -> {Drop[Flatten[Table[Tlist[n],{n,{10,100,1000}}]],-8],Table[10^(n),{n,1,-33,-1}]},
 	GridLinesStyle -> Directive[Line, Lighter[Gray,.8]]]


(* ::Subsection:: *)
(*Number density of the electron - positron gas*)


(* ::Text:: *)
(*In this section, we explore the number (pair) density of electron and positrons in the early universe. An explanation for why increased temperature leads paradoxically to higher magnetization is the overwhelming growth of pairs within the hot dense electron-positron plasma. To prove that however, let's obtain the numerical density and plot the magnetization per electron-positron pair.*)


(* ::Subsection:: *)
(*Energy density of the gas*)


Clear[B]


lnZx[x_,m_,B_]:=x^(3)((m/x)^2 BesselK[2, m/x]+(1/2)(B/x^(2))(m/x)BesselK[1,m/x]+(1/12)(B^(2)/x^(4))BesselK[0,m/x]);


x^(2)D[lnZx[x,m,B], x]/lnZx[x,m,B]//FullSimplify


Series[x^(2)D[lnZx[x,m,B], x]/lnZx[x,m,B],{B,0,1}]//FullSimplify
