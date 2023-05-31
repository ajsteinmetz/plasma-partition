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
(*Define the spin potential and spin fugacity.*)


\[Xi][T_,s_,\[Eta]_]:=Exp[s \[Eta]/T];


(* ::Text:: *)
(*Define the partition function.*)


lnZ[V_, u_, T_, b_, s_, g_,\[Eta]_] := (T^(3) V/(2 Pi)^(2)) (2 Cosh[u/T])\[Xi][T,s,\[Eta]](x[T, b, s, g]^(2) BesselK[2, x[T, b, s, g]] + \
              b x[T, b, s, g] BesselK[1, x[T, b, s, g]]/2 + b^(2) BesselK[0, x[T, b, s, g]]/12)


lnZ[V,u,T,b,-1,2,0]+lnZ[V,u,T,b,1,2,0]//TraditionalForm
lnZ[V,u,T,b,-1,2,0]+lnZ[V,u,T,b,1,2,0]//FullSimplify//TraditionalForm


(* ::Subsection:: *)
(*Zero B-field limit test*)


(T/V)(q/T^(2))D[lnZ[V,u,T,b,1,2,0],b]/.b-> 0/.Sqrt[m^2/T^2]->m/T//FullSimplify//TraditionalForm;
(T/V)(q/T^(2))D[lnZ[V,u,T,b,-1,2,0],b]/.b-> 0/.Sqrt[m^2/T^2]->m/T//FullSimplify//TraditionalForm;


(* ::Text:: *)
(*Check that b=0 limit magnetization is zero.*)


Clear[\[Eta]];
(T/V)(q/T^(2))(D[lnZ[V,u,T,b,1,2,\[Eta]],b]+D[lnZ[V,u,T,b,-1,2,\[Eta]],b])/.b-> 0/.Sqrt[m^2/T^2]->m/T//FullSimplify//TraditionalForm


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


B[T_, b_] := (T/511)^2 b
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
f[T_, b_,eta_] = (T/V) (q/T^(2)) (D[lnZ[V, u, T, b, 1, g, eta], b] + D[lnZ[V, u, T, b, -1, g, eta], b])/Cosh[u/T] /. q -> (4 Pi/137)/m^(2) /. m -> 511;
bValues = {10^(-3), 0};
Manipulate[Plot[Evaluate[{Bc f[T, #,eta],Bc B[T,#]} & /@ bValues], {T, 10, 2000},
 	PlotRange -> {{10,2000},{10^(-29),10^(14)}},
 	Frame -> True,
 	FrameLabel -> {"T [keV]", "B [G]"},
 	PlotStyle -> {Blue,Directive[Dashed,Blue],Red,Directive[Dashed, Red]},
 	PlotLegends -> Placed[LineLegend[{"M(\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-3\)]\))","H(\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-3\)]\))","M(\!\(\*SubscriptBox[\(b\), \(0\)]\)=0)","H(\!\(\*SubscriptBox[\(b\), \(0\)]\)=0)"}], {0.8, 0.3}],
 	Background -> White,
 	ScalingFunctions -> {"Log", "Log"},
 	GridLines -> Automatic, GridLinesStyle -> LightGray],{eta,0,10^(-3),10^(-4),Appearance->"Labeled"}]


(* ::Subsection:: *)
(*Number density of the electron - positron gas*)


(* ::Text:: *)
(*In this section, we explore the number (pair) density of electron and positrons in the early universe. An explanation for why increased temperature leads paradoxically to higher magnetization is the overwhelming growth of pairs within the hot dense electron-positron plasma. To prove that however, let's obtain the numerical density and plot the magnetization per electron-positron pair.*)
