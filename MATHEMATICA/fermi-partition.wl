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
loT = 10; hiT = 1000;
Tlist[n_]:=Table[x,{x,n,9n,n}];
reversal[T_]:=-T+ loT + hiT;
bValues = {0};
Plot[Evaluate[chempot[T, #,g,\[Eta]]/.m->me & /@ bValues], {T, loT, hiT},
 	PlotRange -> {{10,1000},{10^(2),10^(-12)}},
 	Frame -> True,
 	FrameLabel -> {"T [keV]", "\[Mu]/T"},
 	PlotStyle -> {Black,Directive[Dashed,Black],Directive[Dashed,Black],Directive[Dashed,Black],Directive[Dashed,Black],Directive[Dashed,Black]},
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
(*Graphics[{{}, {}, Annotation[{GrayLevel[0], AbsoluteThickness[1.6], Opacity[1.], Line[CompressedData["*)
(*1:eJwV1mk4lI0aB3DzzPI8gxSyzChFioRsg6N4HgxSSWdUKqIUbVKZkDaDKOtB*)
(*KkuJ15Jsr3REUXPL8toyylaRpUSEBk1piON8uK//9buu/4f/x1vN8yzHC5GQ*)
(*kIhZuv9n1RGrW1GZKmAVIRre0cICR4aCs3KGCkBoMsdqyX1tX+Vy01XAv3xr*)
(*ucmSF/H4xOo0FZB2DjyttmRr1f6E2dsqsCfaPPhnMwsaPlyJ945SAalU5u4H*)
(*S+5wroi15qpA2Ehj4LcmFkzYbr4pZqvA79GKDz6NLFBXt92VOcCEjbydHZV1*)
(*LLi1nxowGcoEA/24jNxqFoQ1pGjIaTIhyKDu20wVCwL092X71DCAMpMRNvOU*)
(*BTnTVeFnvBngdI/vc6iEBe8Kd5jlkBiQ2KU7JX7IgqiM4vNP0pWBfP+nwCGd*)
(*BRcnK6VPWCtDKVIR+fIWC04ZT1107lMC6oi06pkIFiR3tNl94SlBbNzEndlA*)
(*FrCLr6c5rlcC3om2PdnHWfAj+HD1eLUiDOZhW/ScWTBQMLzN0EsRUkbrPXUs*)
(*WNAmSMosllCETnzTg3I1FljqKoWQ0hXA8k1EZg7Cgq9DNq/EVgrwPq+WeXPQ*)
(*GO4kJXdo960EFSOFrqgKY+Cl391GClkJsgFuHFGMMfhmvIyVWbUSDpka7t3i*)
(*Ygw9Pc8OTFfIQ2Cd0yV1FWO4ft7jxhU3eeB5hm6DHiPQX9SJsJiVgxqrHbXr*)
(*44xAeSR7ODtVDphuQ83nthpBMDb1+pOlHDB0bvOqPhnCmJpXhEuPLHi6M11v*)
(*+htCwu2wDL9rsjDYSf2+k2wIW0lN0bUMWbjKVxc/uGgAXrlDf8afroC/Tc49*)
(*etevD005j95UHlwBaVfYb4SEPnib606ETiwHbtDmA5WkzVDy0VX8wXA5BLms*)
(*mnP004VfM4YKA5dlwKCzdk8ubAI78YG3nLJlwG52qLZ6uhGWlXVy5OaloeRx*)
(*xbR3iib0emT640bSYDSajvXJrId1/GrXZ4FSkFZov515Rh2OJ0Q5DfAlITpe*)
(*yzdAaw2UOmd9uTBLB/n0VnsLZ2Xwb1qB77Okw16BwWSIrSwQBq3mORcwKP5L*)
(*sik3hwQPG2a0Y89jQKd9jhB6kEDmMEPF3xcDBxZzVJ1Jgt64Y3PsExiQXdZg*)
(*NnESEDQurhxyxWAqqzljdmyBX/pQ00LDBgObDZZpGTwxX0M1BM+SxYDrTzOS*)
(*uT/Fjy7L3Rwtg0EjqWdIfvkUf3pnyxquFAY0bG+3RYiQD5eVJGyoGCj8Ix20*)
(*uWWSf+h9EXyaRWE2cZHh5DXGv5P0wUp9AIWr4wyZ+apBPiZlzM4sRmGx4U+B*)
(*puojvpXI+/jLAhS+qCl4PLG5w7/UnxLVk4cCO6QjNrHQHx9/stimkIXCW+2D*)
(*+XIzRXibW/OhqLsoODSb3nLPrcWTiz0DuTwU7rzb3mgX1YW/Sb6dGn8VBUEc*)
(*3c1R3I1LhjW8KLqEQqVDxjMnr/f4FRc9ytcLKDBsu/18dXrxwxLieLeTKKCu*)
(*VsrGMIBrcRLybTkoiDx7WwNzv+BHtta+9nRCQXWZ179z/jWMp274JQzeiUJv*)
(*8K41bU3DuPScm+lzOxTsj2seYY+M4MIsrVq9LUv7JHX2l9HH8AoR/6OSBgpO*)
(*bGuzT/MT+FT/tARLDYUV84YvdlyYxLWb1mtwVFHQP1dfbjM2id9Ljz4Vo4yC*)
(*XZFnwMKb7zjP3uWXxDIUfiBlLyoahPi2lMkVYyIaDAiJTdqnp3Grmvyxhmka*)
(*nF/5rUA7cxo3n/CuffidBumsEEF81zSuY9Uf6DVKA9/fGy6n4TO47Kigb6CX*)
(*BrbchlSy1A/8g9njwu4aGoheG6dnRorwds8zEU+BBpx7B355lovwlpiNh2+/*)
(*oMFjdM786ZAIfzmQKe9cToM451hFFcufeNbNxEut+TTo2j6+oDj+E/d5x91W*)
(*l0ADN4ll1lyTWdwL0VfPjqNBJFM728x9FnfXGZ8LjaaBUiV2pDB8Ft/NO/Y3*)
(*EU6DsfFLN7LbZ3GW1l7Fqos0aLQ9OzR/+je+EGgyVOpBg4ma2jr3/4jxeKXf*)
(*wRm6NFgo3jIgnzePuya9ivmvNg0YH3uKKurm8Q2yMSkNmkv9Or624ed5vEpy*)
(*zROhGg3+YO3eOav/4F//sIcJRRoU9KyqS4v/gxND8TsHF6jwSPhxhenZBVxY*)
(*osVQb6NCJzL2TpkuQezevr80y48KqWOKpeohJGLkxH5u/VkqtIotOEPRJOLa*)
(*jf3Goz5UiD36IyrtDokorN3/VO84FbzZy9d/KyARGH7g+TNXKkSKksNku0gE*)
(*3+jgKwGbCiWDK/XjtBBCb7Xb2zlFKoiCHRPb6xFimdBj2rmSAk4VAhPzcTJx*)
(*bNfoSbMKCvhpoEWlIjJRVeg3uKqMAgsopX7dIpk4dTK87UsxBRzM3BlCWQpR*)
(*/ym/OPAvCvTWbGoxNKUQ1zpFp9IiKZC72lnZh0chhM+jP39yoYDjut6XFBkq*)
(*8Ta8vOP8DzLs/XFQZq0yjVjedvRSrzkZSizWJYi7UUIuOk3wwJQMvx1td3D6*)
(*UULRrl3jqDEZJFrOuuYOo8SqKmvBqC4ZxNcOJtqIUEIrT03j11oy5F/PG+HI*)
(*YYQ1r/+1LEqG8CT7i407MCJA31XNvh2BkKqj102qMCLoW2KAZBsCGUqWxuo1*)
(*GHElt6n5dQsCuvXzzZJNGBG2yjxgTz0CP89r17ztxoh4jNF85BkC4XaVlpYz*)
(*GJE/0HXhygMEvsTdqI/YSCeK0mSa8HsISL36WmetTydK9tmtIacgsPvDz9E5*)
(*EzpR3lLWGJmIwFo5t2IvNp2orUhSvRuOQFrrkyglDzrxD7eFezAUAccwebkX*)
(*XnSiSY/SuDoYAaZg8b2HD514k83lZl9EgOI3a3Y/iE50eBQ0HPdHQCAq7TLl*)
(*0Yl3zM+rN/khS//IvleCCDrR08nkTvoikB26afFYLJ3oi+c0PD6NwBbehfhf*)
(*t+jE/wBaakZS*)
(*"]]}, "Charting`Private`Tag$108547#1"], Annotation[{GrayLevel[0], AbsoluteThickness[1.6], Opacity[1.], Dashing[{Small, Small}], Line[CompressedData["*)
(*1:eJwV1Hk4lPsXAPAx77xTIstkG2u20mKfUcnP9yXkIhWShptUVFIKUanGkoqh*)
(*CFFKilEhSbI0bnPIFZN9KcmWKArZieL6/XGe83ye85xznvPPUT3o6+BJJpFI*)
(*nKX4fy71MIuPeqAAjdX3hONrmbCDLu0ol6YA1n7lHpFL7moYoGWmKoASNf4F*)
(*e8mLKPZmWYoCTBzWsvBZsrlyd9xsogJc5BXKWS65qv1CrFeUAux7NMqZqmFC*)
(*i2NxjLm/AsyVdV3aveRhS91rcxYKoPv6+rNJARPU1CztH/TIw90jB2yT3zIh*)
(*3gUPHAmTB1lWIbv3DRPCq25r0NbKw+Jrk+r7fCYE6jln+LyhA5Mzn6f/ignc*)
(*8dKIE150wF36lIkCJrTl2G7mCtHB56BmYUEOE6LSck+/SJUDXaHxBCyDCWdH*)
(*eKJHzeUguZT7V9IdJngzxs46dsnCcKs5s/Q6E5JbGqz6Q2RhQ5us5FAoEyxy*)
(*L6fs0JQFue/VRLg/EybZB8qGymRgzugqI9qDCT3ZX60NPGVgq+355072TGio*)
(*T3iQS5KBiCApLtuICabasqFCqdIgGDNrKVBgwkDftvI5M2kw3fqLvXeRAbcS*)
(*klvWd0nBx6jomuAOBoSkJlkLhUrBcpWAXflFDDiZ9jpGTFEK0g0XiIIoBnz6*)
(*VLJvvHgVPK/X4wpcGHD5tPvVC26r4JG92J44DQboLW688r9ZGnQ91PRU6zcE*)
(*uW8ZXzPu0MBwfkdcx0NDYC8fq+01pQGvwnDTXZYhfFf1vLL3kyScc0p6Z4IZ*)
(*QlxieJrfJUkg8Rg6Bx8bgImQgFNBl4TKZc30TGMD8Mzs+zNUKAEfZ0bo28r1*)
(*QcB90shjSUDZasuNPpb64GWsPRw2LA5KvJ+tzbl6kNfpOtduIA79IZn3j9H0*)
(*YGbCQLonWAy87hVRrBp1wGpuX5PDy5XQxNQ4QNusDStftjrQfouCKdRdxWI3*)
(*QIf7gzPIUBRW/JHI+CtuHajzy1xLgkSgUz2HFcJeC0fionb28FfAEXRd+ayZ*)
(*JuQ7pvcHzApDnRPpTugNNTgjkEDOpsIgzdgTrxWtAoR+nTE3YDnsvbZDZrFd*)
(*DpaLMCwe5C6DAuHAAdEbkmB9e0Ti+xQVrDoMA8KkKWD2Jut71TgVDuxWvKzb*)
(*iIHxsFfFo59U2KIu0lcejcFGs+4gz0Eq2Fx6uuUehoHkYH1XTwcVAmjTXM1R*)
(*IWjf/DznwxsqFCcrqYdHLvJ92vyt/42jAs+Ec9KDP833JOupZVyngomURsek*)
(*/jR//8ah+TAOFZ687NYLyZji7wo5/IyIoEII90VndOQkn6m1R6b0LBVOaHWu*)
(*HHUc5y8EGfXluy/N8yvIz+34wY+V/cVO06aCysO2KK+SVr5rQnl0wXoqyMym*)
(*3ulqbuavkYy+XbWWCrUFrxafjTXyS1eovBhVpUJdSI7mB98a/sAfi6+EDBWO*)
(*ah5sfFvN4xN9sXafF3DgKK2YCFTkIREP1r6peRym5+MvjjcCet+p7iX8C4fw*)
(*+sg9lSIVyOdDIVt/Aod+lUPywR4ClCzoyA/9hoNyjbcZyaYZjeZp0dUacLjZ*)
(*nn0nJroL8bTHNY1qcRDWbYla6OpGV7J4BjYCHF4HJsbsv9WDFDPs7fwqcNAY*)
(*u52r1/kZWScFsMuLcWh3220hvqsPpV2E/oMPcGhKtNQpbhpAx+cixwNTcZCZ*)
(*SZgMZg0ioyDHxagUHDbqzNUH9QyimlNf5V4k4nDDh1Fo8u07mj0kaodF4XDg*)
(*by/HmZ4htMvGJT/dD4dNP4f5dlaj6NtRF/9KX3zpb2Tbd18YRZeuujAGfXCg*)
(*0b+0teaPopwKl0KdIzjE97mUzCuOoeVo36sSVxx27s/3D/k+hviGrPJ6Cxyy*)
(*SV9Ce3wmkLMDK3zcbOkezwVbsZQJNHSKZSGNcIip57DeVU0g+VxWJWsLDsm8*)
(*lt55tUkUqOUq6NfGIc9Tb8K3aRLpKLk1zcss9QcK+VmqTKN/t7rFK0vh8O5H*)
(*ycmd26eRG8vNyUwSB2MBl4P7TqOoJLf3V0RwMH3s/S34n2n0VfLvdkkSDgmL*)
(*2leLnGdQKnV/79pBCjzBDKzDg2fRylH3cUceBaD4/GJz4hw6bD94bHMxBf7Y*)
(*7RJzKJ5DpTl+nxVfUuAAVdkmoX0OeR+LaOjPpcDrNVEfjirPo8rerNyghxRY*)
(*uYqx2jt9Hl1qnfJOiaTAqJXE6Yvpv9HoK86X3r0UGGs3dU0LWUDb6dKub50o*)
(*kH4NxMTvLaDUoNSm7N1L9cvR+N6SBWTHeF4WYEsBDeuGlICxBfT46fv7VLS0*)
(*f/KU6IT7IvJIU3Vbv4YC6lWyddoKJKIpoqjl9CQGkP/Mj9dPIuaN7NmcMQzE*)
(*IIRDHyMRGgN967gjGEwFBbQd/00iAm1o7LYBDNpkzEOHaUKEvPiJdagTg51+*)
(*6fSjSIjwSFa7JFqJgVTy7/Ddt4SI0ayYtY+SMHCeeh1/14RMiDccOt9hjEHq*)
(*BxfjY74YQeOk1N/fhIF744ePcecwQsaqWeMQAwOur37v83CMUCw1rx/UxqDO*)
(*g3GmLQkjtB6rasysxuDz51fXy/gYYR7SXSu5DAPDbWNGNDEKEajnqrq9mQyR*)
(*2nXdclwKce7HzcAVDWRw6wnL25FLIS5kCt7V1pBhTTS1NriIQoQrGgc6VZLB*)
(*OGuzQkU1hYhdTn/nUUKGQ1y6B3WEQmT1vA+4cJ8M8E9GvSsTJ56miAnQXTLc*)
(*urVqZL0pTuQ5W6lgt8kgVVsjMWmFE0U1L6sjb5LBNU4nJsAFJyqKE5STIsiw*)
(*f79xmG4wTrz1r/FnhZHh/UnbvvpwnBDoUKqV2GSwuCN11zsaJxoz/P0zzpLB*)
(*Tuu+4Y17ONHinl115AwZsrw4a+iZONEm/0Vpgx8ZfstfSbyXixOfWuX9R06S*)
(*YY/wYKRiEU50xTpUPT9OhqLMadFbfJz4D2L2R9E=*)
(*"]]}, "Charting`Private`Tag$108547#1"], Annotation[{GrayLevel[0], AbsoluteThickness[1.6], Opacity[1.], Dashing[{Small, Small}], Line[CompressedData["*)
(*1:eJwVzXk8lfkCx3Ge59jKeqznHJW9EmM7T+NVjd8kGdmSvTCNijYyCHVtpdSU*)
(*el1bpTSiMKROpRRRz9c1I07kULpiQqJRliKytVz3j+/r83r/9dXfFuEZQklJ*)
(*SUXM7/+tCV6bdbJAAGkyMRPcysCNp+mlky/AxyRX5YB5d0sGucV5AjhtW7fY*)
(*e97fSHpmba4Agq+UleO87Rf3ZEyfEWD0doq+6bwbOhPSQ08KoBm+wWxMwuCZ*)
(*V+Vp+2gBxlq7n8fPe2S9xW+zDgJY50cNpbQwMDBY717Qy0ek//fyes0Msvxl*)
(*YkdT+DhRhaQUMYMjDeeNuEv56MzYtff4Iwaxlr6FYXU8PFj53qqujkHReE1q*)
(*eCgPBkN0znWWQcc1F9siaR6UCtVfpt1ncDJfFHk7TweS7d5bbt9hcGC0WnGX*)
(*vQ7u3mxPiRAx2CMcO+DVrQ2p5OUaxcUMcp5JHAcOaeO6skjim8fAQXQ0181Y*)
(*G7P2fn2CMwwmkn+pHa7VwtfoDlP2BIPesjdO1iFa8CnnqT9MZCBpyS4QSWnB*)
(*0SNFTy6CgZ259mHpPE1UGr7QKf2ZwWD/uv/MrtXEr4qHvxx0Y3A2O+eZabcG*)
(*dvHDVyfYMjiUd85J+rAG3k6W+O40YLAv/+FpZV0NmNjWc0PkGHR1VW0er1RH*)
(*9dLWTZ5vhTgaufV4QqA6DMWL7l5uFMLym9mxH6a5sNVwysVlIXT+KXxTeIGL*)
(*udYu5/MHhUiWH2vus+OiMFHorOkqxDv9kGN+XWqoPIvPRhpCZJw5kh+VpAYV*)
(*h02LvnbYYI20OO1PnhqqXHY1+qbbIKS4/8vwXVWkZ05ZbCQ2EBeVtlZvUUV7*)
(*opud0jtrhK4yH0kZUcGvvgN3/nvYGjdfBsx2WqtAKv7zj/4a1pj6aK3ZG68M*)
(*t1d2c5svWsFxdnObZ4USOuOSLqcusIJSRbsn97MiopWKznpEWeLvrQUxxEYR*)
(*pUFha+ujLWDI1gZUxS2EWGGFtMVHc+zMOLmxl12AtJ6pzEoPM5R7XRnYP62A*)
(*Bsq5osvPFDFiVeJrp4Aa61W5Atdl+NHqyaqi/fKIaBx4EtZuDPmFQocCkRxG*)
(*F3GS3xoawun8qOq7SVmcOnpDLytAD+naM8n55rJIM6Nsw+X58HD2L78SJYPO*)
(*5zuCVuipQ+nD1nGvag7Gik+oa/rKYYf72922lRzYiySBK7lyqLkW9Uq3goOm*)
(*66Z/rXkiiz27UyUDIg4GTZwMZ9bLor7vqijuMgfHht0arYQySGqf3JN7ggPj*)
(*Wk71vYU0PtxPe93nx0F38W9ucx5f2J94mgGPvDk4klNvsrzuM5sXl9dWtomD*)
(*xK7KtTbMZ9ZVeKt2vwsHxx/GNfTpzLEl159fkiUc6L6MWNbfPc0G5+sHmppw*)
(*oEjlvP995wTblnrvWeQEDfaN5cEY7j/s3Er35LQxGvU9PleFym9Yo8H+5UWj*)
(*NLalbZlpWjDAxjpzkzsGaSy0SdR9LPua5auELycvacDs+IN67W42OMcgSbGe*)
(*xo3ye1x6TzN7ckPVMpM6GqWdQ+/P3mxky2c3PiWg4Xj6Q7nCt79YTlDisuj7*)
(*NPYsKfTQ3/2ALdHvaHshotFhtWv7UY1I9sPV00v/OEdDeK8kn0+LCS/QqA3Z*)
(*NMIzg8SGS5uJvVJ1QmcGDb07Zc6RehKSHTHYqnSKRqLMa16N6VNiK1yXsP8Q*)
(*jaw7lxSfSL8ghx5MS9bupqEi2fu8rvcVKd337/iAUBqZcro3NVP6SNsSE5OY*)
(*7TSU1SzvGOm/JkYpXvElQTR2O35neMm/nzQ4ioxVPGlYXl0xpCp6Q1Qk2//1*)
(*9yoartJloZFz7wg3Lbfl0vfz/xcr45RODREtx6dG24U0+mVcAnL4w0S3xr7l*)
(*rTmNsjWnXK5ZjZBlJfpGU3o0Tn9rnONveE/sD/U0q8nR8FmY2+SdNUbWr9Y2*)
(*bKdpyLm7G7HNY8Tpk/uBHCkaP1zU8omRGycbwx4aLJmlIPZz6Gk/ME4C/X+P*)
(*MxuhcLgndUe710cSaxmg/9NTCqa1Yw19IxPk4FBm7AIJBe0Lnq51iydJQrH4*)
(*cXMTBR+LW++MN06SI7qrYr3rKYz7Le6NF02SdHne4+AqConOopgbOz+Rq73P*)
(*9ydcohAruFBhVD9Frucqi8lFChYD2UdXf5giN30dl9DnKeSGtl34xJsm95oq*)
(*Gk9kUhjuyFj/S9g0+bMye/G5VAr0/Ykmc6UZ8ii6KXpLCgXdxIhtWcIZIv6O*)
(*07gomYKDcDCpNGCGtBZGRxceoGAXWlgxXDJDnm0ta9gZQyHdoH+fbssM6eC/*)
(*XrQiat4NV0akJ2dIVzs/enQfhRo7k/A8/izpTvdsuLWXQvEXBf2PZJb8D/EP*)
(*alc=*)
(*"]]}, "Charting`Private`Tag$108547#1"], Annotation[{GrayLevel[0], AbsoluteThickness[1.6], Opacity[1.], Dashing[{Small, Small}], Line[CompressedData["*)
(*1:eJwV1Xk41dkfB/Dca6Ia2/fe773uRZNoKkaJe0qyFLn5mRZPGoRJJaTNmkrZ*)
(*U0KFVJbKrixJTSVDeosoLWQp24jKtBNaJDHn98d5Ps/rOc9zzudznvN5Ppqb*)
(*vde6cyZNmnScrv/Hik3LTsRkqsE88s7j3KcEq0SsnWqGGlaFazlkU3c3vmby*)
(*zqkhcCrpzKCeMI9PrEpTg8dTn2dp1BbTnyWMnFSDlqJZWwL13Y4D8R4xasg6*)
(*VRMTTN1id+Oohb8a3nfof7Sj/mA1P3p0uRq8zDWyxp4QzJxptTqzR4zDKpXJ*)
(*utQnHH8K7I8QQ7HqmOvqVoLIuynazGwxXFbKOjm3EATq2+fsqBbhu1VRqaSZ*)
(*IHeoImqnhwh6cRu0uU0EbUW/G+XK0P0ol6knGwliMop9/zqnihVjc2LkGwj2*)
(*9pf/vNVCFdVfnZP9HxJskwzutesWYlXJms+J9wmSWxqlfWFCiFdV7827R7C8*)
(*+GDaqllC5Hbz40/UEXwK3Vj1vkqAy3E53ip3CHoK/7U2cBfAMrjNY/ttgsaG*)
(*pMziSQJcujYzRO0WgbSd+aA1wSJIqOGmQF3xPN4o5QeL6ybXzSYqCc5/jmuI*)
(*HGXhmZjS/pw6RC3qh+MnFvdUhs/kU//muXs99zWLyJ/MbQn14fE/VBwbWKh7*)
(*zhZZ3SQw0xOGy5xjIZZWxFiWE6j3qHkknGGR7TreZUj9PXHG75ppLAayt2lr*)
(*U98YmcsuS2bRsTEpUZbasMYkPyyRhbLXh7c1fxPMcd78eOIQi9y/U8MtqXlH*)
(*LmqOe9P8o78NLyojeP3S8vboMhat7SLJ9+s0/0VHI7YvZXFe4t7xL7VHzBOL*)
(*LjMWpovNfJqoX8z3qr65hN6/Ye22fOruoGM1YYTFr41Kp+ypW5TbayfPZVGY*)
(*wFZcvEZwy2TXfUaZxXxf+arVVwlOJSW36HTzMb14JMDhMkFl/NZ517v4KHfp*)
(*YJdS98UZRS/r5OPVD52cOdSSqLYljm181FY2R38rIWjerZod1cTHRS9nn1Rq*)
(*ZYdk355aPrzMei+0XSKIFSUrnrrER+7NYAtpMUHYudPWMuF83F4dNPV5AUG0*)
(*cVHwzlA+3i715lVTxz/BlfZgPiabq8vkUKcrvlW/EsQH4o1i3KkrQ0wGNwfw*)
(*wbtbbPsqn2BsQ29KjScf2X2yyX0XCPb8ovvu8Go+pt1UWv8gj2BXRuVRRXU+*)
(*3uR/qjDJInA52q8eIeYjPmvoiYDaJmh60WdVPiySxI0fMwlmrQup72L5GJ+x*)
(*3zWHukvOTK5QiY9+39CBKdT/874ZtoLLh63Ur+JROoG2WYV/2FsejPeHSM3O*)
(*EnR2lq0fusHDxtT57cWnCQz8fCLCSnlYluYiF0kdLT+7UPE6DwLTUsaBeiFJ*)
(*GtP5i4esPScwforWf2xXuttFHoTxRoErqaXLtPpaMnl4IMk07E0iKMmL8ymN*)
(*4SG1tMB6IIHgoK/r4QMuPPxjKLJNjiUoP+8U+86Jh6HgmVEbqIf++eO403oe*)
(*7o/giDa1q83K00b2PLwYWChTEkNgrLX4/Kc1PLjdss+tOUIw0MLU7bDkoWtJ*)
(*Tc7LwwRORrWT/9ThoTdAmjAeSaA/8dsh0xEGjkXN0Yn7aT92F2vVfWEQot9l*)
(*Y0vtXalfZfuZwR8StU4F6uJgyZjbEIPT++rco4MI9MaW+MW8Z0AKs2bs30eg*)
(*881mw9Meel5gKuOwh77n8NaFvvcYSN91Cl/7Eai+yvk3J5XB0KxY3lYvgnNC*)
(*i26SwiBo8K68mFrL+llr7WkGx+r6a+9vJZiXL7rzOonBaE7w2XnUVtuPZf92*)
(*nEG7WUDJoAeB38CejVcjGeyQVVm3cwvBg682HdU7GOgMqVcvdiUIlR98+NyM*)
(*gfN23asz1hG8bGvKWGfKoCd0qmqFHf0v+Vf9a5cwSH9nbOBALbDZKyowYtBg*)
(*OhQdt5agKHZii68Bg+WjXuXDtgQdiopjE7MYXESbytVVtH9YXR0NBQa16ud2*)
(*fpcSvNV0P+TQqYKx9LPyZ4wIEk5GZviFqCC0V27RYjGBiUx9bI1IBapfLmVH*)
(*DEvgnvfyx/vryuCtkP08XCtBfW7+43InZSRxtX90n5LAw1jvQ8QHJcgy5y2a*)
(*nCQo+cd5tMNACellJW5xYgm+DhuwPfsVEdSsfEi+2RDS0fVNa68pIK9nsfu7*)
(*IEMoXGtdy4z9jNiTf35s1TBEl2vmbnPDn7Fgd0FQa6EBtG5VOZftmQbeI80D*)
(*rYYG8EyIWdNzayrCq/RSJSULcMUuuy9gZAoCld/oLVRegN31yub2ZlPgRMb3*)
(*TffSx9IFj4xzA+TxLNtqYIvpfMhPkyzPLJbD4m3vp60L04N1Sr/y28+T4Tj9*)
(*6kqHK7qIF34LzdCbjJc+gjBJ1FzY2jheyfb7CU8jyuYLLGdD4aPrkF25LHS3*)
(*OvRdaNRGU1Rpi+8nLvS+hpsXPNOEUqNbUJcxF13bb5j2F2iAiU1rSF/EhfVg*)
(*wRXPYxoQSJu13SRc1GVzi1p8NaBeYdHwRo+LXWdThEeMNDDngqb21xlcbCE1*)
(*ls531GER9uyhihwXCcNVARlddG7rO2uuaObgTePHoSMyYux7lxg4tZEDg0HD*)
(*Ea8XIhzIq7//8AEH+roHO03uiBCpbhy4rpYDwe0vv5ZHixAvL7q/qYyDzrI/*)
(*l4iURCjoeRJwIJ2Dk9Z2ufYiVVxMU6w3P8NByYZFcptHhSixl/7CTeFAQfPN*)
(*R9dOIUofXLt3JJGDougeK8OzQtTcSJp+OooDE/lel0UzhKjzf+DvFMHBlJWc*)
(*LU8mCVE/T/aeRigHVu6COs9eAR7n+Pvn7OXAx14udlOWAC2uhXc9d3PQOrmt*)
(*516EAG3iFxq6fhwsXRN5QNtNgM5WsX//Lg78foQ4+loK0B2/9u7l7RyM/W7u*)
(*XaIlwH+l62tc*)
(*"]]}, "Charting`Private`Tag$108547#1"], Annotation[{GrayLevel[0], AbsoluteThickness[1.6], Opacity[1.], Dashing[{Small, Small}], Line[CompressedData["*)
(*1:eJwV13k8lN0XAHBGC2HM80zIDF6yvpZmyEy2bAkp8bNlCQmllGxFb5RsFYVI*)
(*KRHRhihLqdAhoRKKIlv2CkmoVNTv+Gs+3898nue595x77zlXducBWx8KDw+P*)
(*Ji8Pz+JvpadxanwOEx5q37pl1M8BKwlRu1XZTLicHKxqiu5r/Uhez2JCwqmp*)
(*YnP0X8PklJoMJpTorqjZijaRfn92Lo0JnoMzfG7oxq7w5F3xTLA86k79D91u*)
(*V3HGJJgJne61G+6hP29knfxlygSdaKfr6wY44J5pLrdnAxOyDM4b6aObZz2q*)
(*O4yZ4PItq88IfSc3ebbMgAl/JiQVLdEhvNM7DmgzYSFJcmg7+ldlmc6wKhN8*)
(*KyvFo9H8WrrjL0kmjF7eH/kavXr1xq05/Qyo2gLzEYMcmF89HKXyngGl6+qN*)
(*o9Fv5aLvl/Yy4Ci3Mf4k+pRCjUx9FwMqrV4op6KnlNfPjLUzINGIkXYTXc3i*)
(*pGs9Y8CLoyYa7WiX9QqDDXcZEDOil7BmiAOpTksPTUYxQDNPbNcU2jb/8ZJj*)
(*xxng2r07+zuanD+cKhLJANFdBu/n0SlXJu+wIxjgMxa3V2CYA2dHO8aDQhnw*)
(*lvH2/mp0Uki+5w8/BmSozbIc0QnJ1lY8DgxouZ2TD+joxovypBIDdvd0BuSN*)
(*4HwfpidGKzDAiJ0eW4iuLbgwNyvHgPVKCdllaEpSWlOHDAMcp5if6tBRDmdD*)
(*MhkM0H3l9HQEHTl4sk6ZyoBj47W3lUc5ED4f6m34TQKW1ha7lKIPsR3z9j2R*)
(*gHcsrmDvBw4ce9F9Y65GAjrFaIqj6JO7PAtiQAI4XK75F/TFy/tKLldJwIGh*)
(*kMu8HzlQyR8NTfclQPU/yVRFNGWguEetUALOGDPKAtGJZwVEP5+TgFrplb1C*)
(*nzhwbboydv8uCZj/rpq4bYwDnYWbta/xSkB9Ktl/8jMH4rOLAkuzVoFBTfg3*)
(*g68cCJt8JORrsgpsTt84Z/aNA3u1vobZ9YnD5YLfvda/OJDe3mo2EikOrunH*)
(*5/P/csC0KCbDSkEcXp85uDR+KRdmj+2omagRg3+n7z5NE+JCf8GohaaPGHj0*)
(*hDRX07nQ2nIup4hHDKIixEqFJblgoC5+nDdLFLQG7A1z5LnwcXhD7S9jUQgt*)
(*yh54p8aF8+fS21X6VsJIWEZbC4cLkVkXLHiPr4SdrrckKAZc8M+uPkOVXAla*)
(*lKiAYDMudHc/cJ6uoIO01MnopVu5oBkUEBV5nw4j22yEZdAn+ZUKqPfocDZW*)
(*+pQumss5N69SSodHesmOB9DJif5XvG7TYdpe1rwTbWYsN9KeQ4d3L+8H3rDm*)
(*wp3rpwPux9OhW2d7lO7/uBAT6HEifDsdRKaCddn2XHh0wyVh3IUOqfevl1qg*)
(*p3sdklyc6VCtOiHuifaw3HJB25EOeq47cs6ideV0bsxa02FIv778K/pLO9mw*)
(*bwMdygI2tt124IKLdv0yNxU6FI9lOYpu48JZ/5oVTcp06DQdEFRFN+ZVUvWU*)
(*6NA6XnnLCM0hSsUk5OnQHDaV5YemjV1RfCNFB7Fkxn5AP804bGZF4Pe2Lznh*)
(*48QF9l+1uPVzJBR8IxUynTE/fUVyDd9JOFOhklCEPlDNrrH5RsKl1KKux+ii*)
(*CK15r2kS5j5WcAbR6vN6QfETJBwhhLoUXLig8tPSvaOfhNSHZi9vouVnfLmB*)
(*z0hwXR89dtmVC3Wvx9p+NZBAo7DcC9BeJfsCY+pJINdyHzxAXw0MKLzwhIRa*)
(*s5UKb9GyU6Grq6tIWFGt/1p4OxekP8dRBUtIMD6RFxyGXvUhbzTvEgkJj8+4*)
(*6bpxIUvcpI9zkYSeS7+8NqLlLN6/qb9AgueLMhsb9JpbEk8/niOh3O9Vjw96*)
(*o19irloSCXoTE1VJ6KAvoTvKokn4r7A4rQ8994+o08YoEg6Wh/z8gI6wKbF+*)
(*G0kCi3rC5Cv65N0Jg7kIEhq228XwuXPhSrCnlH4YCZL5DXWK6KYfll1P9pEg*)
(*P8aj6ou2Vf74yt6PBN2tQj0H0B1Osc9G9pAgemVnYCh68MHjiuW7SagRK7aI*)
(*Q/88onVhsycJobo8N6+ijxa+SuzxIEFiJuPwLfSSXv+4/e4kiIfwsO6gaQa3*)
(*Dia5ktAYP2hahVbkkbZvcyChPs2Y+RZdyH602duehLNrXAR60BqeThu+2WJ+*)
(*oaJrAK1fm6IpbkPCfa1lSp8XxxuznHC1JOHCUTFvigcXOsvy+CcsSPhHLHXv*)
(*crTbiDFPhDkJboqv7YTQvmbhX7JMSUh56NAsij7G//XloAEJG23zPJXQw52v*)
(*s+3Xk+Bbzm5XRVveKguu18P4rrwnx0aLWYZJ5Gtj/unbnLXR4QyXCeY6Euxt*)
(*LAz10QNjeo/PcPD5ml5eI3Rhwl/vQE0SAjdnS5qjye0D64bYJMz0XQ2xRIep*)
(*PVnhwCJh+eXZm1boDS/j7mirkdDsxL5uh76Z6Rudr0JC9p0v/o5oqr+lo+S/*)
(*JOgM3qQ7o7uo1Pm/CiT025ZOuaGN3n9pDpQn4Xv9O+Ud6GvFr3KGVpPQcTRP*)
(*fyd6RWRpiIMsCcMCf1S90QE2aeYN/5DAM9P63Qf9ViaUoSNNQiIsz9yN1v/q*)
(*9DlfkgSv2dPSe9BXa3RBkklCk+DmI3vRy1MkUxMlSDh2TK3UD71v5x8fnlUk*)
(*ZDxh1u9Dv9bs1w4SIyHHnK9sP1qbr1ZweCUJ+w+9DPdHZ7Xl9jnQSVh7z1vm*)
(*AHpJXuzdBgL3q29p1qL3hOyO0aGRIDaW+2PRzaabthVQSZBOlVENQGuJqqpI*)
(*CZMweY+us+hLI0ILiYIkSOXtll00z73JFp4VJITlLRtefH5XXOvVIH6c70JX*)
(*5KJfOJYcHF5GwvvPzT8Wx8NWOmfhuJQEn5Imi0Wf/3GQ2ciH+TleHbw4/vmG*)
(*bZM6FBKsExLDFue3M12npoCHhJEFxW2L82/0ZZ6T+kvAEZ4gkcX4qOss7Epa*)
(*IKC13DnHF50i8F6Hd56AfL1GocX4zr0DoeBfBLTk3/jfYvzd86++H54j4K/Y*)
(*4AEvdN1/MSWOPwj4kui52xOdxLRw0p0lIOlJx8B29Oz4v6qF0wSYxlz1dkG7*)
(*VAr+kfpKwD0nF9iGVnRryeWdJMBTSIdiiz6tfvdQ8AQB++vXftyK/rqQsmlk*)
(*jID2TU+ub0ZXZjl+afxAwGdli9um6NUHtGt1RwlI5RDfFtfzSUNGWuEwAYH7*)
(*xUTXo+36e3WTBwhwjwvv1UKPyfrEbesmIGyjxZnV6IqMl9ND7wjonW4ak0Kf*)
(*EOV6BHQS8OG0mIIEWl5gufapNwSINBazRNAeUzc/PWwhYLn2/KWfuL/X7CUc*)
(*zZsJGDLN75lBzw8drm1rIkDXLuHX4v6/1GGZMfGMgO/qmc8Xz4f26okt0nUE*)
(*BC0sPGhAX9V2eJBfS0Cjg0k3oANLqhTW1RDAKVHsfICmXk/8Y12N8bwW612A*)
(*3nSGfed4BQFuNeefnkZXu4bQRwoJ8DrwLN0CffpNz7HAAgJYGu5+hmhX640T*)
(*C7cwn9f3rOKi50zEnordIOBwxUVeOfRalYpDFjkEyNfZUH/jeXxz7ve7gjQC*)
(*Xntk8uShDwV5m2mfI8Bg/nlg+uJ5P9FUUpdCgKOjceVp9OD7zITeJALipSp7*)
(*QtCSDYbrReIJKOxq27QBfTYtOjvoKAEl3c+fdGK9MSpmf3wfTsD0mpUpz9CT*)
(*jb0sqyMEEH1hGg/RW+bXPVYKIyCgRmPoEnqZ10RvTyAB9cpFIi7ocLYD09yH*)
(*ABtz6dAWrG8qlhSvMi8CzjYf31yFfudVnC+7k4ACreAv+WjuBQHdeXcCFA0y*)
(*rsW6Lq6nKqe7TpgfAaNGbfSuF4rnmVsI8O+vaU3Deis60tZ70pIAyQx2SiS6*)
(*7k+kwncLAmS2pMv6oWU1u8taNxIwajlxyxDdk57UFmtIgDPjndso1vP/7Zqj*)
(*TWkSMC5jqK+C1ud9nlAngfkP9ZRKwP6h1S9FXWcVAavp9z4fQPu8dWm5LUaA*)
(*nNKdeHt0YsE4eYGO602KrSKFHnAQvrRHmIAIw8z5AkcuxOXb3KTy4n7bQnn3*)
(*GPubVruOOqePNGA78FQ9tuOCd/WVXS9HaXBE5e5UJvqnsi+/yQgNkta3fT6C*)
(*Xv1nbrPqIA2276eYrkMfvMloW+imgfnmNW4FtlxgLLgNXG2hgVBdvuop7M98*)
(*rg8vTNyjQcO6D0cZ2L/tfLhWJ7ycBnanDf1nsb/b0RwVIlhGg9wPeqrN6O3f*)
(*Zcb/vUuDLdqNwsfRtmZunbsKaFC1tzZo2Ar71dE3Jf1XaHC+7vHNnC3Yfyg2*)
(*7G4/SQMrwzilb5uwfumK5XqdoEH1Dt7sRjR9q0/fdCwN/IOshzLQ1EN8DkQ0*)
(*Dfb78j83Ri+pNzTZGkGDlWK+ogkWmF+fCsnGQBoIWlxKoptz4fm1W68eudCg*)
(*q1fB//cGXN/O6z+qOtMgIEn80RN0rPCrPxnbaDC/W7A+AW10cE413J4GzLBj*)
(*okz0fVPzOP2t+L+K0QWuCRfyhoZ1q4xpcDZXR9DTCPsJWZm8amUa2GyI647S*)
(*x/m/KX24RokG14RaC03QOqfMX2Up0KD8JVebgp6d8v9zdDUN+Eu7wyP1uOAH*)
(*VU6GkjTo+c338LAuF5w8XIVBhAYRync/uGvj+tRV/xz1WQROcQzUZ9ZyYWbW*)
(*0yp9TAS0enKGitCRxedv3/4gAqohac570Zfkefw7B0TAahVE92tyoYXaNqn+*)
(*VgQGM1Wan2lgfRsKm+qsFoGdh4LkkllcEDxdN7MmWQQ2veujjqtgf97r+qtL*)
(*UwT3n/WTNhkuCOg2UP+yROB9AF0jEr3zvKacnLoIKJrIuaihV1oLbPFTEoFw*)
(*54re6H+4cBjuZf5m4vsbbt9gSXPBJJcwllwqAkIztmv+Y+L559twwq2DCr4D*)
(*gsZdolz4MaMp2n+ECkLlkRoggP1XslL/0zAqrCIvfvwfWkVNsqDgIBXMHDaZ*)
(*DfFjfL2XGoceoEJNesfsEjS0d+yjelFBxzGIYbqMC/Hl4XX6llTYOFfoc5eC*)
(*/W5offDFVVRYN1r+yew3B8x+Ob+2LReGDzpTw/bjHND45Nq6v0QYvJcFLunE*)
(*+51kp1vzyWJhkFzi/tIFPV3u+bz6ljDsq5hx8cD74JXAPbUqWcIgftTK0xvv*)
(*l3Mfw0ooJ4Rh+HfSvBvep/M7LqSUOAlDjX7KYWo3B4TL39iS80IQaqXXuP0Z*)
(*B3o8cg4arhWCT0ZW9o2ZHJB7XOP6IFQQ1p5YkUP358Dus/HW/Y9XQO8JjR9C*)
(*2hwoscsdCZkTgDMtojFHZ7Xg4HOaoaOBABQ84VdiFWqBkUaz7rUQfuilTJ0o*)
(*stUCfkEt05yi5eApH7TFY3ItWFycpI19WwYiW9Iqz8SshWTxn8ey1ZeBdryV*)
(*tgTfWrCxdCrJDVoKLvKTkkJHNEF4ymPa7tESkIrX1V42qgGvY++3B87ygcN6*)
(*h3USGhog0ur1X48uH0gvrXoXEMIGMiGj5co6PnBXU1o9F8AGMbM2eS8tPhD7*)
(*IdIQsZ8NkpUmLZ/U+aAh/o5czC42KN+Ulf8hwwe5Mt6iEU5sMIl8/5JYzgf+*)
(*rj4b/tFn4/3eVda8jQKtxE9eGQobDo+nHFrRSgEh8VxyaJ4F4defv3jZRIFP*)
(*PaFFUd9YEC2pe8i+ngJ+wRMn7o+yIJlf4oXnAwqIShy79KORBfn9b0PCr1BA*)
(*zcyRp/U0C25nUJ8bXqZAqprlOr1YFtxxNPuH7yIFbIis4asRLLjfVP7sVAoF*)
(*nkH25d3+LKirOCd9IZYCD9PPONdbs6AhuCnYJYoCMh68TFELFjxfs+SZ1DEK*)
(*MN/qme8wYsGrvODgvDAKRJDqXeNsFrR7FDTuPkiBpojv69X+ZUEnY0hKNYgC*)
(*tJGCZXtkWdD9hhE86U+BYtk3mlclWNCXbNt4148CdQd/QgfBgv8DlDAtsw==*)
(**)
(*"]]}, "Charting`Private`Tag$108547#1"], Annotation[{GrayLevel[0], AbsoluteThickness[1.6], Opacity[1.], Dashing[{Small, Small}], Line[CompressedData["*)
(*1:eJwV0Hs41Qccx3F+x/1ycORyzqFFomxKxe+cn6coCVOsHDOXNiujVtRyXJs5*)
(*5JKYdRK5HHPJpaZSKlGMj1amU404NSWXFVI4YSo7tey3P77P93k933/ez9di*)
(*1wHfMEJJSUlAz/+7ZefGE1kVXDDmAq9rTzjCm20kMC3n4oJAIDamPdg9zqop*)
(*5cL2XU2kBe0FF3Fuu4SLF0sOOPNouy4eOj6fz8VLB2UqlHbn40RxeBYXunnW*)
(*J1toywRNOa5CLpzc6/VjJh0xtXlVpsKNC9et6/Q+TDnC0nKzT8UwBwlbvquV*)
(*zjjiRIBqrPwwB9ZpU3YLrx2R2llkxbLhILBiY3mvwhGx9v5VEb+x4RPtPfdG*)
(*iUT1bEt6ZDgbva3Sy5XqJPrObeFXK7Ohdir0poo+iazyuoOXS00RUS9QTWWT*)
(*iJc36+xxNQXXRdM8fymJvQ4z8YJBE/woKcgeW0WiUNbtPppsgqDlLqfa15Fw*)
(*q0uTeC8zAcvaJrfai8Sc6Ov2yXZjdIYqRpcFkxg+O+a5JswYtu+HEh7tI9Hd*)
(*lVdRp2SMknhbpo+IhLOdSYpyqRHcnnuIzp8gMT6y6YZioxE6KgODGk+TOJlX*)
(*KLMdXITKNAd5XCuJ5NICT+WURZjeETWdIyOxv7w1h2m2CPPZ4fpbp0j0918L*)
(*nG0yxK56z2+91XhIOxhyJHGHIUYOJVmuXMyD/cInGevnWXirvkOSz+fB9HnV*)
(*WFUxCw9HZn6u2c6DSGPm3lNnFuIuHNuWs5eHlxZhGV/0G6DQRMl0MJ2H4/mp*)
(*5VFJBnhU5aebX8bDOmVp9k22ARb8+B4VzTyE1Yz8O3lVH0w7e47HnzxIq3+5*)
(*3xykj4byzkTODA/hTnZTh6f0wEextpEuHxcHghWP1+ihOEMWz7Xh4+3fa4yG*)
(*v2eicuHRfNwGPtwVgT2+DbrQMn+1rzWYD92GB76s9zpIt5NXPIvm40lIRYzL*)
(*Wh2Mx3utbD7Gx9K29uBrcdpgP2d0bD/Dx+7jWZ8Nt2nhKsH3G23j45KgcjR6*)
(*XhPMQdfHHX18xEj1XfydNVESZChzecXHhtV/OFVHa2CPw/rCWjUKGtoObhV1*)
(*6rhiKM/NM6PgWSTXf/laDeZHu3teraIgNvlHVG6nhuSBY/aFmyhs8wq4VBml*)
(*CuQvWOT4UdCdDpkVNKtgm7r038ZQCj3pjbKDcwy0Dbz2+TKKwjvSR5Q9w4By*)
(*z6S1mLbV+MiKajkDH8qu/3WDdqwXS9Q3zsCu/ecWbIQUOHqRK1wGGMhPruyf*)
(*pr2z0DJJp4OBW8ZnNH+IoTBdm2NzuoCBn/pTJpLiKeh1hx564sTA00TmMusk*)
(*CqxsSVcZj4GvFFcC/Gkbu/dahTrQd8MN32TQNmtx7XphR/e9+ZUYo738jIXV*)
(*2yUMFGlFBFeKKLgmD90zUGdATdM0wjiF7rMPtvDoJaBpnJozkEohYSI3Vqub*)
(*gFBls0QjjUJijfTOvbsErIrPJ62lnWrmFOvXQWBi9YqOTNpiDfadndcIkNJb*)
(*n69Op1A7/DA6sYzA2k+Ddx/KoHBewpS6lBBQ1unvOEX7or/7R4wiAl1L+hRS*)
(*2o13G24fzSVQM/thiHOEws2mvMUF6QT9Z3lJE+3fhXeFQYcJDAR4lg7Rlq5U*)
(*uW0uIiC2yRSqZVK4XyUUVsUTkKRGFgloy0LOdu6OIVAwMvkigXYf55n5x1EE*)
(*0sdVmeW0+x9whPL9BGbDq1U7aA+KfTvr99G99k3dE7T/AxN5kWE=*)
(*"]]}, "Charting`Private`Tag$108547#1"]}, AspectRatio -> GoldenRatio^(-1), Axes -> {True, True}, AxesLabel -> {None, None}, AxesOrigin -> {-Log[1000], -Log[1000000000000]}, Background -> GrayLevel[1], DisplayFunction -> Identity, Frame -> {{True, True}, {True, True}}, FrameLabel -> {{"\[Mu]/T", None}, {"T [keV]", None}}, FrameStyle -> Directive[GrayLevel[0], 20], FrameTicks -> {{{{Log[100], Superscript[10, 2]}, {-Log[100], Superscript[10, -2]}, {-Log[1000000], Superscript[10, -6]}, {-Log[10000000000], Superscript[10, -10]}}, {}}, {{{-Log[1000], Superscript[10, 3]}, {-Log[100], Superscript[10, 2]}, {-Log[10], Superscript[10, 1]}}, {}}}, GridLines -> {{-Log[10], -Log[20], -Log[30], -Log[40], -Log[50], -Log[60], -Log[70], -Log[80], -Log[90], -Log[100], -Log[200], -Log[300], -Log[400], -Log[500], -Log[600], -Log[700], -Log[800], -Log[900], -Log[1000]}, {Log[10], 0, -Log[10], -Log[100], -Log[1000], -Log[10000], -Log[100000], -Log[1000000], -Log[10000000], -Log[100000000], -Log[1000000000], -Log[10000000000], -Log[100000000000], -Log[1000000000000], -Log[10000000000000], -Log[100000000000000], -Log[1000000000000000], -Log[10000000000000000], -Log[100000000000000000], -Log[1000000000000000000], -Log[10000000000000000000], -Log[100000000000000000000], -Log[1000000000000000000000], -Log[10000000000000000000000], -Log[100000000000000000000000], -Log[1000000000000000000000000], -Log[10000000000000000000000000], -Log[100000000000000000000000000], -Log[1000000000000000000000000000], -Log[10000000000000000000000000000], -Log[100000000000000000000000000000], -Log[1000000000000000000000000000000], -Log[10000000000000000000000000000000], -Log[100000000000000000000000000000000], -Log[1000000000000000000000000000000000]}}, GridLinesStyle -> Directive[Line, RGBColor[0.9, 0.9, 0.9]], ImagePadding -> {{74., 15.}, {56., 12.}}, ImageSize -> {564., Automatic}, LabelStyle -> Directive[GrayLevel[0], 14, FontFamily -> "Times"], Method -> {"DefaultBoundaryStyle" -> Automatic, "DefaultMeshStyle" -> AbsolutePointSize[6], "ScalingFunctions" -> None, "CoordinatesToolOptions" -> {"DisplayFunction" -> ({((ConditionalExpression[E^(-#), Inequality[-Pi, LessEqual, Im[#], Less, Pi]]& )[#]& )[Part[#, 1]], (Exp[#]& )[Part[#, 2]]}& ), "CopiedValueFunction" -> ({((ConditionalExpression[E^(-#), Inequality[-Pi, LessEqual, Im[#], Less, Pi]]& )[#]& )[Part[#, 1]], (Exp[#]& )[Part[#, 2]]}& )}}, PlotRange -> {{-6.907755278982137, -2.302585092994046}, {-27.631021115928547`, 4.605170185988092}}, PlotRangeClipping -> True, PlotRangePadding -> Automatic, Ticks -> {Charting`ScaledTicks[{-Log[#]& , ConditionalExpression[E^(-#), Inequality[-Pi, LessEqual, Im[#], Less, Pi]]& }], Charting`ScaledTicks[{Log, Exp}]}]*)


(* ::InheritFromParent:: *)
(**)


g = 2;
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
f[T_, b_,eta_] = (T/V) (q/T^(2)) (D[lnZ[V, chempot[T,0,1,g]T, T, b, 1, g, eta], b] + D[lnZ[V, chempot[T,0,-1,g]T, T, b, -1, g, eta], b]) /. q -> (4 Pi/137)/m^(2) /. m -> 511;
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
eta1 =1 10^(3);
f[T_, b_,eta_] = (T/V) (q/T^(2)) (D[lnZ[V, chempot[T,0,1,g]T, T, b, 1, g, eta], b] + D[lnZ[V, chempot[T,0,-1,g]T, T, b, -1, g, eta], b]) /. q -> (4 Pi/137)/m^(2) /. m -> 511;
bValues = {10^(-3)};
Plot[{Bc B[T,bValues], Bc f[T, bValues,0],Bc f[T, bValues,eta1],Bc f[T, 0,eta1]}, {T, 10, 10000},
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
