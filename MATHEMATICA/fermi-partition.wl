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


g = 2;
Bc = 4.41 10^(13);
Bc = 1;
\[Eta] = {0};
me = 511;
b0 = 10^(-3);
loT = 10; hiT = 1000;
Tlist[n_]:=Table[x,{x,n,9n,n}];
reversal[T_]:=-T+ loT + hiT;
bValues = {0,10,20,50,100};
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


g = 2;
Bc = 4.41 10^(13);
Bc = 1;
\[Eta] = 0;
etaValues = {10^-9,10^-6,10^-3};
loT = 10; hiT = 2000;
Tlist[n_]:=Table[x,{x,n,9n,n}];
reversal[T_]:=-T+ loT + hiT;
ff[T_, b_,eta_] = (T/V) (q/T^(2)) (D[lnZ[V, chempot[T,0,g,\[Eta]]T, T, b, 1, g, eta], b] + D[lnZ[V, chempot[T,0,g,\[Eta]]T, T, b, -1, g, eta], b])/. q -> (4 Pi/137)/m^(2) /. m -> 511;
bValues = {10^(-3),10^(-11)};
Plot[Evaluate[{Bc ff[T, bValues, 0],Bc ff[T, 0, #] & /@ etaValues}], {T, loT, hiT},
 	PlotRange -> {{10,1000},{10^(-30),10^(1)}},
 	Frame -> True,
 	FrameLabel -> {"T [keV]", "\!\(\*OverscriptBox[\(\[ScriptCapitalM]\), \(_\)]\)"},
 	PlotStyle -> {Blue,Red,Directive[Dashed,Black],Directive[Dashed,Black],Directive[Dashed,Black]},
 	PlotLegends -> Placed[LineLegend[{"\[ScriptCapitalM](\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-3\)]\))","\[ScriptCapitalB](\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-3\)]\))","\[ScriptCapitalM](\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-3\)]\)(\!\(\*SubscriptBox[\(T\), \(0\)]\)/T\!\(\*SuperscriptBox[\()\), \(2\)]\))","\[ScriptCapitalB](\!\(\*SubscriptBox[\(b\), \(0\)]\)=\!\(\*SuperscriptBox[\(10\), \(-3\)]\)(\!\(\*SubscriptBox[\(T\), \(0\)]\)/T\!\(\*SuperscriptBox[\()\), \(2\)]\))"}],{0.25, 0.25}],
 	LabelStyle -> Directive[Black,14,FontFamily -> "Times"],
 	FrameStyle -> Directive[Black,20],
 	Background -> White,
 	ScalingFunctions -> {
       {-Log[#]&,InverseFunction[-Log[#]&]},
       {Log,InverseFunction[Log]}},
     FrameTicks -> {{{#, Superscript[10, Log10@#]} & /@ ({10^0, 10^-5, 10^-10, 10^-15, 10^-20, 10^-25,10^-30}), None},
       {{#, Superscript[10, Log10@#]} & /@ ({10^3, 10^2, 10^1}), None}},
 	GridLines -> {Drop[Flatten[Table[Tlist[n],{n,{10,100,1000}}]],-8],Table[10^(n),{n,1,-33,-1}]},
 	GridLinesStyle -> Directive[Line, Lighter[Gray,.8]]]


(* ::InheritFromParent:: *)
(*Graphics[{{{}, {}, Annotation[{RGBColor[0, 0, 1], AbsoluteThickness[1.6], Opacity[1.], Line[CompressedData["*)
(*1:eJwVkns4lIsaxZnRfMOMYtjGJUIGk9xySdhCzDejKFTbdlIoWyKXJFG7U6Hc*)
(*KQq5kyK7J0Uq1VjlUil1XBKRzqGbKMR2jz3nj/d5n/XPen7vWq+Gb4ibH0VM*)
(*TCxGNP/f7SX6P9z81fCkheF6OEIL3m3p0+V+anjzVHp9caQWvi9NLCzuVcPF*)
(*pAaLl1FakPSqIyq91fCUH5qme0ILtmxSleqphqDHHYn9Z7RwI8mbX71FDRkL*)
(*Eo+3Zmsh+XBGgdw6NYxkfjtM3teCI2/WqXNBFfWnMmRSxTnYyx0O905QhbTe*)
(*j5fncznYHxFomaGqCu9dwf/O0NdG5kqK52D1SvSvsAqdEmojiFovsLFbiW6n*)
(*5QMH9+hgjbvUm41dKigJ4LTSZnSQOXDMdIuvCuRvBK56VaQL7VbqcMSYMqJG*)
(*nMxcnbgoUFfeyDujDMv5qWiDES4eVe/cqyGnjLbJ6uD1RWvQckODU3pVCaOn*)
(*bnVNb9HDPI1qa2+jBJmu7HbuvB6686X0Rp4rQms+01BXuBa5+UXV1j6KcJjt*)
(*8XCN1ofku13M7G9sNG+WXnrubIDCSEOhzHE2ZG+fV+zlGGLlsXRpIxYb0z5B*)
(*0XXThnAQsvU9CxXQGx8IHDKC2wH9r2tMFZAa+cBY/LMRmDsDImn4BW8HWxZ9*)
(*Ao3xPe++15LzL0hy8ag5P2mM2i15/ysYlMd9o2vND4+tA9e61EwQKo9qhcdC*)
(*KtsEpzMtdC5JyOPyV5dB9k0TUPkfwui5cqjY6P2j2NUURla1yad05XCcEe2V*)
(*wzDDpiwFx+KHLNhP81kRTWYQN8lvfS5goSlMPfBdvjk+cy1vyHXLYlV464J5*)
(*zHrUJPgZUIJlsSKu4ehEgAXWcu/8OT8ng/LZ2xaFRzcg0TjrQG6aDCLD+8yv*)
(*JFpC5kJKlYe2DJyfltkkXbdCTOj2FLv9K+A/o3N8ROJXVMdfm1TLWI6c+syb*)
(*1r422KDuvCD7lzT4eZlz0mK2SE5gbnTrZOKtjebIjQ5bGFbMmLSKMbFDGNa0*)
(*u9YOXL7xMy8VBjoYbQZ6jfYYvXsmx8pRCgraoQET3ZvArNCVMvSSBCtnP523*)
(*2hFRW6fXekXRwd2RUGKawIPnz2BlxhE6FmsVi66k8mBZ+Ym4d4iOrbJ7ypQy*)
(*eVigdQ3IB9FxTn35EUohDyeFNdkvdtPB5o5l9dXwEK9/SMLagQ6XeuHni//l*)
(*IaB3aHzIlg6f6XuJjE88OCX4vM/6lQ6VNK+BE8M8MD9uq5swp8N4OLLMe4qH*)
(*tDzD0EouHSa/X1elM0lcZHzrVV5Bh1dNCdLNSUTe2/fsKYOOZ5pkT5MVCQ//*)
(*vtojdDqY+VotP2xJKDc8P9cuTkeQmdpZ480k8qOvkYmTBGJ9mzUM9pA4oatp*)
(*ZjFOwPLk/JD0PhK7u3I0P30j8OB0f8z7/STU1yX8tPtMwEtqtmX7IRKXh/yr*)
(*Z3sInPOqeBISSyIu633R1S4CPzlqe5vjSfzh+Fvqjg4CO417NOgpJHSLeQE3*)
(*XxDIPaYn73SBRKUHZ1VAvci/5c+3ty+TSKblM9kPCPi2nz7gXk4iqEZ+rvEu*)
(*AaZq1PU3lSQMZCReq98iUCnkhmfcInHryUBi92UCra++D3UJSaRf7X1pU0xA*)
(*s0xS9sMjEsFnO1ll+QSoR07pdDaS0CObL4VdJLDN/bcTe1pISOnU9785T+D3*)
(*2fechRckvtDuatqkEVhtp3jmxCvRfc0V16TiRfw/rKgWnSROXyn5HhpLiHLV*)
(*oId1kfA+k7vuzUkCo5Em6RndJFR5KXWlUQQa1j/j5fSJ+mkKbbEOJBBg7mHX*)
(*OEjicFnA8lJ/AmsrCzl/fCThFufrJrmPwAbrgZLZTyRWOG5/+3oXAQ1zT6nZ*)
(*IRLxjeu/BrsQMNpUNBQ6KsrzspHBaycC/qVf9e+MkXCI5R6yIgm451xw/Xuc*)
(*hLiDyhxhS+BknHjevyZJRDcsShYbE9jhEDvZPyP6h9JpZ8KAQJ+younELAnz*)
(*mLFzB9cQqGibM6XOk5iwH1CyXE0gadvSa/ZPUf6Pm7gd8iLehsR7KuJ87H6U*)
(*LFg2Q4NkQXeNkM6Hfb1yQoCQBuz+sH0zm4/9LDGK9H0a3OY8eBqKfKT4fYyu*)
(*ukNDYHPTsmmR7mZUBU1X0VBh2LxUrMxHsIfjtrhSGmJ7xHlTqnzkjIewi+Np*)
(*iJNKDrqkxcfY6qaybjcaGkJlrXqN+SiID3lMfl6GZgennx5b+HCij1/3ll+G*)
(*vHu+7jtP8uFM1G9IM5RA5IxEZPtdPlS7K8QYnlTIWHPasr/xsU7uP7muO6lo*)
(*k3l0KfY7HzyXKfNsNyoUxM/6hY6KeBrtD3I2UxE4mr5IjvNRX9Xbs9GaisL4*)
(*eJvpST58EqSrD6tRsXAh++8dC3yUWYb59Q9Q0DhZbi3PEKAuIkuM856CRMNz*)
(*a8SYArysepgb2EtB2ZcMpRGRntKWap/poEDX4cNcw3IBSLlSG7kmCoSuf7WE*)
(*swQYGu5k869SQDHzLOhQEmBRe/5WaikFOjMV+fXKArB8NVxeF1LAOsYsqFQR*)
(*wLL7YKxvtsh/UaIkRlWApEba+PEECgJZ9bUmGgIULa1NboijoC41Q7hKU4Aa*)
(*S3cdqdMUDEsdfcJYLUBfVaFXVjQFRiEBfYNaAowNN828O0JBk3TEl1ccASR0*)
(*RjK0wimYL0/5+762AIq+LMPAEAoydGuo5ToC6OdbtNwMpKDA/CMrU1eAfwAk*)
(*NfA/*)
(*"]]}, "Charting`Private`Tag$154708#1"], Annotation[{RGBColor[1, 0, 0], AbsoluteThickness[1.6], Opacity[1.], Line[CompressedData["*)
(*1:eJwV1Xk0VWsfB3D25uxjypyhSNKRkpSUOJmnTO0jQxFJCVEkwy3qllBKUbkR*)
(*jkIS3Uo0KDe+RINwu0KkW5cGGVLSoBLv8/6x116f9X3Wen7PPDswwj2IEhER*)
(*aSff//9thQs/uQdrIl7k9T2j83wE/JPx7UKQJlJEazLqSvgYmRqbmNykiQT/*)
(*kLv1F/hIMPCmPAM1cTSj+FdjKR8SfreZiwGa+FsvXTqtjA9LFQcN2kcTWmUH*)
(*FG2v8HHlSIBjpYsmKl7EeAiq+EiLPpmvuEQT1qM+aZYP+bCz/+7UPqGBsWfH*)
(*2uSH+dikN7QzIFUD7XkPOLKmKxESE2Z6UkMDBu0hXrzElcicSfm8qpwJyQ/T*)
(*e+hHKxFO164yt5oJxYRfL7VmmGP+GsmnFp0zYP2qauJusDky++KXugTOAOf7*)
(*7bCJanPwWuihmI/q2JJ8fPkNOQvka6lb2KeoI/zXdN/hjRaoq/TaNFtRHSHX*)
(*vLj91RZoujJ7blGJGpb0U1c1JSzxk0NbWpurYfVPqad1Sy3RJZRcMPxIFT2x*)
(*ucZ+QZbIFZ6t5G9URV1eBxOXaQmJf9dLZ79XQUBr2Cqfh5Y4E7eoRi5BBVfe*)
(*TBYIpiwxMz5DxlBBBR77PdvqjKxgW6Oy0OfMdITLbJU2i7CC+9aFg/OXTkfj*)
(*w2O3fEqtIO0VGseBMnSCo63ThqwwklftN+WqjF15gzt5hta44ZLXm/9KCeKu*)
(*avTVndbQ4xcZr4pUQkzJrn1eddZIzDTRzRFTQnega7LvNBvQjq93cHMV4Rob*)
(*513obwNDsxtp++cpwurTsWeF1TawyZpuV3BHAcLezro0VVuIGglbHq1SgEvv*)
(*vtSmJFv065leUeySh/eKj0lp32xxLTXIgNoujytZKnurttpBX+/mnp8/5BD/*)
(*YfmUySc7HF6ctTU3XQ5NC25EfNltD7k/jpav5cnh07KbG1qlHXAg0uOoVYgs*)
(*oq7LiIZVOaDyUNlnzZPTMPxbi6ZesCNWaLlOyP8pAw/f9TW8ZauQlipt4d4u*)
(*jbak31s6uU5YVDpu1CIijc0qfalmw07Qc1z80G+GFI5vMR76Me6MD1Upp83s*)
(*JCESK6f1eporpEvnSS7yk8CHvlaL6HVu2LX6m77fLi5q4+1MULsaPr+2q0vF*)
(*ciGX5IofDathevEtcyuKC15MfI5n02pMcDr7lMK5WLHu+eGo9tXYV3Mtu9mf*)
(*C6OjSpF7B1fj0MIoMb4tFxYHo5+lSbII7RkYHbDk4m3h/VFzGRZOqRtfZq3k*)
(*oudIxdJHciwWLO9utjPlwma9qu8iRRbSb9jbY8u4sI/QR4syi/S8RZEX9bj4*)
(*T6bZs1uLxSmp9z3qslzUl6iicBmLuFubHz6Q4qL41pi3pwmLtcHPb8RyuVDz*)
(*2SA8asZihbLHOR1xLh6P3PtNfSUL9buPjreJcuFslLSZb8VCuLvM4fBnBhk+*)
(*X1bdcmKxd562sckoAz/vS87PXVj4d57WfvueATdEfKk0y8IiSV4uc5DBryiJ*)
(*Yhl3FlpLUn9Z9TOozDa1iFzD4txAcOX3bgb7Lztk3lzPIjnr5dmSTgabY0+H*)
(*Fm1gscXO+5jnEwZqj9YVtQawcBhrjacfM+BNXfbM3cRiXoF96NVmBo8VBT8y*)
(*g1hIrK7x2vCQQbXsqQ9bg1kMThjbytxjMF/EqmkklEVT2aXF1fUMtDI1PiOM*)
(*xcW1c2eF1jKw+X7PNTmcRRpHKK3yF4O2sL/GuBEswq8p/WioYnD5wjPdsztY*)
(*uAam9UddZ3C2O8TfNoqFgZxYh1YFg/4H1LzLO1lU3O873HWOgeyJ1LnGu1lk*)
(*lPS0mhcwqMn3T0iNZ7H9YLtCsZCB4ubX+Y0JLFyCW7ykchjU77xQZ7yXrK/D*)
(*vZwdpxhE3lW/ye5jIalb++LpCQb0tNKYlP0s3nGqtM3TGYwGjZtFJLK497Z8*)
(*y7kjDA5O/3JdkETm815pmeQhBuXVH7WyiRPPF45EJjE4HbGb9k1hEZCSu+Tp*)
(*PgbGakPhRgdZmG/JjF25h0Gr/vp6tVQWGvZHbxftYtCu+Uu+g3hibsqkRCyD*)
(*iN7OlHTiZ+K/W0dGMViazR3xOcyi6k1cSud2Bn95zypRTSP7rTGyiR/GwMNQ*)
(*7Oc14uji0GlFwQxWLtAIjzzKwj050F1iM4PM707cz8SGQb6nIgIYqPhE/JeZ*)
(*zkLWzuNZx3oGk7naYiszWIzouGry1zFQry32vEXcLGYfWOjJwNHbduOP4yzK*)
(*Xpuf57qT+pU/tOcTH2pYPrjdjUFRJm/O1xNk/5wzNOhwYiB0avRLPMnCNkkv*)
(*ysyBQf6S+OJhkmtv1r5RYMPgUEWrgX8mC1HbGT8YSwbPXaQc7xD/N0fJfDuf*)
(*rG/IR57DHyxqaJnEdhMGLabdM8+cYpH3SvyeqTEDg5LkLIMsFrvvTkoULGaQ*)
(*o7wGicRri765MgYMJNfJBtwhXnbg4/Ft8xmk+P58r5DNQmnTQMcTHoM9t44k*)
(*eROPWfepmc5hkOvhNll5mkWbdo/f2VkMkp6oZZrlsCin2gs4MxnI2f+WtYM4*)
(*vL5R74kSg/dN5y21csl9UVizbYU8cbHr+6155Hwk3rx6RoZBln983QjJOYHl*)
(*X8QlGfIeQU1cyOKNVemKcA6D31fqRu7OZ9Ewu3BPG0XqEadfXCR5oWhunckU*)
(*eVGmTicOknxf70nxMz85eLymt8mZ2L8ubZX4OAf3r+s83kvML0g+GvaZAz3M*)
(*9xA9S+6P/Xv/+ecjB3t/012iTTweEKds8p4D+VkLAhIKWDy1jFyXP8ABb2w4*)
(*v5fkLNaELu/hwN3tjPxIIQvrWvXU0BoOthjKDI6fZxGiIELJVHNQrdEcFlfC*)
(*4mjQm93lNzl4J2d6+wpxRVXT2JprHMRLTSv9Qdp3SZWHfyvnYOKwY7rEBRaT*)
(*/n+8ybnEgaThtitHiOdU7PY3LyPjW2YSWkDsKB7wtPc8B3Gtlw60lZLzutaO*)
(*TS7iQPhd0s+1jMXp0QiVgkOkvwNv65wvsfg4p7G4y50DZ5t0pq6SRf6hiHqH*)
(*fnGcy8l2m0FG48QdvRSgJI7+V8Vqas/J/cHUrkhfJIYbvrOD8r+Q89ZVKiLl*)
(*Q8N2VfX7JHkBlig+zhV40fD2kxL/Smzv9nVZtjuNzsZPs7kKAvikzmx74Uaj*)
(*+M7fZpnKAmxvsN4215nGQCb3v2KSJ06FcMMdaOQMK59SVxTglGl6UYUNDcvg*)
(*VC9l4rKY6+bfLWhYJPqq+hDXlvd0W/Bp3H/DMaKVBGgfEo1JMaGREeAiq0P8*)
(*jjdPrmUpjfHHJ5qDVASY2Oh2UXExjRmaO9xqSS4njLb3WUjDz1Va0oDUo9OV*)
(*03tWj8bjtf2e+aoCmCjWJfTPpTHVdVXuJ7GLW7+KgTYNoW5H3t/EG1NlKqM1*)
(*aZg1DEyPUhMgtsHIrVqdxtMTVcVHiA9PrRsQVSHz88Fj2nfifNN9SY6KNBoS*)
(*Q/qGpgtQEXN+VrosjX/042Sq1QW4V958u0OKxudEi3f3iZ8NffKcySX1vtvk*)
(*b6spwAhPbTRQjEa6lkG4/QwB6ECLtFIRGv9GN++pIbmqMEj34wSF6ynHv6nP*)
(*EkC/60j9su8UOjfd6NWfKYClYoXfni8UDryowkliD7eu8bujFLb5KgV6agkQ*)
(*kjp5UnKEwktn4Utd4oQGnUWCQQrfqp78eZ84Y8qpKestBQftMz8KNAQoNt0R*)
(*9KKPQuD0bj0Z0v/tmCyRuS8p1LZ2q24lbi2/kxvWQ6H4oaxO/WwB+oZeLat4*)
(*SqFjsuMNT1uArzzJtvEnFC4ZXlGXIJYMNNxm8ZhCn0L/u+3EmkIvbkozhXe5*)
(*J+Y0Ei/pSihqfkCh4drwZ/05AjgoFpkrNlIw05S/WcETYL3bw+51dRQ8WsS6*)
(*/iV5ZOqH6LN3KNiXu6Sb6QqQ1KAs13+LgvYsHTW+jgDZU2YXF96gsFE3489k*)
(*4j9NA+2jKyj46CY92zBPAMQc6r19mcLBluD7R4jbyy8niF6kUC5lsreYeGCo*)
(*XcWxhELiz+6C18STvJ8Vx4rIfN8dd0nUE0AhcLZbxxkKeFEwX0jMEzoMzMij*)
(*0D5mnJFL6jXt2pYUmE3BWfKruNF8AY40cEYTUik0Xq4Z3bJAgLNT+ml3kym0*)
(*hKxrukh8zXSNrmQihVUv/oiOWSTAg5hd9exeCjovBt2d9QV4Xn7GL2s3Bb0N*)
(*vPoWkn8cahz/N5ZCvljfk/ckF9MdPqmzkwLvTdMn/4VkvwQqLAqLoOB9iycT*)
(*YyjAQqFJ09UwCkbzJk/FkPx/+CtSTg==*)
(*"]]}, "Charting`Private`Tag$154708#1"], Annotation[{GrayLevel[0], AbsoluteThickness[1.6], Opacity[1.], Dashing[{Small, Small}], Line[CompressedData["*)
(*1:eJwVV3k8Vfkb5sq91shaOCdE2Ze0yHaQPfdkKWsilIw1SyjayFLKtEwqW9xk*)
(*nUlpWph4Ck0qGiWhZFLSqsSMFqPf9/fX/Tyf532fdznv+55zNULjvDZzhISE*)
(*NISFhP7/+6DC8LNXBI20iJGvbmIMQnp+nq7eTCOrVnFzqwSD8R+TM7NhNDSr*)
(*BWeUZBmkG/ly1ofSkKoI2TQgzUA8qIlXF0Ij6mXMrnzCnzhISQsH00g+e57T*)
(*KMlAq2mPnG8QDctj0fPlxRlcfP1C+ddAGrtuzLQ+JfFslZ0pkQAa4wvK7AVS*)
(*DLodazX9/Wj4q467FBD/DUnSOud9aKzem3RIk8R7WxFvKLqehiDgu2L7PAap*)
(*PQ+XBnrTQIuyzDXCc4VWml/wpLFsi0wHn+gdNzptzfOgkfu7rNN9eQaaQf/Z*)
(*B7E0IvecmPWUY3D+YIhLozuNFSqKllcJtmpq44uvoTEwGvMHn+A7rxd7B7vS*)
(*COI89PWSYeCnfMDvd2cahy2Cnj8j9b5y/BAk6UTDjKvLhBP7pCSPsE0ONK5K*)
(*aaUHEHuOoHHrFXsadsvipO3mMvi5RylW2o5GYOr86bvKDGihHYlhDA2ZtVa+*)
(*/ooM6o2GUq9Z0/ion3fZgehbBNnukrGi4aMqveqoCoPbBwWZmy1ovDm3c0ZD*)
(*icH6Jm5eszmNsFH3tlSCR15HHp63kkZTt8mafAUG8cpdxyKW07j1RKzuNIk3*)
(*62hy6roZDd3+e/RW0o/8pGOl8ktpvHPPMdci9iqCfwWRJjS6UwJfnCb9re7x*)
(*r2k1ouGSLLA8ReKvELr+m6IhjXkpZdWVJF6bkfqlKH0aodpXVXkUA4+gzGs3*)
(*dGlcmlULuTyfwbODr1qUdWhsStXg7FdjEN3k2h6zmMaE76LMRsJ/e13f2aZF*)
(*Q8+lJ0id+Ocqy95fsIjGea8punoBAyWnxN44DRo9A+8PpKozOJvUN9CxkAa/*)
(*x/bRFU0GpoJVw6o0jUwX/a4pot/SU/xymxoNt9/m2j8j2F1I6O2fKjREv3Mk*)
(*soj/oFHYR2oBwU+kO0RJvK1Bt6YSlWn8fcarNZLo/XtQ91unIo0U6QA9TcJn*)
(*NeX/WKhAI+ee+BtnmoHcm49ztsvRGJLbXVylxaBM2VviniwNh5mNEn0LGRg4*)
(*XZbRlKHhrhNY3Ejsm5IWKKZK09jPSNXEEd5FkK7SLUlDdqw0/YA2g76e4YVa*)
(*EjS05+yL8yZ8uNBq7R1i5PmOyDZUknw+G53T+4tLQyJzwr2TxNsTJG6yWJSG*)
(*mt2iliyC5+ZHL08XoVHn+cH6lS6D4qb7Fg+EaegYhRhyib7em6W2OkLEP712*)
(*uG0Rg6vKJxx3zVIISxptP094R6evbr0zFBbySl2MCP9AumQn9zsF73VDVyZJ*)
(*vzY+sq03/0ph+Y2U+b2Lyf4Vv3z60zQF9crCkVR9BtvDc6VL/qGgv8Pa6I4e*)
(*mXcDA5v7kxQmPVIS6kh+BZP3Y4U/E/vurU9+EH+15sQys08UrlJPz2QakPna*)
(*p/zX5nEKB1UHF48S/2VuzUIn31NQkveWkTJkcGNesOmdtxSUTcMb35J47AAn*)
(*dOY1hZ6HCRm3SL2DZ84dNRqjENwxcUvMlMGWrW5tIaMUZOJ6fhwl/KTx+OTR*)
(*FxROb70z+5no75k+otXxnILEsff/fddhINW6fP30MIWKNHueG7E/mT2wX/cZ*)
(*BavxA2eciJ42m3E58CmFQZVrHjUkn4uKGmOHBikwrz0rnxLeZqhdGf0Udiz1*)
(*t/AzI/fi7FaXz30URDcIBsZMGPhGS6VpPaKQlqPvU7OUwUuzhhqfhxQmliW2*)
(*bCD28d+9B3N7KGgWHWp1N2cwc3Naovk+habqJ9u1ljHIO1Bk+aGLwqlpr0uy*)
(*RE/Ji4leeI9CgPCcZ+uNGFQseFHseYdCwvsLDwuJntHz7K7M20T/UXkYZxWD*)
(*5mq92d9vUfDide42Ifou8d1Gr9sp2Jb+nrWO6PeuTAhWaaOQtK3NxnsFgzDd*)
(*d4kheRT8bn9Q93ch91NPlLmQQ8H8J+MNIQS36y+U4GQTvRWDtWvdGEwZrHrk*)
(*lUXhYtK3Mlsncu+NvM8I9lEYzgz+OET4dcYxUVN7KNwYUTves4bsk0nOCsfd*)
(*xN/E5PFzZwaXTMuFT2RQoHRK21SJ/sulzfde7aQg4BoHaxN7hWWPClfuoFDl*)
(*pLjalmCH5R9Dc1Mp7OoVWjfhSu7tCnGjge0UMq3bTzcR/uzKRV91kymYmfs/*)
(*mCH6vebW7TsSKVxxdn0tyWcwx8K34O42Co2MVEWTB5k3y20BavEU5njcdvdn*)
(*yf5ZHdSOiaWwxf9d9KA7qd+68tP1aAoNLyZ6rIh/u01r89woCsXyGa7xxH+K*)
(*GcgOjqRQvzu7yI5gLbtJz4YIChcW9WUeInidvTQlvIWC+2VBfOlaUv/qJa89*)
(*w8k81nh3/OVN6newa6wIpZD71OZDhheDUcfAXZMhZH/efBu7QfyVnJNdHYIp*)
(*vCi1iw0g9k4uBQq/BJH5MNqs+png7a41w6OBFAzPuxRH+jA459ZWuyKAQv6y*)
(*zs7rnuS+rBlKzvGj4L/sxkFrX/J+5E/b9vtQ6J81fWGxjtxzdp6U7npSX43t*)
(*9A6Ct6zVf5zmTeHtc15Msj95n3s4VtzxJPOmpq6hRvAtz+AYVQ9yD7giwjnr*)
(*yb30SjOPZilwuVHtHwhevO6YyHV3CqN630tWkXx81v/aLb2GzK/l2fA1fgyy*)
(*ff48tdGVgsv2uuT5AQwu+z4PP+9MYWD79QZ7gl/5fTcWcqJwaPrngV6ClQIU*)
(*v3s4kH5o2epZkfhOgca3yu1J/zOWO+UHkfo3uB75bEvhjzHmxkfCnwsK27Ca*)
(*oRBzTK5g1SZS/8aMJcetKXwKM7Hct4HUH1L4+aUlheRR7UT9YFL/pgvXl1tQ*)
(*8FTYsuRfgiNC7+Zmm5P7E2j5cojonQwb9X68ggJ/1UuBJvG/Hf6D1llO+v/f*)
(*D+dZYv9l84K3qWYUStseVA+TeDoRZr93mlI4MqC/02oj+V7Yyt+jYkJhWlKi*)
(*8xjBuZERa6KMiJ5r7Yp3IeQe/7RX6Q8DCnV7lrHJRO91VNFzKX1yb6rK6h+G*)
(*M5gf83t9kC6FjvTpHxeJvUvs/ZTfllAoD/k6dGoz+f6Je2P/Q5uCz5awvRpb*)
(*yD2NF5nroUXy8Thy6BGx799GDZzRJPsXG566P4KBWOLKsxPq5P7G+kfGEn3z*)
(*JM84+4UUxvqe0AXEf2tylMUxiuxHr6i8OdE/rsYJeNGohi0O37Ju7yfvb5FW*)
(*Vxs7NbTGplnnVZH3i7fEY6ZPFSX7ij3F7hH7kZ3L3ENVMbyqfp7hRzIPXSLv*)
(*kj+poKri8Lk3ErYoVVdhnLJV0PXBIl2PssWNRp8wDXkVpFhkNwyb2uLOeQ1t*)
(*QdUCxB8webDLyRbfuSK29jYL4BoYM7w90Bb9JRL67+/OR8WFJTaHE21RVHKm*)
(*0WrTfHjuCtYLOmgL8aENUic/KKMosVai9ZwtylKMW2TTlfFg9s/E4222UNv5*)
(*s7SJnDKOvZMP8PzbFg4tyoYBZUpI7PafeiZsB6+fDN/qLVOCtYfIP6cW2UHK*)
(*JzKFC0WwEfPStdfYYby4OegHXxEj8V6f4pLtcNm9+HnpCwVkJacZe56xg66V*)
(*YLlrvAJGYjP2zD6ww77j5ktOz1FA9diM91WuPURcXm4TK5IH9nNjpW3sYWJ5*)
(*OX+vjjx8px40Dey2x+pCJcfy63KQXWKr79ZmD2Gzkq67rnKI2qxTcFN+NcZ0*)
(*Lc7L98/DEd+c9fOjVuNS3mYjTuw8fKm+9GGofTUMdK9kfP8mC72ze5wCTRxw*)
(*wLTwp6ICWQwVVaXzyxwg+8uhBr/FsljpdyZAXdkRmfHrDtltlcHRS5eEJ845*)
(*ojG3doo+Nhcns8b1cy2dsEqdPzOvXhp/J6y2XzzuhPw8KcarVwpH3yY4ulU5*)
(*w7jmi1mXkBRcpmzqp+JcoOti2hmkKokrHc+cMja54uPV7FOWjhJ4mh4yPhjo*)
(*BqkaHQnjIHFkJ9SbZ+atQdraaYOgNDHwihWaHR66I+C/WBXJ7WLQ/q6bwe93*)
(*h0XdK961BDFYuDQl+g65Y4bbN6IQLQbJU0q/xIy5Y0/LpZP3NophWOPtMsF3*)
(*d+QaJsyxchDDcU5elrwWH5FP3ky8sRXDM+bmZ1qXD7e8TcOF1mJYGR6rqWvE*)
(*h9SoR9PkCjEkpI2Vm5vzUVBsHF+nK4bKy6a9Fu58nJD88ERFRgzl2TfvpCXy*)
(*kXItvPO2pBj+sTxc6pjKh1/E08vbxUh+7mme4hl8qLTdPfJAWAxJo0U5ydl8*)
(*lOyodT4wxcOvl0MNuk7ysUtHc7n5BA/2kV82ri7hY2PfKc1XH3hQsUqcU1/O*)
(*h/rSvP/sxnh4UmXcx9bycfZNROPXAR7irSv2tjfxsb9w+ExVHw+rPv42daGF*)
(*jy2OvofXP+QhbepixoGbfOiUO0VeuMeDQLyq72snH+JrW3yCO3lI6NCKKuji*)
(*4+3McgfpW4R/v8hHuoePOj/thZGtPNxopvl4zEc+t0RK+Q8evHf13Joa5CP6*)
(*ksK39qs8fG5svy7zjA8j2TmP1C/yULh/qPvLCz4u/jlyoP8sD36bPJcUjvPx*)
(*c9WTbptyHlL+9hJdOsFHbE6vXGUJDxf6gyUbJvnQd751etsJHpqvdIn5fuFD*)
(*Yknrs8dHeWD1k8r2fuPjNfeqpk0Bj/zf2HP/yAyp/1ZNrUQuD+0fy15sEWKx*)
(*71zFeHwWD5e+NnQYcliEZBctfbyHhzH5Kt4zERaU06EmQRrJN7wl5geXxYx2*)
(*9qz4dh5E/zBwSRBjMSi62z4+gYe7UWvVH4qzONERf8cqiofb1bktm6RYJFVG*)
(*zhVE8KC996PuEWkWXvtDvcTDeTCTEE28OJeFjOO6wUcbeKBX7SzskmUxrsWn*)
(*rfxJPinanM55LO7NcQqtWM9DuPI/66/JschtX/k2luUh6mlKTroCiy1nTYwe*)
(*ufGwiRd21UORhUOWboKlMw+/1VlUqymxEHZQ/caz5aG0IMeiVJnFjrZZ8XJT*)
(*HtzCFIf9VFj4Cab5PCMSP2Lg2jeCV2R+OhKjx8PBv7urTqqymLQfWWCxiIe/*)
(*JkuabquxiL7ZoftQgYfoLNEakYUsNt7IdxX9woX0U/cTCotY2Leq5EW2cPEu*)
(*4jadocdiq5wQR7qZixLZ9vDXBB/aPLqj4QoXa28E7vTWZ9Ev2RA93cDFAcVd*)
(*MzoGLGL9HD32C7j48Vwz7oEhi1MTccrluVw4jWh3q5uy+LSoo7Lfi4v7n/oX*)
(*Xl7JojQ37qbzmChaqmsqnziwcBOb+DVEQRQbh5VC6zay4PNaVxUYz0H6H35B*)
(*zTvI8+yvEZIMEEFJb9/g9HEWS+X/KvL0EYHG0b9FYn5h4cT+u+Kklwi+dA2l*)
(*jBAc224fo71GhOxRw9i9EyxaG54MMFYicB65OFN+ksWmPOnGJFoEE/IG4e7F*)
(*LCottm1+NsJBUV1GfbmARVNyoZD2MAf7an3/UzrLorvhelHUEw60q3+OyCf4*)
(*38USD7485KDURopOqWThLC+wke/gIKfw8D73KhZv3vUqu1RxoE/F3ZquZTG7*)
(*+PvFwwIOIt8kmcXUsZAL1WAflXEQu9v+3AjBFv0xWaEnOfBUc/ipq57FwXbu*)
(*RHoeB9a/5CsLfmNx5odBftt+DqYfHs+Zf57FJQvvJRL7OLjeeU3hMMFPG8qC*)
(*CndwwDU3G01tIP1+1/FlaDsHj7cqh40TPGfJ+2NaiRxU+r67H3aBxfxQOeOo*)
(*OA5s1DjmAwQblpjfuRDFgeqdf6rYiyz+B0L/aSM=*)
(*"]]}, "Charting`Private`Tag$154708#2"], Annotation[{GrayLevel[0], AbsoluteThickness[1.6], Opacity[1.], Dashing[{Small, Small}], Line[CompressedData["*)
(*1:eJwV0nk41NsfB3DNMN+x3YSMdbImypbSGJN9mzO2VK6EZKm02CouaVG4FKls*)
(*iVHTKEmLpauFeHdt0XJvCymVX/qVKKJNRfWb3x/nOc/rOZ/nnM953h+98Fj/*)
(*KJqUlFS2ZP1/v3/S7IP/ejYeG2QFR+hyEHbv0OSZKDZ+C7xQ/VjisV8fp39G*)
(*sFERwOvw0+NANuQaVR3Gxvvg4U08fQ4cWR469CA2ahbxniobcnDxQJhnvRcb*)
(*K/8zo7zNmIOcbfnlKgvZSGvsfGlqwYGb+zfycFoHr8983jDLnoMIk7dbw7J1*)
(*QC0d/OdBKAcbtm/i5uvooDCkfYC5i4MCbVrQy3ptbGe6LfcVcrCZ3sK3d9LG*)
(*48KmL5rgwHS53COHXi0MFHxtzXghqR/cscgrXAsGJbFmwwxbzL1Df7t9XBOn*)
(*q5t9fiywRbmupoN7piYeB9RvfbPMFjfqAyL0VDSxW+SWcHGHLbov6hmJKzUg*)
(*rqi+mVNhiykG3dHZXgNT1vusiu/aok8oN//dLXW4nquLuvzTFqXCE/W8teoY*)
(*eFFK8c25kH0WrHB0lIWBmYRpEMXF8SSLZqVUFszZif4t5Vxo7zikaKnMQrSI*)
(*3uXSx4VrM8ss6LgaBsN3rvLRsoP/RrMR00VqGPqgP7871A4KAdFJDMxG+9gN*)
(*v/fn7TBW1hjyy3s2uD4H7qTSeGjwKntR/lIVan77xuwCeTDhiRfz41Tx9cjB*)
(*uJLrPOwt4Bgfk1bFoNwFs5lGS0H3/G88s1QFRfX0m1fyl8LSriEnbZ4KGhc1*)
(*tfux7OFSrOYmuq4Mx/NH7liX22OGtfDOLb4yxM1F0voODhgy4V5U6ZuFzqSc*)
(*xtujDriUHWVOi5kFLQUv2nk/Rywwubxz6rsSOi9wuk2bHLHfqnhjaZ4Srjjt*)
(*cNae7wSlwtyawLlKmDokHBk76YR9cStynTbMRMLqwAstC51Rn3X2Ezv/N8wh*)
(*7KqATmfY6npPzzqniPaTvrM+J7sgJ1vBwf+hAirpaQPq1q6wqPpqfUdKAV5v*)
(*fDe3TrnCxNOqK0RLHuJML3WvUTe8v5JZYucmhwrBXzurhtyhUDVPziJEFnVr*)
(*l71bNtsTyb6TC0KSmVi14tzvKmv4CPoRoymfyIT6B20D7Ug+uNWvqasJTFC1*)
(*Q5+NovmYZvQOqm5mwjDCMoy3lY89zZeO3g5losH+j5r4P/nIMkuQ5rkyUd3V*)
(*4Dl4kY/o/uGJYUcmdBsH60Yv8UGy1w4UL2VCcZLcm7zKh8Irv2sfbZhYYd6u*)
(*S7XxkVdmEVdtIjlvFp2T7eOjSH60X3MmEz2XN0R4/uIj6Wpk1015Jni+9TxV*)
(*aYLA9U8bEplMPDvaNvqISaDZeuvw/RlM/KqSM3VVJhCmnPXY/4lCsW9Mc6kh*)
(*wa55+os5ExTuK/50UTYhCO0t0X89SqEo//OxFDMC3YXZP5yGKEQc/hakY0NQ*)
(*Mby+/ttjCpHJha757gQZxQMnKnspKKf7H9tDCNa5/X5w5QMKz+uGlVf5EMwT*)
(*uUfX3pbUa9xIurmSQNa3OWBNFwXLkkBG6CqCkenFroodFMaWyUg9DSaoDjSa*)
(*E91CISMgnJ0XQZDDECqwmiiEmRVGta4j2HxJ9XvbFQrqhg5rnkcTmCtJ9+jW*)
(*UdC0tjrcGUtQ1zm4v6+CQtPUu/6oFIJDlf137UUUwhWKe6V2EsT8+VD5lJDC*)
(*rqD52LWbYL5Hx7H4IgpZf400GqUTyBm3PH90hMITddvFfpkEbxhX9O3zKKxd*)
(*h3VhWZL/d1SdlcuiYHrcPGNxDsHe0yfH4tIp7FBXCZrOJQjLLF34aA8FQ+Oq*)
(*2nN5BDruudfEyRT2Xv0uc/cIwbRR5k/ZRAovQuXVnAoInsjsdo5LkPS3O+GW*)
(*qJCgqD2um7eJgnP9ORvLowTbTkX/Jl5PYc84f2hNCYF/Rri/bCQF5OVYph0j*)
(*mOm24klPMIVnrZfMy8oIxgy92bxVkvuncoXFQoLb0u7hJ1dK8sudYGWVE2S1*)
(*LRmJ8aFQ0/5tseCEJK8KS/MeQmG5n43LHBGBa7pJgp0HhQMyYRZvJZ7hqvWd*)
(*cqRg8zlCdouYIKX1p6zIisIUSzHV+LRk3sST3pQ5hfoljLQuiW32jR/eYkph*)
(*44Np9/WVBB+dBzW4BhTOx7TEFZ+R5Pd3u8kDVQrJbTVh989K5u1GDl/mKwMH*)
(*RElhLy8QOLdoZkc3M5B6MviqoIFgg7IUTbGRAZPakm8NEudGvUqpucxApOgP*)
(*Sv8yQZ98zebJGgbS8zjCSYljAt38MsQMHPa6qSa6SlAyEcsSZTGQ+G5W/vsm*)
(*gnGD9lN9/gzMs+nvTm8lKM+K/dtjSAbPFxWEnP2XgDAnzoepyoCl0fMleYTA*)
(*m2qxzbOQRu+97l/ZDAF0+qqk5IPoONZ9cF+HrgALVf4tXRZAh4tpapudngDu*)
(*Pl9sjvrTkekXwqiVOKbNeYuRgI66HydySvUFaKnpf+zAo+NHb4k4zlCAtdmK*)
(*9dvYdHiPv/quOU+AU9z4qOeDNJARXa84CwGubS+WMhqgoXu115lXEt+tuV66*)
(*qZ+GUh8is9pSgC9z5e5/fUDDmoutN9ysBPBQEdurtNOQ7nPAVctagOG3D1me*)
(*lTQkJTjHd9gI8HPuVN1BMQ0WH3c+sVsigHK4nk/PcRpSPTNcayXm9m1JDz9K*)
(*Q/TIF80yjgAH2hgTqdk0dLxJvBfPFeDErwU5rRk0fPpnptOQxJe4y43l9tLw*)
(*JnhnbbCdAE9rjocUp9Dw0u5sgQdPgPG37V+fJdIQOhhPXZdY2vhdvuFWGnI1*)
(*6DsWLhVAPVzZYlMsDWnflr+vlNhMyOmu3SR5z3d7pI69AP8D89IS/Q==*)
(*"]]}, "Charting`Private`Tag$154708#2"], Annotation[{GrayLevel[0], AbsoluteThickness[1.6], Opacity[1.], Dashing[{Small, Small}], Line[CompressedData["*)
(*1:eJwV0Xk41GsbB3DNML8xiMYUpkwkyxAaWwzZ5/dIpKM4jhpZ2uxK6YhTvcma*)
(*4rRRUsnSok5zjLTzLUtpOR1rSitvSUWpTqkc3nn/eK77+lzPcz/XfX1vw8jE*)
(*wFUMJSWlXMX5f20/bvkxcI0APYOaHUqOIoS3FX49uUqAx/LXlSkKD098GhuP*)
(*EqBwYGzrO4VVpZep6nABfFWrfHqcRHDXIfrMUAGWW7ZZyZxFOLcz3EfuJ0Co*)
(*v8/SCDcR8jfsPaJtI8CXyIoZrbQIEvqbb+eYPmSJjsXngkSIEr5NDs/Vx0w7*)
(*aW9SighrN8aK9+rrw3eFX3bVfhH2zWCE9stnIH5deNW7WhHimA0LXD1mID3s*)
(*ZN7FByKYL+E8cOuejt7+9sJZ3xTv+9Ls/CKnIzt+akGWwAYm95hvN37gI0+p*)
(*KfeQtw2OGPDd6Cw++KtPt2ZG2+C6PDjKUJuPpPP0c7+9Nrh9ztC4/IQeuJpO*)
(*AvOrNvjBYrp7uuqh71+dArtXNugp5Vi8u6OLuwvT65dOs0VJ6TG5S4Qu3A65*)
(*Sl9JbKH6ZLl68ZAOnut9jUS6LY5usq7XStcBjrSdXl5nixlphRpzuToQd3GM*)
(*Xry3hXe9jmXo0WmY6P/73FsbOwTGWL4xt5sGx9GAjDUpdlAPjt7EwlS8+vGi*)
(*Or/ZDsOHr0gn/Kfi54xKUz7fHnV+h18c6edB9XggeZhoD6FLuf2CJB5Miu61*)
(*OnTZY/s+R9NDyjzsteBmn3FxANPnv+vYJdqIkxcaBlc7YK5zXf5/zLSx0DV0*)
(*5rDVPHgVTZOUXeOioqTCtO38PEyyLb13ZwEXBekfpQ3BjhgQis9p90xBk+45*)
(*XoyKE2pzV1kxEqag+Ipz2NeLTpgjvPDbj+9aUJ4SM/yxQIw8UVFMSYEWHhhd*)
(*f9oV7wyt/btkISZasOg/ZJUmdUFG0tJdHms1sTHLsm7Gb/Mhzzn9WbB3MvQS*)
(*eqt9q1zhZOA/NuWMBvx7jesbVdyRn6vuFtipjmsdJH+g3h3Wp0Zt7ympY06G*)
(*hd6ebA8IfUSt0ulqsGAZNfNTPfH+YtZBZwkHB0rTzX9N9oL6KTOOtVQVXn6C*)
(*gwZl3kgN+DpHmsqGhXqyct8bCUL/TeCrpbCR3/b6/pv3EoirX1GX1rORkXby*)
(*+OfPEoyxuvt4cWwELWNpqE1IsK2+tvhuGBs4by9116aRY7le2cWbjbqs7H9k*)
(*zjSiewdHBt3Z4H2Rr7jqTsM3N+JZ0Xw2jEcEW1okNNRfLr78yYGNkda8e10B*)
(*NAoOWydVC9mYy03+pSeKxgG1oV6+JhvfIqy1dHbS2HRpZestNTY6I1MH+wpo*)
(*hKx5XJfCZiPVYd+hqn00+I13fm+fpLjf1dHNPUKjdPNpkveZAs/cPjZGRmOL*)
(*2Sx7xxEKOlE9H57X0gjrPjjr1RCFleaNZMElGgY2uf96DFC4WXGFenedRsXg*)
(*Gvm3h4p+ycqRiHYamUXPjp3opjCUo0cCummslvy8O6iDQpNM1mD6iIZZGR39*)
(*510KC7bnzy56QUM1oD54RSsFU+WZ8jkvabwZs/fWaKEQNtF25vRrGtUhxjOj*)
(*GyhYH3a5I31PI59Vqq5zlYJ3Rxh3z0cacbW8700XKQxsmqx+9h8aVlrKXQY1*)
(*FCaGKdcDP2jU3OzL66mgENNyIceORVB4ovcv1zIK/O5VGffZBAnZndzKUgpk*)
(*eCzOX43AgrQcWneAQmBm3drPmgQc04anD/ZQYHUw+wRcgtesi7NcCygsM9it*)
(*bc0jqGg5dZqTQ0H8eUykrEuwver4cNIOCkH3e1l39AjCs0psHmyjEMrYGJY2*)
(*nUCf3nW5PJXCt7wbNYcFBGPGWeOqKRQMPZ7f0DQgeKSy1TNpPQWNAbdfkwwJ*)
(*DjQn3XaJpTDCT+gfNyLYUBk9uXwNhae5hvusjIlizshA1ZUUTMbCngSYEGhK*)
(*lj7qWk5h+vy03hgzguHZ/gKXXxT52Qb5rhYS3FWmI48HUYhNbmwKNifIaZr3*)
(*JmERheS8m316cwhWV8y16vKl0Dr+4/2Qwt47hOudCYUSzqzHlywJJnlP/065*)
(*U5D+YV/nbE2wuXFctUyk+K/1gdkXEUFI+Vd/yopC5pZnlsU2BA4ZH36PN6dQ*)
(*ufjWJ3tbgk+efXpiI0W+wkL91XYEcTeahR08CmmB2zXKHAjCrucvUBllwa6T*)
(*aNSICTwb+LnR9SzMmyxczfckWMtVYmhcYUHcGVWVqfCuVS83yy6wsHmeS+0H*)
(*hXvUZHFfZSzsT2L6tngp9hsiWZxZzsK2hILGRAnBwZFEnbIcFi7Ev3Sr9yH4*)
(*YNRc2RPIwrLyik1LAgiO5CTeIAMqyLoWxkldTuDLHjkbzlNBfdfRcp8UAn+q*)
(*wanAWhkuL/ecCdmj2GfPKSW1UCY2WATN3X+GwEb775Kfgplw7fddxzlLQC/6*)
(*4lAcyITQ3UC+VeGEJs9444VMOC5e6xTzB0GDrPehmwsTnipRi91kBBG5GvIN*)
(*AiZ45WdODMoJKsXrVj3tY+BUwOtR1ysElzcWKRk/YyBriBdUq/Bfsmslsb0M*)
(*/DRlao3wKsEXE077aAcDH9fvTuBdIyDa5a7azQxsCFk2MlhPMPi2U8fnBANK*)
(*6uqC/TcIxk1+1OwuZ6DGQZLNaSTgRhou6jrKgLOu38hWhcU98TsiixkYyrp+*)
(*M6aJYGcTayQ9lwEvE48t7i0Exybm5DdmMnB29MnQeYVrxUtMOdsZmFcpCbO4*)
(*SfBYdlRatJmBKqMNHtNuKfJ+2zz6JIUBXS9x7U6FlU3f7Z2dzMCL9lbTSa0E*)
(*upFc69hEBmJTZx1OUdiy1PH2n7EMmD31mvJO4f8Bc6ULKw==*)
(*"]]}, "Charting`Private`Tag$154708#2"]}, {}, {}}, AspectRatio -> GoldenRatio^(-1), Axes -> {True, True}, AxesLabel -> {None, None}, AxesOrigin -> {-Log[1000], -Log[1000000000000000000000000000000]}, Background -> GrayLevel[1], DisplayFunction -> Identity, Frame -> {{True, True}, {True, True}}, FrameLabel -> {{"\!\(\*OverscriptBox[\(\[ScriptCapitalM]\), \(_\)]\)", None}, {"T [keV]", None}}, FrameStyle -> Directive[GrayLevel[0], 20], FrameTicks -> {{{{0, Superscript[10, 0]}, {-Log[100000], Superscript[10, -5]}, {-Log[10000000000], Superscript[10, -10]}, {-Log[1000000000000000], Superscript[10, -15]}, {-Log[100000000000000000000], Superscript[10, -20]}, {-Log[10000000000000000000000000], Superscript[10, -25]}, {-Log[1000000000000000000000000000000], Superscript[10, -30]}}, {}}, {{{-Log[1000], Superscript[10, 3]}, {-Log[100], Superscript[10, 2]}, {-Log[10], Superscript[10, 1]}}, {}}}, GridLines -> {{-Log[10], -Log[20], -Log[30], -Log[40], -Log[50], -Log[60], -Log[70], -Log[80], -Log[90], -Log[100], -Log[200], -Log[300], -Log[400], -Log[500], -Log[600], -Log[700], -Log[800], -Log[900], -Log[1000]}, {Log[10], 0, -Log[10], -Log[100], -Log[1000], -Log[10000], -Log[100000], -Log[1000000], -Log[10000000], -Log[100000000], -Log[1000000000], -Log[10000000000], -Log[100000000000], -Log[1000000000000], -Log[10000000000000], -Log[100000000000000], -Log[1000000000000000], -Log[10000000000000000], -Log[100000000000000000], -Log[1000000000000000000], -Log[10000000000000000000], -Log[100000000000000000000], -Log[1000000000000000000000], -Log[10000000000000000000000], -Log[100000000000000000000000], -Log[1000000000000000000000000], -Log[10000000000000000000000000], -Log[100000000000000000000000000], -Log[1000000000000000000000000000], -Log[10000000000000000000000000000], -Log[100000000000000000000000000000], -Log[1000000000000000000000000000000], -Log[10000000000000000000000000000000], -Log[100000000000000000000000000000000], -Log[1000000000000000000000000000000000]}}, GridLinesStyle -> Directive[Line, RGBColor[0.9, 0.9, 0.9]], ImagePadding -> All, ImageSize -> {553., Automatic}, LabelStyle -> Directive[GrayLevel[0], 14, FontFamily -> "Times"], Method -> {"DefaultBoundaryStyle" -> Automatic, "DefaultMeshStyle" -> AbsolutePointSize[6], "ScalingFunctions" -> None, "CoordinatesToolOptions" -> {"DisplayFunction" -> ({((ConditionalExpression[E^(-#), Inequality[-Pi, LessEqual, Im[#], Less, Pi]]& )[#]& )[Part[#, 1]], (Exp[#]& )[Part[#, 2]]}& ), "CopiedValueFunction" -> ({((ConditionalExpression[E^(-#), Inequality[-Pi, LessEqual, Im[#], Less, Pi]]& )[#]& )[Part[#, 1]], (Exp[#]& )[Part[#, 2]]}& )}}, PlotRange -> {{-Log[10], -Log[1000]}, {-Log[1000000000000000000000000000000], Log[10]}}, PlotRangeClipping -> True, PlotRangePadding -> {{0, 0}, {0, 0}}, Ticks -> {Charting`ScaledTicks[{-Log[#]& , ConditionalExpression[E^(-#), Inequality[-Pi, LessEqual, Im[#], Less, Pi]]& }], Charting`ScaledTicks[{Log, Exp}]}]*)


(* ::InheritFromParent:: *)
(*Graphics[{{{}, {}, Annotation[{RGBColor[0, 0, 1], AbsoluteThickness[1.6], Opacity[1.], Line[CompressedData["*)
(*1:eJwVkns4lIsaxZnRfMOMYtjGJUIGk9xySdhCzDejKFTbdlIoWyKXJFG7U6Hc*)
(*KQq5kyK7J0Uq1VjlUil1XBKRzqGbKMR2jz3nj/d5n/XPen7vWq+Gb4ibH0VM*)
(*TCxGNP/f7SX6P9z81fCkheF6OEIL3m3p0+V+anjzVHp9caQWvi9NLCzuVcPF*)
(*pAaLl1FakPSqIyq91fCUH5qme0ILtmxSleqphqDHHYn9Z7RwI8mbX71FDRkL*)
(*Eo+3Zmsh+XBGgdw6NYxkfjtM3teCI2/WqXNBFfWnMmRSxTnYyx0O905QhbTe*)
(*j5fncznYHxFomaGqCu9dwf/O0NdG5kqK52D1SvSvsAqdEmojiFovsLFbiW6n*)
(*5QMH9+hgjbvUm41dKigJ4LTSZnSQOXDMdIuvCuRvBK56VaQL7VbqcMSYMqJG*)
(*nMxcnbgoUFfeyDujDMv5qWiDES4eVe/cqyGnjLbJ6uD1RWvQckODU3pVCaOn*)
(*bnVNb9HDPI1qa2+jBJmu7HbuvB6686X0Rp4rQms+01BXuBa5+UXV1j6KcJjt*)
(*8XCN1ofku13M7G9sNG+WXnrubIDCSEOhzHE2ZG+fV+zlGGLlsXRpIxYb0z5B*)
(*0XXThnAQsvU9CxXQGx8IHDKC2wH9r2tMFZAa+cBY/LMRmDsDImn4BW8HWxZ9*)
(*Ao3xPe++15LzL0hy8ag5P2mM2i15/ysYlMd9o2vND4+tA9e61EwQKo9qhcdC*)
(*KtsEpzMtdC5JyOPyV5dB9k0TUPkfwui5cqjY6P2j2NUURla1yad05XCcEe2V*)
(*wzDDpiwFx+KHLNhP81kRTWYQN8lvfS5goSlMPfBdvjk+cy1vyHXLYlV464J5*)
(*zHrUJPgZUIJlsSKu4ehEgAXWcu/8OT8ng/LZ2xaFRzcg0TjrQG6aDCLD+8yv*)
(*JFpC5kJKlYe2DJyfltkkXbdCTOj2FLv9K+A/o3N8ROJXVMdfm1TLWI6c+syb*)
(*1r422KDuvCD7lzT4eZlz0mK2SE5gbnTrZOKtjebIjQ5bGFbMmLSKMbFDGNa0*)
(*u9YOXL7xMy8VBjoYbQZ6jfYYvXsmx8pRCgraoQET3ZvArNCVMvSSBCtnP523*)
(*2hFRW6fXekXRwd2RUGKawIPnz2BlxhE6FmsVi66k8mBZ+Ym4d4iOrbJ7ypQy*)
(*eVigdQ3IB9FxTn35EUohDyeFNdkvdtPB5o5l9dXwEK9/SMLagQ6XeuHni//l*)
(*IaB3aHzIlg6f6XuJjE88OCX4vM/6lQ6VNK+BE8M8MD9uq5swp8N4OLLMe4qH*)
(*tDzD0EouHSa/X1elM0lcZHzrVV5Bh1dNCdLNSUTe2/fsKYOOZ5pkT5MVCQ//*)
(*vtojdDqY+VotP2xJKDc8P9cuTkeQmdpZ480k8qOvkYmTBGJ9mzUM9pA4oatp*)
(*ZjFOwPLk/JD0PhK7u3I0P30j8OB0f8z7/STU1yX8tPtMwEtqtmX7IRKXh/yr*)
(*Z3sInPOqeBISSyIu633R1S4CPzlqe5vjSfzh+Fvqjg4CO417NOgpJHSLeQE3*)
(*XxDIPaYn73SBRKUHZ1VAvci/5c+3ty+TSKblM9kPCPi2nz7gXk4iqEZ+rvEu*)
(*AaZq1PU3lSQMZCReq98iUCnkhmfcInHryUBi92UCra++D3UJSaRf7X1pU0xA*)
(*s0xS9sMjEsFnO1ll+QSoR07pdDaS0CObL4VdJLDN/bcTe1pISOnU9785T+D3*)
(*2fechRckvtDuatqkEVhtp3jmxCvRfc0V16TiRfw/rKgWnSROXyn5HhpLiHLV*)
(*oId1kfA+k7vuzUkCo5Em6RndJFR5KXWlUQQa1j/j5fSJ+mkKbbEOJBBg7mHX*)
(*OEjicFnA8lJ/AmsrCzl/fCThFufrJrmPwAbrgZLZTyRWOG5/+3oXAQ1zT6nZ*)
(*IRLxjeu/BrsQMNpUNBQ6KsrzspHBaycC/qVf9e+MkXCI5R6yIgm451xw/Xuc*)
(*hLiDyhxhS+BknHjevyZJRDcsShYbE9jhEDvZPyP6h9JpZ8KAQJ+younELAnz*)
(*mLFzB9cQqGibM6XOk5iwH1CyXE0gadvSa/ZPUf6Pm7gd8iLehsR7KuJ87H6U*)
(*LFg2Q4NkQXeNkM6Hfb1yQoCQBuz+sH0zm4/9LDGK9H0a3OY8eBqKfKT4fYyu*)
(*ukNDYHPTsmmR7mZUBU1X0VBh2LxUrMxHsIfjtrhSGmJ7xHlTqnzkjIewi+Np*)
(*iJNKDrqkxcfY6qaybjcaGkJlrXqN+SiID3lMfl6GZgennx5b+HCij1/3ll+G*)
(*vHu+7jtP8uFM1G9IM5RA5IxEZPtdPlS7K8QYnlTIWHPasr/xsU7uP7muO6lo*)
(*k3l0KfY7HzyXKfNsNyoUxM/6hY6KeBrtD3I2UxE4mr5IjvNRX9Xbs9GaisL4*)
(*eJvpST58EqSrD6tRsXAh++8dC3yUWYb59Q9Q0DhZbi3PEKAuIkuM856CRMNz*)
(*a8SYArysepgb2EtB2ZcMpRGRntKWap/poEDX4cNcw3IBSLlSG7kmCoSuf7WE*)
(*swQYGu5k869SQDHzLOhQEmBRe/5WaikFOjMV+fXKArB8NVxeF1LAOsYsqFQR*)
(*wLL7YKxvtsh/UaIkRlWApEba+PEECgJZ9bUmGgIULa1NboijoC41Q7hKU4Aa*)
(*S3cdqdMUDEsdfcJYLUBfVaFXVjQFRiEBfYNaAowNN828O0JBk3TEl1ccASR0*)
(*RjK0wimYL0/5+762AIq+LMPAEAoydGuo5ToC6OdbtNwMpKDA/CMrU1eAfwAk*)
(*NfA/*)
(*"]]}, "Charting`Private`Tag$152783#1"], Annotation[{GrayLevel[0], AbsoluteThickness[1.6], Opacity[1.], Dashing[{Small, Small}], Line[CompressedData["*)
(*1:eJwVV3k8Vfkb5sq91shaOCdE2Ze0yHaQPfdkKWsilIw1SyjayFLKtEwqW9xk*)
(*nUlpWph4Ck0qGiWhZFLSqsSMFqPf9/fX/Tyf532fdznv+55zNULjvDZzhISE*)
(*NISFhP7/+6DC8LNXBI20iJGvbmIMQnp+nq7eTCOrVnFzqwSD8R+TM7NhNDSr*)
(*BWeUZBmkG/ly1ofSkKoI2TQgzUA8qIlXF0Ij6mXMrnzCnzhISQsH00g+e57T*)
(*KMlAq2mPnG8QDctj0fPlxRlcfP1C+ddAGrtuzLQ+JfFslZ0pkQAa4wvK7AVS*)
(*DLodazX9/Wj4q467FBD/DUnSOud9aKzem3RIk8R7WxFvKLqehiDgu2L7PAap*)
(*PQ+XBnrTQIuyzDXCc4VWml/wpLFsi0wHn+gdNzptzfOgkfu7rNN9eQaaQf/Z*)
(*B7E0IvecmPWUY3D+YIhLozuNFSqKllcJtmpq44uvoTEwGvMHn+A7rxd7B7vS*)
(*COI89PWSYeCnfMDvd2cahy2Cnj8j9b5y/BAk6UTDjKvLhBP7pCSPsE0ONK5K*)
(*aaUHEHuOoHHrFXsadsvipO3mMvi5RylW2o5GYOr86bvKDGihHYlhDA2ZtVa+*)
(*/ooM6o2GUq9Z0/ion3fZgehbBNnukrGi4aMqveqoCoPbBwWZmy1ovDm3c0ZD*)
(*icH6Jm5eszmNsFH3tlSCR15HHp63kkZTt8mafAUG8cpdxyKW07j1RKzuNIk3*)
(*62hy6roZDd3+e/RW0o/8pGOl8ktpvHPPMdci9iqCfwWRJjS6UwJfnCb9re7x*)
(*r2k1ouGSLLA8ReKvELr+m6IhjXkpZdWVJF6bkfqlKH0aodpXVXkUA4+gzGs3*)
(*dGlcmlULuTyfwbODr1qUdWhsStXg7FdjEN3k2h6zmMaE76LMRsJ/e13f2aZF*)
(*Q8+lJ0id+Ocqy95fsIjGea8punoBAyWnxN44DRo9A+8PpKozOJvUN9CxkAa/*)
(*x/bRFU0GpoJVw6o0jUwX/a4pot/SU/xymxoNt9/m2j8j2F1I6O2fKjREv3Mk*)
(*soj/oFHYR2oBwU+kO0RJvK1Bt6YSlWn8fcarNZLo/XtQ91unIo0U6QA9TcJn*)
(*NeX/WKhAI+ee+BtnmoHcm49ztsvRGJLbXVylxaBM2VviniwNh5mNEn0LGRg4*)
(*XZbRlKHhrhNY3Ejsm5IWKKZK09jPSNXEEd5FkK7SLUlDdqw0/YA2g76e4YVa*)
(*EjS05+yL8yZ8uNBq7R1i5PmOyDZUknw+G53T+4tLQyJzwr2TxNsTJG6yWJSG*)
(*mt2iliyC5+ZHL08XoVHn+cH6lS6D4qb7Fg+EaegYhRhyib7em6W2OkLEP712*)
(*uG0Rg6vKJxx3zVIISxptP094R6evbr0zFBbySl2MCP9AumQn9zsF73VDVyZJ*)
(*vzY+sq03/0ph+Y2U+b2Lyf4Vv3z60zQF9crCkVR9BtvDc6VL/qGgv8Pa6I4e*)
(*mXcDA5v7kxQmPVIS6kh+BZP3Y4U/E/vurU9+EH+15sQys08UrlJPz2QakPna*)
(*p/zX5nEKB1UHF48S/2VuzUIn31NQkveWkTJkcGNesOmdtxSUTcMb35J47AAn*)
(*dOY1hZ6HCRm3SL2DZ84dNRqjENwxcUvMlMGWrW5tIaMUZOJ6fhwl/KTx+OTR*)
(*FxROb70z+5no75k+otXxnILEsff/fddhINW6fP30MIWKNHueG7E/mT2wX/cZ*)
(*BavxA2eciJ42m3E58CmFQZVrHjUkn4uKGmOHBikwrz0rnxLeZqhdGf0Udiz1*)
(*t/AzI/fi7FaXz30URDcIBsZMGPhGS6VpPaKQlqPvU7OUwUuzhhqfhxQmliW2*)
(*bCD28d+9B3N7KGgWHWp1N2cwc3Naovk+habqJ9u1ljHIO1Bk+aGLwqlpr0uy*)
(*RE/Ji4leeI9CgPCcZ+uNGFQseFHseYdCwvsLDwuJntHz7K7M20T/UXkYZxWD*)
(*5mq92d9vUfDide42Ifou8d1Gr9sp2Jb+nrWO6PeuTAhWaaOQtK3NxnsFgzDd*)
(*d4kheRT8bn9Q93ch91NPlLmQQ8H8J+MNIQS36y+U4GQTvRWDtWvdGEwZrHrk*)
(*lUXhYtK3Mlsncu+NvM8I9lEYzgz+OET4dcYxUVN7KNwYUTves4bsk0nOCsfd*)
(*xN/E5PFzZwaXTMuFT2RQoHRK21SJ/sulzfde7aQg4BoHaxN7hWWPClfuoFDl*)
(*pLjalmCH5R9Dc1Mp7OoVWjfhSu7tCnGjge0UMq3bTzcR/uzKRV91kymYmfs/*)
(*mCH6vebW7TsSKVxxdn0tyWcwx8K34O42Co2MVEWTB5k3y20BavEU5njcdvdn*)
(*yf5ZHdSOiaWwxf9d9KA7qd+68tP1aAoNLyZ6rIh/u01r89woCsXyGa7xxH+K*)
(*GcgOjqRQvzu7yI5gLbtJz4YIChcW9WUeInidvTQlvIWC+2VBfOlaUv/qJa89*)
(*w8k81nh3/OVN6newa6wIpZD71OZDhheDUcfAXZMhZH/efBu7QfyVnJNdHYIp*)
(*vCi1iw0g9k4uBQq/BJH5MNqs+png7a41w6OBFAzPuxRH+jA459ZWuyKAQv6y*)
(*zs7rnuS+rBlKzvGj4L/sxkFrX/J+5E/b9vtQ6J81fWGxjtxzdp6U7npSX43t*)
(*9A6Ct6zVf5zmTeHtc15Msj95n3s4VtzxJPOmpq6hRvAtz+AYVQ9yD7giwjnr*)
(*yb30SjOPZilwuVHtHwhevO6YyHV3CqN630tWkXx81v/aLb2GzK/l2fA1fgyy*)
(*ff48tdGVgsv2uuT5AQwu+z4PP+9MYWD79QZ7gl/5fTcWcqJwaPrngV6ClQIU*)
(*v3s4kH5o2epZkfhOgca3yu1J/zOWO+UHkfo3uB75bEvhjzHmxkfCnwsK27Ca*)
(*oRBzTK5g1SZS/8aMJcetKXwKM7Hct4HUH1L4+aUlheRR7UT9YFL/pgvXl1tQ*)
(*8FTYsuRfgiNC7+Zmm5P7E2j5cojonQwb9X68ggJ/1UuBJvG/Hf6D1llO+v/f*)
(*D+dZYv9l84K3qWYUStseVA+TeDoRZr93mlI4MqC/02oj+V7Yyt+jYkJhWlKi*)
(*8xjBuZERa6KMiJ5r7Yp3IeQe/7RX6Q8DCnV7lrHJRO91VNFzKX1yb6rK6h+G*)
(*M5gf83t9kC6FjvTpHxeJvUvs/ZTfllAoD/k6dGoz+f6Je2P/Q5uCz5awvRpb*)
(*yD2NF5nroUXy8Thy6BGx799GDZzRJPsXG566P4KBWOLKsxPq5P7G+kfGEn3z*)
(*JM84+4UUxvqe0AXEf2tylMUxiuxHr6i8OdE/rsYJeNGohi0O37Ju7yfvb5FW*)
(*Vxs7NbTGplnnVZH3i7fEY6ZPFSX7ij3F7hH7kZ3L3ENVMbyqfp7hRzIPXSLv*)
(*kj+poKri8Lk3ErYoVVdhnLJV0PXBIl2PssWNRp8wDXkVpFhkNwyb2uLOeQ1t*)
(*QdUCxB8webDLyRbfuSK29jYL4BoYM7w90Bb9JRL67+/OR8WFJTaHE21RVHKm*)
(*0WrTfHjuCtYLOmgL8aENUic/KKMosVai9ZwtylKMW2TTlfFg9s/E4222UNv5*)
(*s7SJnDKOvZMP8PzbFg4tyoYBZUpI7PafeiZsB6+fDN/qLVOCtYfIP6cW2UHK*)
(*JzKFC0WwEfPStdfYYby4OegHXxEj8V6f4pLtcNm9+HnpCwVkJacZe56xg66V*)
(*YLlrvAJGYjP2zD6ww77j5ktOz1FA9diM91WuPURcXm4TK5IH9nNjpW3sYWJ5*)
(*OX+vjjx8px40Dey2x+pCJcfy63KQXWKr79ZmD2Gzkq67rnKI2qxTcFN+NcZ0*)
(*Lc7L98/DEd+c9fOjVuNS3mYjTuw8fKm+9GGofTUMdK9kfP8mC72ze5wCTRxw*)
(*wLTwp6ICWQwVVaXzyxwg+8uhBr/FsljpdyZAXdkRmfHrDtltlcHRS5eEJ845*)
(*ojG3doo+Nhcns8b1cy2dsEqdPzOvXhp/J6y2XzzuhPw8KcarVwpH3yY4ulU5*)
(*w7jmi1mXkBRcpmzqp+JcoOti2hmkKokrHc+cMja54uPV7FOWjhJ4mh4yPhjo*)
(*BqkaHQnjIHFkJ9SbZ+atQdraaYOgNDHwihWaHR66I+C/WBXJ7WLQ/q6bwe93*)
(*h0XdK961BDFYuDQl+g65Y4bbN6IQLQbJU0q/xIy5Y0/LpZP3NophWOPtMsF3*)
(*d+QaJsyxchDDcU5elrwWH5FP3ky8sRXDM+bmZ1qXD7e8TcOF1mJYGR6rqWvE*)
(*h9SoR9PkCjEkpI2Vm5vzUVBsHF+nK4bKy6a9Fu58nJD88ERFRgzl2TfvpCXy*)
(*kXItvPO2pBj+sTxc6pjKh1/E08vbxUh+7mme4hl8qLTdPfJAWAxJo0U5ydl8*)
(*lOyodT4wxcOvl0MNuk7ysUtHc7n5BA/2kV82ri7hY2PfKc1XH3hQsUqcU1/O*)
(*h/rSvP/sxnh4UmXcx9bycfZNROPXAR7irSv2tjfxsb9w+ExVHw+rPv42daGF*)
(*jy2OvofXP+QhbepixoGbfOiUO0VeuMeDQLyq72snH+JrW3yCO3lI6NCKKuji*)
(*4+3McgfpW4R/v8hHuoePOj/thZGtPNxopvl4zEc+t0RK+Q8evHf13Joa5CP6*)
(*ksK39qs8fG5svy7zjA8j2TmP1C/yULh/qPvLCz4u/jlyoP8sD36bPJcUjvPx*)
(*c9WTbptyHlL+9hJdOsFHbE6vXGUJDxf6gyUbJvnQd751etsJHpqvdIn5fuFD*)
(*Yknrs8dHeWD1k8r2fuPjNfeqpk0Bj/zf2HP/yAyp/1ZNrUQuD+0fy15sEWKx*)
(*71zFeHwWD5e+NnQYcliEZBctfbyHhzH5Kt4zERaU06EmQRrJN7wl5geXxYx2*)
(*9qz4dh5E/zBwSRBjMSi62z4+gYe7UWvVH4qzONERf8cqiofb1bktm6RYJFVG*)
(*zhVE8KC996PuEWkWXvtDvcTDeTCTEE28OJeFjOO6wUcbeKBX7SzskmUxrsWn*)
(*rfxJPinanM55LO7NcQqtWM9DuPI/66/JschtX/k2luUh6mlKTroCiy1nTYwe*)
(*ufGwiRd21UORhUOWboKlMw+/1VlUqymxEHZQ/caz5aG0IMeiVJnFjrZZ8XJT*)
(*HtzCFIf9VFj4Cab5PCMSP2Lg2jeCV2R+OhKjx8PBv7urTqqymLQfWWCxiIe/*)
(*JkuabquxiL7ZoftQgYfoLNEakYUsNt7IdxX9woX0U/cTCotY2Leq5EW2cPEu*)
(*4jadocdiq5wQR7qZixLZ9vDXBB/aPLqj4QoXa28E7vTWZ9Ev2RA93cDFAcVd*)
(*MzoGLGL9HD32C7j48Vwz7oEhi1MTccrluVw4jWh3q5uy+LSoo7Lfi4v7n/oX*)
(*Xl7JojQ37qbzmChaqmsqnziwcBOb+DVEQRQbh5VC6zay4PNaVxUYz0H6H35B*)
(*zTvI8+yvEZIMEEFJb9/g9HEWS+X/KvL0EYHG0b9FYn5h4cT+u+Kklwi+dA2l*)
(*jBAc224fo71GhOxRw9i9EyxaG54MMFYicB65OFN+ksWmPOnGJFoEE/IG4e7F*)
(*LCottm1+NsJBUV1GfbmARVNyoZD2MAf7an3/UzrLorvhelHUEw60q3+OyCf4*)
(*38USD7485KDURopOqWThLC+wke/gIKfw8D73KhZv3vUqu1RxoE/F3ZquZTG7*)
(*+PvFwwIOIt8kmcXUsZAL1WAflXEQu9v+3AjBFv0xWaEnOfBUc/ipq57FwXbu*)
(*RHoeB9a/5CsLfmNx5odBftt+DqYfHs+Zf57FJQvvJRL7OLjeeU3hMMFPG8qC*)
(*CndwwDU3G01tIP1+1/FlaDsHj7cqh40TPGfJ+2NaiRxU+r67H3aBxfxQOeOo*)
(*OA5s1DjmAwQblpjfuRDFgeqdf6rYiyz+B0L/aSM=*)
(*"]]}, "Charting`Private`Tag$152783#2"], Annotation[{GrayLevel[0], AbsoluteThickness[1.6], Opacity[1.], Dashing[{Small, Small}], Line[CompressedData["*)
(*1:eJwV0nk41NsfB3DNMN+x3YSMdbImypbSGJN9mzO2VK6EZKm02CouaVG4FKls*)
(*iVHTKEmLpauFeHdt0XJvCymVX/qVKKJNRfWb3x/nOc/rOZ/nnM953h+98Fj/*)
(*KJqUlFS2ZP1/v3/S7IP/ejYeG2QFR+hyEHbv0OSZKDZ+C7xQ/VjisV8fp39G*)
(*sFERwOvw0+NANuQaVR3Gxvvg4U08fQ4cWR469CA2ahbxniobcnDxQJhnvRcb*)
(*K/8zo7zNmIOcbfnlKgvZSGvsfGlqwYGb+zfycFoHr8983jDLnoMIk7dbw7J1*)
(*QC0d/OdBKAcbtm/i5uvooDCkfYC5i4MCbVrQy3ptbGe6LfcVcrCZ3sK3d9LG*)
(*48KmL5rgwHS53COHXi0MFHxtzXghqR/cscgrXAsGJbFmwwxbzL1Df7t9XBOn*)
(*q5t9fiywRbmupoN7piYeB9RvfbPMFjfqAyL0VDSxW+SWcHGHLbov6hmJKzUg*)
(*rqi+mVNhiykG3dHZXgNT1vusiu/aok8oN//dLXW4nquLuvzTFqXCE/W8teoY*)
(*eFFK8c25kH0WrHB0lIWBmYRpEMXF8SSLZqVUFszZif4t5Vxo7zikaKnMQrSI*)
(*3uXSx4VrM8ss6LgaBsN3rvLRsoP/RrMR00VqGPqgP7871A4KAdFJDMxG+9gN*)
(*v/fn7TBW1hjyy3s2uD4H7qTSeGjwKntR/lIVan77xuwCeTDhiRfz41Tx9cjB*)
(*uJLrPOwt4Bgfk1bFoNwFs5lGS0H3/G88s1QFRfX0m1fyl8LSriEnbZ4KGhc1*)
(*tfux7OFSrOYmuq4Mx/NH7liX22OGtfDOLb4yxM1F0voODhgy4V5U6ZuFzqSc*)
(*xtujDriUHWVOi5kFLQUv2nk/Rywwubxz6rsSOi9wuk2bHLHfqnhjaZ4Srjjt*)
(*cNae7wSlwtyawLlKmDokHBk76YR9cStynTbMRMLqwAstC51Rn3X2Ezv/N8wh*)
(*7KqATmfY6npPzzqniPaTvrM+J7sgJ1vBwf+hAirpaQPq1q6wqPpqfUdKAV5v*)
(*fDe3TrnCxNOqK0RLHuJML3WvUTe8v5JZYucmhwrBXzurhtyhUDVPziJEFnVr*)
(*l71bNtsTyb6TC0KSmVi14tzvKmv4CPoRoymfyIT6B20D7Ug+uNWvqasJTFC1*)
(*Q5+NovmYZvQOqm5mwjDCMoy3lY89zZeO3g5losH+j5r4P/nIMkuQ5rkyUd3V*)
(*4Dl4kY/o/uGJYUcmdBsH60Yv8UGy1w4UL2VCcZLcm7zKh8Irv2sfbZhYYd6u*)
(*S7XxkVdmEVdtIjlvFp2T7eOjSH60X3MmEz2XN0R4/uIj6Wpk1015Jni+9TxV*)
(*aYLA9U8bEplMPDvaNvqISaDZeuvw/RlM/KqSM3VVJhCmnPXY/4lCsW9Mc6kh*)
(*wa55+os5ExTuK/50UTYhCO0t0X89SqEo//OxFDMC3YXZP5yGKEQc/hakY0NQ*)
(*Mby+/ttjCpHJha757gQZxQMnKnspKKf7H9tDCNa5/X5w5QMKz+uGlVf5EMwT*)
(*uUfX3pbUa9xIurmSQNa3OWBNFwXLkkBG6CqCkenFroodFMaWyUg9DSaoDjSa*)
(*E91CISMgnJ0XQZDDECqwmiiEmRVGta4j2HxJ9XvbFQrqhg5rnkcTmCtJ9+jW*)
(*UdC0tjrcGUtQ1zm4v6+CQtPUu/6oFIJDlf137UUUwhWKe6V2EsT8+VD5lJDC*)
(*rqD52LWbYL5Hx7H4IgpZf400GqUTyBm3PH90hMITddvFfpkEbxhX9O3zKKxd*)
(*h3VhWZL/d1SdlcuiYHrcPGNxDsHe0yfH4tIp7FBXCZrOJQjLLF34aA8FQ+Oq*)
(*2nN5BDruudfEyRT2Xv0uc/cIwbRR5k/ZRAovQuXVnAoInsjsdo5LkPS3O+GW*)
(*qJCgqD2um7eJgnP9ORvLowTbTkX/Jl5PYc84f2hNCYF/Rri/bCQF5OVYph0j*)
(*mOm24klPMIVnrZfMy8oIxgy92bxVkvuncoXFQoLb0u7hJ1dK8sudYGWVE2S1*)
(*LRmJ8aFQ0/5tseCEJK8KS/MeQmG5n43LHBGBa7pJgp0HhQMyYRZvJZ7hqvWd*)
(*cqRg8zlCdouYIKX1p6zIisIUSzHV+LRk3sST3pQ5hfoljLQuiW32jR/eYkph*)
(*44Np9/WVBB+dBzW4BhTOx7TEFZ+R5Pd3u8kDVQrJbTVh989K5u1GDl/mKwMH*)
(*RElhLy8QOLdoZkc3M5B6MviqoIFgg7IUTbGRAZPakm8NEudGvUqpucxApOgP*)
(*Sv8yQZ98zebJGgbS8zjCSYljAt38MsQMHPa6qSa6SlAyEcsSZTGQ+G5W/vsm*)
(*gnGD9lN9/gzMs+nvTm8lKM+K/dtjSAbPFxWEnP2XgDAnzoepyoCl0fMleYTA*)
(*m2qxzbOQRu+97l/ZDAF0+qqk5IPoONZ9cF+HrgALVf4tXRZAh4tpapudngDu*)
(*Pl9sjvrTkekXwqiVOKbNeYuRgI66HydySvUFaKnpf+zAo+NHb4k4zlCAtdmK*)
(*9dvYdHiPv/quOU+AU9z4qOeDNJARXa84CwGubS+WMhqgoXu115lXEt+tuV66*)
(*qZ+GUh8is9pSgC9z5e5/fUDDmoutN9ysBPBQEdurtNOQ7nPAVctagOG3D1me*)
(*lTQkJTjHd9gI8HPuVN1BMQ0WH3c+sVsigHK4nk/PcRpSPTNcayXm9m1JDz9K*)
(*Q/TIF80yjgAH2hgTqdk0dLxJvBfPFeDErwU5rRk0fPpnptOQxJe4y43l9tLw*)
(*JnhnbbCdAE9rjocUp9Dw0u5sgQdPgPG37V+fJdIQOhhPXZdY2vhdvuFWGnI1*)
(*6DsWLhVAPVzZYlMsDWnflr+vlNhMyOmu3SR5z3d7pI69AP8D89IS/Q==*)
(*"]]}, "Charting`Private`Tag$152783#2"], Annotation[{GrayLevel[0], AbsoluteThickness[1.6], Opacity[1.], Dashing[{Small, Small}], Line[CompressedData["*)
(*1:eJwV0Xk41GsbB3DNML8xiMYUpkwkyxAaWwzZ5/dIpKM4jhpZ2uxK6YhTvcma*)
(*4rRRUsnSok5zjLTzLUtpOR1rSitvSUWpTqkc3nn/eK77+lzPcz/XfX1vw8jE*)
(*wFUMJSWlXMX5f20/bvkxcI0APYOaHUqOIoS3FX49uUqAx/LXlSkKD098GhuP*)
(*EqBwYGzrO4VVpZep6nABfFWrfHqcRHDXIfrMUAGWW7ZZyZxFOLcz3EfuJ0Co*)
(*v8/SCDcR8jfsPaJtI8CXyIoZrbQIEvqbb+eYPmSJjsXngkSIEr5NDs/Vx0w7*)
(*aW9SighrN8aK9+rrw3eFX3bVfhH2zWCE9stnIH5deNW7WhHimA0LXD1mID3s*)
(*ZN7FByKYL+E8cOuejt7+9sJZ3xTv+9Ls/CKnIzt+akGWwAYm95hvN37gI0+p*)
(*KfeQtw2OGPDd6Cw++KtPt2ZG2+C6PDjKUJuPpPP0c7+9Nrh9ztC4/IQeuJpO*)
(*AvOrNvjBYrp7uuqh71+dArtXNugp5Vi8u6OLuwvT65dOs0VJ6TG5S4Qu3A65*)
(*Sl9JbKH6ZLl68ZAOnut9jUS6LY5usq7XStcBjrSdXl5nixlphRpzuToQd3GM*)
(*Xry3hXe9jmXo0WmY6P/73FsbOwTGWL4xt5sGx9GAjDUpdlAPjt7EwlS8+vGi*)
(*Or/ZDsOHr0gn/Kfi54xKUz7fHnV+h18c6edB9XggeZhoD6FLuf2CJB5Miu61*)
(*OnTZY/s+R9NDyjzsteBmn3FxANPnv+vYJdqIkxcaBlc7YK5zXf5/zLSx0DV0*)
(*5rDVPHgVTZOUXeOioqTCtO38PEyyLb13ZwEXBekfpQ3BjhgQis9p90xBk+45*)
(*XoyKE2pzV1kxEqag+Ipz2NeLTpgjvPDbj+9aUJ4SM/yxQIw8UVFMSYEWHhhd*)
(*f9oV7wyt/btkISZasOg/ZJUmdUFG0tJdHms1sTHLsm7Gb/Mhzzn9WbB3MvQS*)
(*eqt9q1zhZOA/NuWMBvx7jesbVdyRn6vuFtipjmsdJH+g3h3Wp0Zt7ympY06G*)
(*hd6ebA8IfUSt0ulqsGAZNfNTPfH+YtZBZwkHB0rTzX9N9oL6KTOOtVQVXn6C*)
(*gwZl3kgN+DpHmsqGhXqyct8bCUL/TeCrpbCR3/b6/pv3EoirX1GX1rORkXby*)
(*+OfPEoyxuvt4cWwELWNpqE1IsK2+tvhuGBs4by9116aRY7le2cWbjbqs7H9k*)
(*zjSiewdHBt3Z4H2Rr7jqTsM3N+JZ0Xw2jEcEW1okNNRfLr78yYGNkda8e10B*)
(*NAoOWydVC9mYy03+pSeKxgG1oV6+JhvfIqy1dHbS2HRpZestNTY6I1MH+wpo*)
(*hKx5XJfCZiPVYd+hqn00+I13fm+fpLjf1dHNPUKjdPNpkveZAs/cPjZGRmOL*)
(*2Sx7xxEKOlE9H57X0gjrPjjr1RCFleaNZMElGgY2uf96DFC4WXGFenedRsXg*)
(*Gvm3h4p+ycqRiHYamUXPjp3opjCUo0cCummslvy8O6iDQpNM1mD6iIZZGR39*)
(*510KC7bnzy56QUM1oD54RSsFU+WZ8jkvabwZs/fWaKEQNtF25vRrGtUhxjOj*)
(*GyhYH3a5I31PI59Vqq5zlYJ3Rxh3z0cacbW8700XKQxsmqx+9h8aVlrKXQY1*)
(*FCaGKdcDP2jU3OzL66mgENNyIceORVB4ovcv1zIK/O5VGffZBAnZndzKUgpk*)
(*eCzOX43AgrQcWneAQmBm3drPmgQc04anD/ZQYHUw+wRcgtesi7NcCygsM9it*)
(*bc0jqGg5dZqTQ0H8eUykrEuwver4cNIOCkH3e1l39AjCs0psHmyjEMrYGJY2*)
(*nUCf3nW5PJXCt7wbNYcFBGPGWeOqKRQMPZ7f0DQgeKSy1TNpPQWNAbdfkwwJ*)
(*DjQn3XaJpTDCT+gfNyLYUBk9uXwNhae5hvusjIlizshA1ZUUTMbCngSYEGhK*)
(*lj7qWk5h+vy03hgzguHZ/gKXXxT52Qb5rhYS3FWmI48HUYhNbmwKNifIaZr3*)
(*JmERheS8m316cwhWV8y16vKl0Dr+4/2Qwt47hOudCYUSzqzHlywJJnlP/065*)
(*U5D+YV/nbE2wuXFctUyk+K/1gdkXEUFI+Vd/yopC5pZnlsU2BA4ZH36PN6dQ*)
(*ufjWJ3tbgk+efXpiI0W+wkL91XYEcTeahR08CmmB2zXKHAjCrucvUBllwa6T*)
(*aNSICTwb+LnR9SzMmyxczfckWMtVYmhcYUHcGVWVqfCuVS83yy6wsHmeS+0H*)
(*hXvUZHFfZSzsT2L6tngp9hsiWZxZzsK2hILGRAnBwZFEnbIcFi7Ev3Sr9yH4*)
(*YNRc2RPIwrLyik1LAgiO5CTeIAMqyLoWxkldTuDLHjkbzlNBfdfRcp8UAn+q*)
(*wanAWhkuL/ecCdmj2GfPKSW1UCY2WATN3X+GwEb775Kfgplw7fddxzlLQC/6*)
(*4lAcyITQ3UC+VeGEJs9444VMOC5e6xTzB0GDrPehmwsTnipRi91kBBG5GvIN*)
(*AiZ45WdODMoJKsXrVj3tY+BUwOtR1ysElzcWKRk/YyBriBdUq/Bfsmslsb0M*)
(*/DRlao3wKsEXE077aAcDH9fvTuBdIyDa5a7azQxsCFk2MlhPMPi2U8fnBANK*)
(*6uqC/TcIxk1+1OwuZ6DGQZLNaSTgRhou6jrKgLOu38hWhcU98TsiixkYyrp+*)
(*M6aJYGcTayQ9lwEvE48t7i0Exybm5DdmMnB29MnQeYVrxUtMOdsZmFcpCbO4*)
(*SfBYdlRatJmBKqMNHtNuKfJ+2zz6JIUBXS9x7U6FlU3f7Z2dzMCL9lbTSa0E*)
(*upFc69hEBmJTZx1OUdiy1PH2n7EMmD31mvJO4f8Bc6ULKw==*)
(*"]]}, "Charting`Private`Tag$152783#2"]}, {}, {}}, AspectRatio -> GoldenRatio^(-1), Axes -> {True, True}, AxesLabel -> {None, None}, AxesOrigin -> {-Log[1000], -Log[1000000000000000000000000000000]}, Background -> GrayLevel[1], DisplayFunction -> Identity, Frame -> {{True, True}, {True, True}}, FrameLabel -> {{"\!\(\*OverscriptBox[\(\[ScriptCapitalM]\), \(_\)]\)", None}, {"T [keV]", None}}, FrameStyle -> Directive[GrayLevel[0], 20], FrameTicks -> {{{{0, Superscript[10, 0]}, {-Log[100000], Superscript[10, -5]}, {-Log[10000000000], Superscript[10, -10]}, {-Log[1000000000000000], Superscript[10, -15]}, {-Log[100000000000000000000], Superscript[10, -20]}, {-Log[10000000000000000000000000], Superscript[10, -25]}, {-Log[1000000000000000000000000000000], Superscript[10, -30]}}, {}}, {{{-Log[1000], Superscript[10, 3]}, {-Log[100], Superscript[10, 2]}, {-Log[10], Superscript[10, 1]}}, {}}}, GridLines -> {{-Log[10], -Log[20], -Log[30], -Log[40], -Log[50], -Log[60], -Log[70], -Log[80], -Log[90], -Log[100], -Log[200], -Log[300], -Log[400], -Log[500], -Log[600], -Log[700], -Log[800], -Log[900], -Log[1000]}, {Log[10], 0, -Log[10], -Log[100], -Log[1000], -Log[10000], -Log[100000], -Log[1000000], -Log[10000000], -Log[100000000], -Log[1000000000], -Log[10000000000], -Log[100000000000], -Log[1000000000000], -Log[10000000000000], -Log[100000000000000], -Log[1000000000000000], -Log[10000000000000000], -Log[100000000000000000], -Log[1000000000000000000], -Log[10000000000000000000], -Log[100000000000000000000], -Log[1000000000000000000000], -Log[10000000000000000000000], -Log[100000000000000000000000], -Log[1000000000000000000000000], -Log[10000000000000000000000000], -Log[100000000000000000000000000], -Log[1000000000000000000000000000], -Log[10000000000000000000000000000], -Log[100000000000000000000000000000], -Log[1000000000000000000000000000000], -Log[10000000000000000000000000000000], -Log[100000000000000000000000000000000], -Log[1000000000000000000000000000000000]}}, GridLinesStyle -> Directive[Line, RGBColor[0.9, 0.9, 0.9]], ImagePadding -> All, ImageSize -> {554., Automatic}, LabelStyle -> Directive[GrayLevel[0], 14, FontFamily -> "Times"], Method -> {"DefaultBoundaryStyle" -> Automatic, "DefaultMeshStyle" -> AbsolutePointSize[6], "ScalingFunctions" -> None, "CoordinatesToolOptions" -> {"DisplayFunction" -> ({((ConditionalExpression[E^(-#), Inequality[-Pi, LessEqual, Im[#], Less, Pi]]& )[#]& )[Part[#, 1]], (Exp[#]& )[Part[#, 2]]}& ), "CopiedValueFunction" -> ({((ConditionalExpression[E^(-#), Inequality[-Pi, LessEqual, Im[#], Less, Pi]]& )[#]& )[Part[#, 1]], (Exp[#]& )[Part[#, 2]]}& )}}, PlotRange -> {{-Log[10], -Log[1000]}, {-Log[1000000000000000000000000000000], Log[10]}}, PlotRangeClipping -> True, PlotRangePadding -> {{0, 0}, {0, 0}}, Ticks -> {Charting`ScaledTicks[{-Log[#]& , ConditionalExpression[E^(-#), Inequality[-Pi, LessEqual, Im[#], Less, Pi]]& }], Charting`ScaledTicks[{Log, Exp}]}]*)


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
