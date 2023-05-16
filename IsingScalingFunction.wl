BeginPackage["IsingScalingFunction`"]

InverseCoordinates::usage = "Numerically convert Schofield coordinates to t and h."

g::usage = "g[θ0, gs][θ] gives the Schofield coordinate transformation defined in (14)."

ut::usage = "ut[θ] gives the scaling field u_t as a function of Schofield coordinates."

uh::usage = "uh[θ0, gs][θ] gives the scaling field u_h as a function of Schofield coordinates."

η::usage = "η[θ0, gs][θ] gives the invariant scaling combination η."

ξ::usage = "ξ[θ0, gs][θ] gives the invariant scaling combination ξ."

ReScriptF::usage = "ReScriptF[θ0, θYL, B, C0, CYL, Gs][θ] gives the free energy scaling function defined in (19)."

ScriptF::usage = "ScriptF[θ0, θYL, B, C0, CYL, Gs][θ] gives the free energy scaling function defined in (35)."

DScriptFPlusMinusDξθ0List::usage =
  "DScriptFPlusMinusDξθ0List computes the first m derivatives of the scaling function F_- evaluated at θ_0."

DScriptFPlusMinusDξList::usage =
  "DScriptFPlusMinusDξList computes the first m derivatives of the scaling function F_+/-."

DScriptF0DηList::usage =
  "DScriptF0DηList computes the first m derivatives of the scaling function F_0."

DScriptF0Dη::usage =
  "DScriptF0Dη computes the mth derivative of the scaling function F_0."

DScriptMCasDξList::usage = "Computes the first m derivatives of the scaling function M given by Caselle et al."

uf::usage = "uf computes the singular free energy u_f."

DufDut::usage =
  "DufDut computes derivatives of the singular free energy u_f with respect to the temperature-like scaling field u_t."

DufDuh::usage =
  "DufDuh computes derivatives of the singular free energy u_f with respect to the temperature-like scaling field u_h."

ruleB::usage = "Fixes B given other data as in (38)."

ruleC0::usage = "Fixes C0 given other data as in (39)."

Data::usage = "Data[n] gives data from the fit to nth order from Table II."

PrepareArgument::usage = "Converts scaling function data into appropriate argument to function interfaces."

θ0Cas::usage = "θ0 from Caselle et al."

gsCas::usage = "g function coefficients from Caselle et al."

Φs::usage = "List of numeric coefficients for the scaling function F_0"

Gls::usage = "List of numeric coefficients for the scaling function F_-"

Ghs::usage = "List of numeric coefficients for the scaling function F_+"

Begin["Private`"]

β := 1/8

δ := 15

Δ := β δ

OverlineS := 2^(1/12) Exp[-1/8] Glaisher^(3/2)

Φs := {
  -Gamma[1/3]Gamma[1/5]Gamma[7/15]/(2 π Gamma[2/3]Gamma[4/5]Gamma[8/15])(4 π^2 Gamma[13/16]^2 Gamma[3/4]/(Gamma[3/16]^2 Gamma[1/4]))^(8/15),
  -0.31881012489061,
  Around[0.110886196683, 2.0 10^-12],
  Around[0.01642689465,  2.0 10^-11],
  Around[-2.639978 10^-4, 1.0 10^-10],
  Around[-5.140526 10^-4, 1.0 10^-10],
  Around[2.08865 10^-4,   1.0 10^-9],
  Around[-4.4819 10^-5,   1.0 10^-9],
  Around[3.16 10^-7,      1.0 10^-9],
  Around[4.31 10^-6,      0.01 10^-6],
  Around[-1.99 10^-6,     0.01 10^-6]
}

Gls := {
  0,
  -OverlineS,
  −1.000960328725262189480934955172097320572505951770117 Sqrt[2]/((2 )^(-7/8) (2^(3/16)/OverlineS)^2)/2/(12 \[Pi]),
  Around[ 0.038863932, 3.0 10^(-9)],
  Around[−0.068362119, 2.0 10^(-9)],
  Around[ 0.18388370,  1.0 10^(-8)],
  Around[-0.6591714,   1.0 10^(-7)],
  Around[ 2.937665,    3.0 10^(-6)],
  Around[-15.61,       1.0 10^(-2)],
  Around[ 96.76,       1.0 10^(-2)],
  Around[-6.79 10^2,   1.0],
  Around[ 5.34 10^3,    10.],
  Around[-4.66 10^4,   0.01 10^4],
  Around[ 4.46 10^5, 0.01 10^5],
  Around[-4.66 10^6, 0.01 10^6]
}

Ghs := {
  0,
  0,
  -1.000815260440212647119476363047210236937534925597789 Sqrt[2]/((2 )^(-7/8) (2^(3/16)/OverlineS)^2)/2,
  0,
  Around[  8.333711750, 5.0 10^(-9)],
  0,
  Around[-95.16896,     1.0 10^(-5)],
  0,
  Around[1457.62,       3.0 10^(-2)],
  0,
  Around[-2.5891 10^4,  2.0],
  0,
  Around[5.02 10^5, 0.01 10^5],
  0,
  Around[-1.04 10^7, 0.01 10^7]
}

Data[2] = Rationalize[#, 10^-20] & /@ <|
  "θ0"  -> 1.148407773492004`,
  "θYL" -> 0.9896669889911205`,
  "CYL" -> -0.172823989504767`,
  "Gs"  -> {-0.31018352388662596`, 0.2474537923130002`},
  "gs"  -> {0.37369093055254343`, -0.021636313152585823`}
|>

Data[3] = Rationalize[#, 10^-20] & /@ <|
  "θ0"  -> 1.2542120477507488`,
  "θYL" -> 0.6020557328641167`,
  "CYL" -> -0.38566364361428684`,
  "Gs"  -> {-0.3527514794812415`, 0.2582430860166863`},
  "gs"  -> {0.4483788209731592`, -0.022032295172535358`, 0.00022200608228654115`}
|>

Data[4] = Rationalize[#, 10^-20] & /@ <|
  "θ0"  -> 1.3164928721109121`,
  "θYL" -> 0.6400189996493497`,
  "CYL" -> -0.3563974694580203`,
  "Gs"  -> {-0.3550547624920263`, 0.23465947408509413`, -0.0019083731028066697`},
  "gs"  -> {0.4410742751152714`, -0.034817777358116885`, 0.000678172648789985`, -0.00004305140578834467`}
|>

Data[5] = Rationalize[#, 10^-20] & /@ <|
  "θ0"  -> 1.3403205742656135`,
  "θYL" -> 0.6238113973493433`,
  "CYL" -> -0.38002950945224295`,
  "Gs"  -> {-0.35127522582179693`, 0.23704589676915347`, -0.007319731639727028`},
  "gs"  -> {0.44371885415894785`, -0.04609943321005163`, -0.0007458341071947777`, 0.00005966875622885447`, -4.403083529955303`*^-6}
|>

Data[6] = Rationalize[#, 10^-20] & /@ <|
  "θ0"  -> 1.3626103817690176`,
  "θYL" -> 0.6462147447024515`,
  "CYL" -> -0.35576386594103865`,
  "Gs"  -> {-0.3520586281920383`, 0.23316561297622435`, -0.006649030656179257`, -0.0016899077640685814`},
  "gs"  -> {0.43845335615925396`, -0.05312704168994819`, -0.003914782631377569`, -0.0004080160912692574`, 0.000026262906640471588`, -1.0974538440746828`*^-6}
|>

PrepareArgument[data_] := With[
  {
    θ0 = data["θ0"],
    gs = data["gs"]
  },
  {
    θ0,
    data["θYL"],
    ruleB[θ0, gs],
    ruleC0[θ0, gs],
    data["CYL"],
    data["Gs"],
    gs
  }
]

t[θ_] := θ^2 - 1

g[θ0_, gs_][θ_] := (1 - (θ/θ0)^2) Total[MapIndexed[Function[{gi, i}, gi θ^(2*i[[1]]-1)], gs]]

ut[R_, θ_] := R t[θ]

uh[θ0_, gs_][R_, θ_] := R^Δ g[θ0, gs][θ]

InverseCoordinates[\[Theta]0_, gs_, wp_:20][tn_, hn_] :=
 ({Exp[logR], \[Theta]0 Tanh[x]} /.
    FindRoot[{
      Rationalize[tn, 10^-30] == ut[Exp[logR], \[Theta]0 Tanh[x]],
      Rationalize[hn, 10^-30] == uh[\[Theta]0, gs][Exp[logR], \[Theta]0 Tanh[x]]
     }, {{logR, 2}, {x, Sign[hn]/2}}, WorkingPrecision -> wp]) /;
  NumericQ[tn] && NumericQ[hn]

η[θ0_, gs_][θ_] := t[θ] / RealAbs[g[θ0, gs][θ]]^(1 / Δ)

ξ[θ0_, gs_][θ_] := g[θ0, gs][θ] / RealAbs[t[θ]]^Δ

ScriptR[θc_, B_][θ_] := (θc Exp[1/(B θc)] ExpIntegralEi[-1/(B θc)] + (θ - θc) Exp[-1/(B (θ - θc))] ExpIntegralEi[1/(B (θ - θc))]) / π

ReScriptF0[C0_, θc_, B_][θ_] := C0 (ScriptR[θc, B][θ] + ScriptR[θc, B][-θ])

ScriptFYL[θYL_, CYL_][θ_] := CYL ((-I θ + θYL)^(5/6) + (I θ + θYL)^(5/6)  - 2 θYL^(5/6))

ReScriptFRegular[θ0_, θYL_, B_, C0_, CYL_, Gs_][θ_] := C0 ScriptR[θ0, B][-θ] + ScriptFYL[θYL, CYL][θ] + Total[MapIndexed[Function[{G, i}, G θ^(2*i[[1]])], Gs]]

ReScriptF[θ0_, θYL_, B_, C0_, CYL_, Gs_][θ_] := ReScriptFRegular[θ0, θYL, B, C0, CYL, Gs][θ] + C0 ScriptR[θ0, B][θ]

DReScriptFIrregular[θ0_, B_, C0_][m_] := Piecewise[{{C0 m! Gamma[m - 1] B^(m - 1) / π, m > 1}, {C0 θ0 Exp[1/(B θ0)] ExpIntegralEi[-1/(B θ0)] / π, m == 0}}, 0]

ScriptF[θ0_, θYL_, B_, C0_, CYL_, Gs_][θ_] := ReScriptF[θ0, θYL, B, C0, CYL, Gs][θ] + C0 I Sign[Im[θ]] ((θ-θ0)Exp[-1/(B(θ-θ0))]-(-θ-θ0)Exp[-1/(B(-θ-θ0))])

ScriptFPlusMinus[ScriptF_][θ_] := ScriptF[θ] / t[θ]^2 - 1/(8 \[Pi]) Log[t[θ]^2]

ScriptF0[θ0_, gs_][ScriptF_][θ_] := RealAbs[g[θ0, gs][θ]]^(-2 / Δ) ScriptF[θ] - η[θ0, gs][θ]^2 Log[g[θ0, gs][θ]^2] / (8 π Δ)

uf[params__][R_, θ_] := R^2 ReScriptF[params][θ] + t[θ]^2 R^2 / (8 π) Log[R^2]

EfficientDerivativeList[n_][f_][x_] := Module[
  {xp}, NestList[Function[g, D[g, xp]], f[xp], n] /. xp -> x
]

InverseDerivativeList[n_][f_][x_] := Module[
    {xp, dfs, fp, Pns},
  dfs = Rest[EfficientDerivativeList[n][f][x]];
  Pns = FoldList[Function[{Pm, m},
        fp'[xp] D[Pm, xp] - (2 m - 1) fp''[xp] Pm], 1, Range[n - 1]] /.
    Derivative[m_][fp][xp] :> dfs[[m]];
  MapIndexed[{Pn, i} \[Function] Pn/dfs[[1]]^(2 i[[1]] - 1), Pns]
]

CompositeFunctionDerivativeList[G_, F_, X_, FSupp_:(0&)][m_, θ_] := Module[
  { ds, dF, df, fp },
  ds = InverseDerivativeList[m+1][X][θ];
  dF = EfficientDerivativeList[m][F][θ] + FSupp /@ Range[0, m];
  df = EfficientDerivativeList[m][G[fp]][θ] /.
   Map[Derivative[#][fp][θ] -> dF[[# + 1]] &, Range[0, m]];
  Table[Sum[df[[k+1]] BellY[j, k, ds[[;; j - k + 1]]], {k, 0, j}]/(j!), {j, 0, m}]
]

DScriptFPlusMinusDξθ0List[θ0_, θYL_, B_, C0_, CYL_, Gs_, gs_][m_] := CompositeFunctionDerivativeList[
    ScriptFPlusMinus, ReScriptFRegular[θ0, θYL, B, C0, CYL, Gs],
    ξ[θ0, gs], DReScriptFIrregular[θ0, B, C0]
  ][m, θ0]

DScriptFPlusMinusDξList[θ0_, θYL_, B_, C0_, CYL_, Gs_, gs_][m_, θ_] := CompositeFunctionDerivativeList[
    ScriptFPlusMinus, ReScriptF[θ0, θYL, B, C0, CYL, Gs], ξ[θ0, gs]
  ][m, θ]

DScriptF0DηList[θ0_, θYL_, B_, C0_, CYL_, Gs_, gs_][m_, θ_] := CompositeFunctionDerivativeList[
    ScriptF0[θ0, gs], ReScriptF[θ0, θYL, B, C0, CYL, Gs], η[θ0, gs]
  ][m, θ]

DScriptFPlusMinusDξθ0[params__][m_] := Last[DScriptFPlusMinusDξθ0List[params][m]]

DScriptFPlusMinusDξ[params__][m_, θ_] := Last[DScriptFPlusMinusDξList[params][m, θ]]

DScriptF0Dη[params__][m_, θ_] := Last[DScriptF0DηList[params][m, θ]]

DufDut[θ0_, θYL_, B_, C0_, CYL_, Gs_, gs_][m_][R_, θ_] := m! RealAbs[uh[θ0, gs][R, θ]]^(2 / Δ - m / Δ) DScriptF0Dη[θ0, θYL, B, C0, CYL, Gs, gs][m, θ] + Log[uh[θ0, gs][R, θ]^2] / (8 π Δ) Derivative[m][Function[utp, utp^2]][ut[R, θ]]

DufDuh[θ0_, θYL_, B_, C0_, CYL_, Gs_, gs_][m_][R_, θ_] := m! RealAbs[ut[R, θ]]^(2-m Δ) DScriptFPlusMinusDξ[θ0, θYL, B, C0, CYL, Gs, gs][m, θ] + Derivative[m][1&][θ] ut[R, θ]^2 / (8 π) Log[ut[R, θ]^2]

ruleB[θ0_, gs_] := (2 * OverlineS / π) * (- g[θ0, gs]'[θ0] / t[θ0]^Δ)

ruleC0[θ0_, gs_] := Exp[Δ t[θ0]^(Δ - 1) t'[θ0] / (2 OverlineS / π g[θ0, gs]'[θ0]) - t[θ0]^Δ g[θ0, gs]''[θ0] / (4 OverlineS / π g[θ0, gs]'[θ0]^2)] t[θ0]^(1/8) OverlineS / (2 π) * g[θ0, gs]'[θ0]

θ0Cas := Sqrt[1.16951];

h0Cas := a b ρ /. {
    a -> c2^2/c4,
    b -> (-c4/c2^3)^(1/2),
    ρ -> 2.00881
  } /. {
    c2->Ghs[[3]] 2!,
    c4->Ghs[[5]]["Value"] 4!
  }

gsCas := h0Cas * {
    1,
    -0.222389,
    -0.043547,
    -0.014809,
    -0.007168
  }

m0Cas := -Ghs[[3]]2! h0Cas

DScriptMCasDξList[m_, θ_] := CompositeFunctionDerivativeList[
    Identity, Function[θp, m0Cas * θp / RealAbs[θp^2 - 1]^β], ξ[θ0Cas, gsCas]
  ][m, θ]

End[]

EndPackage[]

