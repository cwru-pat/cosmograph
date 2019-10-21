SetDirectory[NotebookDirectory[]];
dirs = DirectoryName /@ FileNames["*/log.txt", "", \[Infinity]];
dirs

statdata = Import[# <> "calculated.dat.gz"] & /@ dirs;
frwdata = Import[# <> "calculated.frwdat.gz"] & /@ dirs;
r = 1;

NaturalSort = #[[
    Ordering@
     PadRight@
      StringSplit[#, x : DigitCharacter .. :> FromDigits@x]]] &;
ImportAggregate[f_] := 
  Block[{file = f, raydata, dirraydata, dirraydatasets},
   raydata = {};
   Do[
    If[
     Length[FileNames[dir <> file <> ".agg_values.h5.gz"]] == 1,
     
     Print["Importing aggregate data for " <> file <> " in " <> dir];
     dirraydata = 
      Import[dir <> file <> ".agg_values.h5.gz", {"Datasets", 
        "Dataset1"}];
     Print["...Done."],
     
     Print[
      "Importing non-aggregate data (and exporting for future) for " <>
        file <> " in " <> dir];
     dirraydatasets = 
      NaturalSort@Import[dir <> file <> ".values.h5.gz"];
     dirraydata = 
      Import[dir <> file <> ".values.h5.gz", {"Datasets", 
        dirraydatasets}];
     Export[dir <> file <> ".agg_values.h5.gz", dirraydata];
     Print["...Done."]
     ];
    AppendTo[raydata, dirraydata];,
    {dir, dirs}
    ];
   Return[raydata];
   ];
rayavgs = Import[# <> "raytracedata_avg.dat.gz"] & /@ dirs;
raydata = ImportAggregate["raytracedata"];
Dimensions /@ raydata
\[CapitalPhi]data = ImportAggregate["bardeen_phi_data"];
Dimensions /@ \[CapitalPhi]data
\[CapitalPsi]data = ImportAggregate["bardeen_psi_data"];
Dimensions /@ \[CapitalPsi]data
dtBdata = ImportAggregate["bardeen_dt_B_data"];
Dimensions /@ dtBdata

\[Psi]\[CapitalPhi]s = 
 Table[Integrate[
   Interpolation[{raydata[[n, ;; , r, 1]] - 1, \[CapitalPhi]data[[
        n, ;; , r, 1]]}\[Transpose]][z], {z, 0, 0.25}], {n, 4}]
\[Psi]\[CapitalPsi]s = 
 Table[Integrate[
   Interpolation[{raydata[[n, ;; , r, 1]] - 1, \[CapitalPsi]data[[
        n, ;; , r, 1]]}\[Transpose]][z], {z, 0, 0.25}], {n, 4}]

resolutions = {128, 160, 192, 256};
pairs = Sort /@ Subsets[Range[4], {2}]
triads = Sort /@ Subsets[Range[4], {3}]

{(\[Psi]\[CapitalPhi]s[[#[[1]]]] - \[Psi]\[CapitalPhi]s[[#[[
      2]]]])/(\[Psi]\[CapitalPhi]s[[#[[
      2]]]] - \[Psi]\[CapitalPhi]s[[#[[3]]]]), 
   GetConvergenceRate[resolutions[[#[[1]]]], resolutions[[#[[2]]]], 
    resolutions[[#[[3]]]], 1.]} & /@ triads
{(\[Psi]\[CapitalPsi]s[[#[[1]]]] - \[Psi]\[CapitalPsi]s[[#[[
      2]]]])/(\[Psi]\[CapitalPsi]s[[#[[
      2]]]] - \[Psi]\[CapitalPsi]s[[#[[3]]]]), 
   GetConvergenceRate[resolutions[[#[[1]]]], resolutions[[#[[2]]]], 
    resolutions[[#[[3]]]], 1.]} & /@ triads

RECoeffs = 
 Solve[{A/#[[1]] + B/#[[2]] == 0, A + B == 1}, {A, B}][[
    1]] & /@ (resolutions[[#]] & /@ pairs)

Table[A \[Psi]\[CapitalPhi]s[[pairs[[n, 1]]]] + 
   B \[Psi]\[CapitalPhi]s[[pairs[[n, 2]]]] /. RECoeffs[[n]], {n, 
  Length[pairs]}]
Table[A \[Psi]\[CapitalPsi]s[[pairs[[n, 1]]]] + 
   B \[Psi]\[CapitalPsi]s[[pairs[[n, 2]]]] /. RECoeffs[[n]], {n, 
  Length[pairs]}]

As = A /. RECoeffs;
Bs = B /. RECoeffs;

REData = Table[
   As[[n]] raydata[[pairs[[n, 1]], ;; ;; -As[[n]]]] + 
    Bs[[n]] raydata[[pairs[[n, 2]], ;; ;; Bs[[n]]]], {n, {1, 2, 4, 
     6}}];
Dimensions /@ REData

DAfrw[e_] := 3 /(2/(2/3 + 400 0.00078125))*2 1/e (1 - 1/Sqrt[e]);

\[Kappa]s = {};

Timing[Do[
  midpt = (Length[frwdata[[n]]] - 1)/2 + 1;
  
  aFRWs = Exp[2 frwdata[[n, midpt ;;, 1]]]/
   Exp[2 frwdata[[n, midpt, 1]]];
  zFRWs = 1/aFRWs - 1;
  \[Rho]FRWs = frwdata[[n, midpt ;;, 3]];
  KFRWs = frwdata[[n, midpt ;;, 2]];
  
  Hinvs = 1/(-frwdata[[n, midpt ;;, 2]]/3);
  
  HinvofZ = Interpolation[{zFRWs, Hinvs}\[Transpose]];
  RInt[zf_] := NIntegrate[HinvofZ[z], {z, 0, zf}];
  RFRWs = Table[{z, RInt[z]}, {z, zFRWs}];
  R = Interpolation[RFRWs];
  Rsrc = R[.25];
  
  AppendTo[\[Kappa]s, Table[
    \[CapitalPhi]Ray = \[CapitalPhi]data[[n, ;; ;; 10, rayn, 1]];
    \[Rho]Ray = raydata[[n, ;; ;; 10, rayn, 9]];
    \[Epsilon]0Ray = aFRWs^2/2 dtBdata[[n, ;; ;; 10, rayn, 1]];
    zRay = raydata[[n, ;; ;; 10, rayn, 1]] - 1;
    Rs = R /@ zRay;
    
    dr2\[CapitalPhi] = 
     D[Interpolation[{Rs, \[CapitalPhi]Ray}\[Transpose]][x], {x, 2}];
    dr2\[CapitalPhi]s = (dr2\[CapitalPhi] /. x -> #) & /@ Rs;
    
    \[Delta]\[Rho] = (\[Rho]Ray - \[Rho]FRWs) + \[Epsilon]0Ray KFRWs \
\[Rho]FRWs;
    \[Delta]u = -\[Epsilon]0Ray;
    
    Uobs = Interpolation[{Rs, \[Delta]u}\[Transpose]]'[0];
    Usrc = Interpolation[{Rs, \[Delta]u}\[Transpose]]'[Rsrc];
    
    Ifn = 
     Interpolation[{Rs, (Rsrc - Rs) Rs/
         Rsrc (4 \[Pi] aFRWs^2 \[Delta]\[Rho] + 
           4 \[Pi] aFRWs^2 KFRWs \[Rho]FRWs \[Delta]u - 
           dr2\[CapitalPhi]s)}\[Transpose]];
    \[Kappa]GR = 
     Interpolation[{Rs, (
        DAfrw[raydata[[n, ;; ;; 10, rayn, 1]]] - 
         raydata[[n, ;; ;; 10, rayn, 6]])/(
        DAfrw[raydata[[n, ;; ;; 10, rayn, 1]]] + 10^-10)}\[Transpose]];
    {NIntegrate[Ifn[r], {r, 0, Rsrc}], \[Kappa]GR[Rsrc], Uobs, Usrc},
    {rayn, 12 16 16}
    ]];
  
  , {n, 4} 
  ]]

RE\[Kappa]s = 
  Table[As[[n]] \[Kappa]s[[pairs[[n, 1]]]] + 
    Bs[[n]] \[Kappa]s[[pairs[[n, 2]]]], {n, {1, 2, 3, 4, 5, 6}}];

Dimensions /@ RE\[Kappa]s

Hobs = -1/3 (-2.04255319148937`);
hofz = Interpolation[{Exp[2 frwdata[[4, 161, 1]]]/
       Exp[2 frwdata[[4, 161 ;;, 1]]] - 1, 
     1/(-1/3 frwdata[[4, 161 ;;, 2]])}\[Transpose]];
\[Chi] = Integrate[hofz[z], {z, 0, zsrc}];

Needs["ErrorBarPlots`"]
REerrpltdata = 
  Table[{{Mean[#[[;; , 1]] - 
         1./(Hobs \[Chi]) (#[[;; , 4]] - #[[;; , 3]]) + #[[;; , 3]]], 
       Mean[#[[;; , 2]]]}, 
      ErrorBar[StandardDeviation[#[[;; , 1]]], 
       StandardDeviation[#[[;; , 2]]]]} &@RE\[Kappa]s[[;; , r]], {r, 
    Length[RE\[Kappa]s[[1]]]}];

Export["grf_Nside16_kappx.txt", 
 REerrpltdata[[;; 12 16^2, 1, 1]], "Table"]
Export["grf_Nside16_kGR.txt", 
 REerrpltdata[[;; 12 16^2, 1, 2]], "Table"]
Export["grf_Nside16_kGRminuskappx.txt", 
 REerrpltdata[[;; 12 16^2, 1, 2]] - 
  REerrpltdata[[;; 12 16^2, 1, 1]], "Table"]
