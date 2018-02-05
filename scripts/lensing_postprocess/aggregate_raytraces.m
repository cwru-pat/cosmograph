
dirs = DirectoryName /@ FileNames["*/log.txt", "", \[Infinity]];
dirs

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

ImportAggregate["raytracedata"];
ImportAggregate["bardeen_phi_data"];
ImportAggregate["bardeen_psi_data"];
ImportAggregate["bardeen_dt_B_data"];
