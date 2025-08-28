




(*****************************)
(*     Write Feynman Rules   *)
(*****************************)

(* Substitutions for writing Feynman Rules in FORM language  *)



SubstitutionsLorentzFRtoForm={



  (* Product of Gamma Matrices *)

  (* L$FV(3, Index(Lorentz, Int(1)))*L$Ga(5, Index(Spin, Int(2)), Index(Spin, Ext(2)))*L$Ga(Index(Lorentz, Ext(3)), Index(Spin, Int(3)), Index(Spin, Int(2)))*L$Ga(Index(Lorentz, Int(1)), Index(Spin, Ext(1)), Index(Spin, Int(3)))) *)


  L$FV[ w_?(NumericQ[#]&), Index[Lorentz, Int[ a_ ] ] ]*L$Ga[ Index[Lorentz, Ext[ z_ ] ], Index[Spin, Int[ b_ ] ], Index[Spin, Int[c_ ] ] ]*L$Ga[Index[Lorentz, Int[a_ ] ], Index[Spin, Ext[x_ ] ], Index[Spin, Int[b_ ] ] ]*L$ProjP[Index[Spin, Int[c_ ] ], Index[Spin, Ext[y_ ] ] ]  :> 1/2 (FO$Gamma[MomentumFO[[w]],LorentzIndexFO[[z]],SpinIndexFO[[x]],SpinIndexFO[[y]]]+FO$Gamma[MomentumFO[[w]],LorentzIndexFO[[z]],5,SpinIndexFO[[x]],SpinIndexFO[[y]]]),

  L$FV[w_?(NumericQ[#]&), Index[Lorentz, Int[a_]]]*L$Ga[Index[Lorentz, Ext[z_]], Index[Spin, Int[b_]], Index[Spin, Int[c_]]]*L$Ga[Index[Lorentz, Int[a_]], Index[Spin, Ext[x_]], Index[Spin, Int[b_]]]*L$ProjM[Index[Spin, Int[c_ ] ], Index[Spin, Ext[y_] ] ] :> 1/2 (FO$Gamma[MomentumFO[[w]],LorentzIndexFO[[z]],SpinIndexFO[[x]],SpinIndexFO[[y]]]-FO$Gamma[MomentumFO[[w]],LorentzIndexFO[[z]],5,SpinIndexFO[[x]],SpinIndexFO[[y]]]),

  L$FV[ w_?(NumericQ[#]&), Index[Lorentz, Int[ a_ ] ] ]*L$Ga[ Index[Lorentz, Ext[ z_ ] ], Index[Spin, Ext[x_ ] ], Index[Spin, Int[b_ ] ] ]*L$Ga[Index[Lorentz, Int[a_ ] ], Index[Spin, Int[ b_ ] ], Index[Spin, Int[c_ ] ]  ]*L$ProjP[Index[Spin, Int[c_ ] ], Index[Spin, Ext[y_ ] ] ]  :> 1/2 (FO$Gamma[LorentzIndexFO[[z]],MomentumFO[[w]],SpinIndexFO[[x]],SpinIndexFO[[y]]]+FO$Gamma[LorentzIndexFO[[z]],MomentumFO[[w]],5,SpinIndexFO[[x]],SpinIndexFO[[y]]]),

  L$FV[ w_?(NumericQ[#]&), Index[Lorentz, Int[ a_ ] ] ]*L$Ga[ Index[Lorentz, Ext[ z_ ] ], Index[Spin, Ext[x_ ] ], Index[Spin, Int[b_ ] ] ]*L$Ga[Index[Lorentz, Int[a_ ] ], Index[Spin, Int[ b_ ] ], Index[Spin, Int[c_ ] ]  ]*L$ProjM[Index[Spin, Int[c_ ] ], Index[Spin, Ext[y_ ] ] ]  :> 1/2 (FO$Gamma[LorentzIndexFO[[z]],MomentumFO[[w]],SpinIndexFO[[x]],SpinIndexFO[[y]]]-FO$Gamma[LorentzIndexFO[[z]],MomentumFO[[w]],5,SpinIndexFO[[x]],SpinIndexFO[[y]]]),

  L$FV[ w_?(NumericQ[#]&), Index[Lorentz, Int[a_] ] ]*L$Ga[5, Index[Spin, Int[c_] ], Index[Spin, Ext[y_] ] ]*L$Ga[Index[Lorentz, Ext[z_] ], Index[Spin, Int[b_] ], Index[Spin, Int[c_] ] ]*L$Ga[Index[Lorentz, Int[a_] ], Index[Spin, Ext[x_] ], Index[Spin, Int[b_] ] ] :> FO$Gamma[MomentumFO[[w]],LorentzIndexFO[[z]],5,SpinIndexFO[[x]],SpinIndexFO[[y]]],

  L$FV[ w_?(NumericQ[#]&), Index[Lorentz, Int[a_] ] ]*L$Ga[5, Index[Spin, Int[c_] ], Index[Spin, Ext[y_] ] ]*L$Ga[Index[Lorentz, Ext[z_] ], Index[Spin, Ext[x_] ], Index[Spin, Int[b_] ] ]*L$Ga[Index[Lorentz, Int[a_] ],  Index[Spin, Int[b_] ], Index[Spin, Int[c_] ] ] :> FO$Gamma[LorentzIndexFO[[z]],MomentumFO[[w]],5,SpinIndexFO[[x]],SpinIndexFO[[y]]],

  L$FV[ w_?(NumericQ[#]&), Index[Lorentz, Int[ a_ ] ] ]*L$Ga[ Index[Lorentz, Ext[ z_ ] ], Index[Spin, Ext[x_ ] ], Index[Spin, Int[b_ ] ] ]*L$Ga[Index[Lorentz, Int[a_ ] ], Index[Spin, Int[ b_ ] ], Index[Spin,  Ext[y_ ] ]  ]  :>  FO$Gamma[LorentzIndexFO[[z]],MomentumFO[[w]],SpinIndexFO[[x]],SpinIndexFO[[y]] ],

  L$FV[ w_?(NumericQ[#]&), Index[Lorentz, Int[ a_ ] ] ]*L$Ga[ Index[Lorentz, Ext[ z_ ] ],  Index[Spin, Int[ b_ ] ], Index[Spin,  Ext[y_ ] ] ]*L$Ga[Index[Lorentz, Int[a_ ] ],  Index[Spin, Ext[x_ ] ], Index[Spin, Int[b_ ] ] ]  :>  FO$Gamma[MomentumFO[[w]],LorentzIndexFO[[z]],SpinIndexFO[[x]],SpinIndexFO[[y]] ],

  L$Ga[Index[Lorentz, Ext[z_] ], Index[Spin, Ext[x_] ], Index[Spin, Int[a_] ] ]*L$Ga[Index[Lorentz, Ext[w_] ], Index[Spin, Int[a_] ], Index[Spin, Int[b_] ] ]*L$ProjM[Index[Spin, Int[b_] ], Index[Spin, Ext[y_] ] ] :> 1/2 ( FO$Gamma[LorentzIndexFO[[z]],LorentzIndexFO[[w]],SpinIndexFO[[x]],SpinIndexFO[[y]] ] - FO$Gamma[LorentzIndexFO[[z]],LorentzIndexFO[[w]],5,SpinIndexFO[[x]],SpinIndexFO[[y]] ] ),

  L$Ga[Index[Lorentz, Ext[z_] ], Index[Spin, Ext[x_] ], Index[Spin, Int[a_] ] ]*L$Ga[Index[Lorentz, Ext[w_] ], Index[Spin, Int[a_] ], Index[Spin, Int[b_] ] ]*L$ProjP[Index[Spin, Int[b_] ], Index[Spin, Ext[y_] ] ] :> 1/2 ( FO$Gamma[LorentzIndexFO[[z]],LorentzIndexFO[[w]],SpinIndexFO[[x]],SpinIndexFO[[y]] ] + FO$Gamma[LorentzIndexFO[[z]],LorentzIndexFO[[w]],5,SpinIndexFO[[x]],SpinIndexFO[[y]] ] ),

  L$Ga[Index[Lorentz, Int[z_]], Index[Spin, Ext[x_]], Index[Spin, Int[ii_]]]*L$ProjP[Index[Spin, Int[ii_]], Index[Spin, Ext[y_]]]:> 1/2 (FO$Gamma[iLorentzIndexFO[[z]],SpinIndexFO[[x]],SpinIndexFO[[y]]]+FO$Gamma[iLorentzIndexFO[[z]],5,SpinIndexFO[[x]],SpinIndexFO[[y]]]),

  L$Ga[Index[Lorentz, Int[z_]], Index[Spin, Ext[x_]], Index[Spin, Int[ii_]]]*L$ProjM[Index[Spin, Int[ii_]], Index[Spin, Ext[y_]]]:> 1/2 (FO$Gamma[iLorentzIndexFO[[z]],SpinIndexFO[[x]],SpinIndexFO[[y]]]-FO$Gamma[iLorentzIndexFO[[z]],5,SpinIndexFO[[x]],SpinIndexFO[[y]]]),
  
  L$Ga[Index[Lorentz, Ext[z_]], Index[Spin, Ext[x_]], Index[Spin, Int[ii_]]]*
   L$ProjP[Index[Spin, Int[ii_]], Index[Spin, Ext[y_]]]:> 1/2 (FO$Gamma[LorentzIndexFO[[z]],SpinIndexFO[[x]],SpinIndexFO[[y]]]+FO$Gamma[LorentzIndexFO[[z]],5,SpinIndexFO[[x]],SpinIndexFO[[y]]]),

  L$Ga[Index[Lorentz, Ext[z_]], Index[Spin, Ext[x_]], Index[Spin, Int[ii_]]]*
   L$ProjM[Index[Spin, Int[ii_]], Index[Spin, Ext[y_]]]:> 1/2 (FO$Gamma[LorentzIndexFO[[z]],SpinIndexFO[[x]],SpinIndexFO[[y]]]-FO$Gamma[LorentzIndexFO[[z]],5,SpinIndexFO[[x]],SpinIndexFO[[y]]]),



  (* Gamma Matrix  *)
  L$Ga[Index[Lorentz,Ext[z_]],Index[Spin,Ext[x_]],Index[Spin,Ext[y_]]]:>FO$Gamma[LorentzIndexFO[[z]],SpinIndexFO[[x]],SpinIndexFO[[y]]],
  L$Ga[Index[Lorentz,Ext[z_]],Index[Spin,Ext[x_]],Index[Spin,Int[y_]]]:>FO$Gamma[LorentzIndexFO[[z]],SpinIndexFO[[x]],iSpinIndexFO[[y]]],

  (* Delta *)

  L$IndexDelta[Index[Spin,Ext[x_]],Index[Spin,Ext[y_]]]:>FO$Delta[SpinIndexFO[[x]],SpinIndexFO[[y]]],

  (* Momentum *)
  L$FV[x_?(NumericQ[#]&),Index[Lorentz,Ext[y_]]]:>  MomentumFO[[x]][LorentzIndexFO[[y]]] ,

  L$FV[ x_?(NumericQ[#]&), Index[Lorentz, Int[z_] ] ]*L$FV[ y_?(NumericQ[#]&), Index[Lorentz, Int[z_] ] ] :> FO$Dot[ MomentumFO[[x]], MomentumFO[[y]] ],

  (* Gamma 5 *)

  L$Ga[5,Index[Spin,Ext[x_]],Index[Spin,Ext[y_]]]:>FO$Gamma[5,SpinIndexFO[[x]],SpinIndexFO[[y]]],

  (* Metric tensor *)
  L$ME[Index[Lorentz,Ext[x_]],Index[Lorentz,Ext[y_]]]:> FO$Metric[LorentzIndexFO[[x]],LorentzIndexFO[[y]]],

  (* Projectors *)
  L$ProjP[Index[Spin,Ext[x_]],Index[Spin,Ext[y_]]]:> 1/2 (  FO$Gamma[1,SpinIndexFO[[x]],SpinIndexFO[[y]]] + FO$Gamma[5,SpinIndexFO[[x]],SpinIndexFO[[y]]]),
  L$ProjM[Index[Spin,Ext[x_]],Index[Spin,Ext[y_]]]:> 1/2 (  FO$Gamma[1,SpinIndexFO[[x]],SpinIndexFO[[y]]] - FO$Gamma[5,SpinIndexFO[[x]],SpinIndexFO[[y]]]),
  (* L$ProjP[Index[Spin,Int[x_]],Index[Spin,Ext[y_]]]:> 1/2 ( 1+FO$Gamma[5,iSpinIndexFO[[x]],SpinIndexFO[[y]]]),*)
  (* L$ProjM[Index[Spin,Int[x_]],Index[Spin,Ext[y_]]]:> 1/2 ( 1-FO$Gamma[5,iSpinIndexFO[[x]],SpinIndexFO[[y]]])*)

  YY_[Index[Generation, x_], Index[Generation, y_]]:> YY[x,y]

};


FO$GaugeIndex[GG_,x_]:=If[x>0,IndicesType[GG][[x]],IntIndicesType[GG][[x]]]

SubstitutionsColorFRtoForm={

(* Structure constant *)
  f[x_,y_,z_]:> fabc[FO$GaugeIndex[Gluon,x],FO$GaugeIndex[Gluon,y],FO$GaugeIndex[Gluon,z] ],
 
(* Generator *)
  T[x_,y_,z_]:> TM[FO$GaugeIndex[Colour,y],FO$GaugeIndex[Colour,z],FO$GaugeIndex[Gluon,x] ],

(* Delta *)

  FO$Delta[x_,y_] :> FO$Delta[FO$GaugeIndex[Colour,x],FO$GaugeIndex[Colour,y] ]

};



(* Writing of the Feynman rules in terms of GC couplings *)

Get[$FeynRulesPath<>"/Interfaces/UFO/PYIntVertices.m"];
Get[$FeynRulesPath<>"/Interfaces/UFO/PYIntParticles.m"];

(* LorentzStructures provides substitution rules for the YYYY lorentz names *)
LorentzStructures[iLorentzObjects_]:=iLorentzObjects/.{LorentzObject[name_,_,struct_]:>name-> struct }

(* FO$GetLorentzStruct takes the format "'YYYY'" where L.YYYY as input and gives the corresponding mathematical expression of the Lorentz Structure *)
FO$GetLorentzStruct[lorentzObjects_List,YYYY_String]:=Module[{Structures},

  Structures=(LorentzStructures[#]&)/@lorentzObjects;
  YYYY/.Structures

]


(* Write vertex in terms of GC couplings *)
(*FO$Vertex[list_List,lorentzStruct_]:=Module[{lorentzList,couplingList,colorList},

  lorentzList=StringCases[Cases[list,{"lorentz",str_String}:>str][[1]],"L."~~name:WordCharacter..:>name];
  couplingList=StringCases[Cases[list,{"couplings",str_String}:>str][[1]],"C.GC_"~~n:DigitCharacter..:>n];
  colorList=Cases[list,{"color",str_String}:>str][[1]];
  colorList=StringReplace[colorList,{"'"->"", "["->"", "]"->"", "("->"[", ")"->"]", "Identity" -> "FO$Delta"}];
  colorList=(ToExpression/@StringSplit[colorList,", "])/.SubstitutionsColorFRtoForm;

  If[Length[colorList]==1,
    colorList[[1]]*Plus@@MapThread[Function[{lorentz,coupling},FO$GetLorentzStruct[lorentzStruct,lorentz]*Symbol["GC"<>coupling]],{lorentzList,couplingList}],
    Plus@@MapThread[Function[{lorentz,coupling,color},color*FO$GetLorentzStruct[lorentzStruct,lorentz]*Symbol["GC"<>coupling]],{lorentzList,couplingList,colorList}]
    ]
]*)

FO$Vertex[list_List,lorentzStruct_]:=Module[{lorentzList,couplingList,colorList,Positions,terms},

lorentzList=StringCases[Cases[list,{"lorentz",str_String}:>str][[1]],"L."~~name:WordCharacter..:>name];
lorentzList=FO$GetLorentzStruct[lorentzStruct,#]&/@lorentzList;

couplingList=StringSplit[Cases[list,{"couplings",str_String}:>str][[1]],",("];
couplingList=Map[(Symbol["GC"<>#]&),StringCases[#,"C.GC_"~~n:DigitCharacter..:>n]&/@couplingList,{2}];
couplingList=Total[#]&/@couplingList;

colorList=Cases[list,{"color",str_String}:>str][[1]];
colorList=StringReplace[colorList,{"'"->"","["->"","]"->"","("->"[",")"->"]","Identity"->"FO$Delta"}];
colorList=(ToExpression/@StringSplit[colorList,", "])/.SubstitutionsColorFRtoForm;

Positions=1+ToExpression[StringCases[Cases[#,{"couplings",str_String}:>str][[1]],"("~~i:DigitCharacter..~~","~~j:DigitCharacter..~~")":>{i,j}]&@list];

terms={};
Do[ terms=Append[terms,colorList[[Positions[[i,1]]]]*lorentzList[[Positions[[i,2]]]]*couplingList[[i]]],{i,1,Length[Positions]}];
terms=Total[terms];
Return[terms];

];

FO$GCVertices[FR_]:=Block[{SplittedVertices},

  SplittedVertices=MapIndexed[CreateVertexObjectEntry,PYSplitVertices[FR] ];
  FO$LorentzObjects=PY$LorentzObjects//.SubstitutionsLorentzFRtoForm;
  FO$FRw=FO$Vertex[#,FO$LorentzObjects]&/@SplittedVertices
]


(* GCouplingsWithIm for finding GCouplings with an Imaginary *)

GCouplingsWithIm[list_List] := (Cases[list, HoldPattern[_ -> __*im]])[[;; , 1]];

(* GCouplingsEquivalencies for finding GCouplings that are equal *)

rulesGC[list_List]:=Module[{firstElement,restElements},

  firstElement = First[list];
  restElements = Rest[list];
  Thread[restElements->firstElement]

];

GCouplingsEquivalencies[couplingsRules_List]:=Module[{GroupedbyRHS,firstGCi},

  GroupedbyRHS = GatherBy[couplingsRules,Last];
  GroupedbyRHS = Select[GroupedbyRHS,Length[#]>1&];

  GroupedbyRHS = rulesGC[#]&/@(GroupedbyRHS[[;;,;;,1]]);
  GroupedbyRHS = Flatten[GroupedbyRHS]

];


(* Function to add the repeat lines to vertices with dummy indices *)

IdentifyiGluons[line_String] := DeleteDuplicates[StringCases[line, "iGluon" ~~ DigitCharacter ..] ];

(* modifyVerticesFile[file_]:=Module[{content,modifiedContent,lines,stringQ,modifiedLines,newFile},

  content=Import[file,"Text"];

  (*Identify lines containing "iGluon999"*)

  lines=StringSplit[content,"\n"];
  stringQ=StringContainsQ[lines,"iGluon999"];
  modifiedLines=Position[stringQ,True];

  (*Modify the lines*)

  Do[
    modifiedContent=StringInsert[lines[[modifiedLines[[ii]]]],"repeat;\n",1];
    modifiedContent=StringReplace[modifiedContent,"id "-> "id, once "];
    lines[[modifiedLines[[ii]]]]=StringInsert[modifiedContent,"\nsum iGluon999;\nendrepeat;\n.sort",-1],

  {ii,1,Length@modifiedLines}];

  newFile=StringRiffle[lines,"\n"];
  newFile=StringInsert[newFile,"\n",-1];


  Export[file,newFile,"Text"]
] *)

(* modifyVertices adds the repeat and id, once statements whenever there is an iGluon index *)

modifyVertices[file_]:=Module[{content,modifiedContent,lines,stringQ,modifiedLines,newFile,InternalGluonIndices,currentIndex},

  content=Import[file,"Text"];

  (*Identify lines containing iGluon*)

  lines=StringSplit[content,"\n"];

  (* Get the iGluon indices *)
  InternalGluonIndices=IdentifyiGluons[#]&/@lines;

  (* Find the lines to be modified *)
  modifiedLines=Flatten@Position[InternalGluonIndices,Except[{}],{1}];
  modifiedLines=modifiedLines[[2;;-1]];

  (* Clean the list of iGluon indices*)
  InternalGluonIndices=Flatten[InternalGluonIndices];
  InternalGluonIndices=DeleteDuplicates[InternalGluonIndices];

  Do[

    modifiedContent=StringInsert[lines[[modifiedLines[[ii]]]],"repeat;\n",1];
    modifiedContent=StringReplace[modifiedContent,"id "->"id, once "];
    lines[[modifiedLines[[ii]]]]=StringInsert[modifiedContent,"\nendrepeat;\n.sort",-1],

  {ii,1,Length@modifiedLines}];

  newFile=StringRiffle[lines,"\n"];
  newFile=StringInsert[newFile,"\n",-1];

  Export[file,newFile,"Text"];

];

(* modifyVerticesGluon adds sum iGluon whenever needed *)

(* add a substring to iGluons *)
ChangeGluons[str_String,add_String] := StringReplace[str, RegularExpression["iGluon(\\d+)"] :> ("iGluon" <> "$1" <> add)];

modifyVerticesGluon[file_]:=Module[{content,modifiedContent,lines,stringQ,modifiedLines,newFile,InternalGluonIndices,currentIndex},


  content=Import[file,"Text"];

  (*Identify lines containing iGluon*)

  lines=StringSplit[content,"\n"];

  (* Get the iGluon indices *)
  InternalGluonIndices=IdentifyiGluons[#]&/@lines;

  (* Find the lines to be modified *)
  modifiedLines=Flatten@Position[InternalGluonIndices,Except[{}],{1}];
  modifiedLines=modifiedLines[[2;;-1]];

  (* Clean the list of iGluon indices*)
  InternalGluonIndices=Flatten[InternalGluonIndices];
  InternalGluonIndices=DeleteDuplicates[InternalGluonIndices];

  Do[

    currentIndex=InternalGluonIndices[[jj]];

    Do[

      If[StringContainsQ[lines[[modifiedLines[[ii]]]],currentIndex],

      	  modifiedContent=lines[[modifiedLines[[ii]]]];
      	  lines[[modifiedLines[[ii]]]]=StringInsert[modifiedContent,"\nsum "<>ToString[currentIndex]<>";",-1]

        ];

      (* modifiedContent=StringInsert[lines[[modifiedLines[[ii]]]],"repeat;\n",1]; *)

      ,{ii,1,Length@modifiedLines}];

    ,{jj,1,Length[InternalGluonIndices]}];

  (*Modify the lines*)

  newFile=StringRiffle[lines,"\n"];
  newFile=StringInsert[newFile,"\n",-1];

  newFile=ChangeGluons[newFile,"00"];

  Export[file,newFile,"Text"]

];


(* Options of WriteVertices *)

Options[FO$WriteVertices]={InternalParameters-> {}};




(*  Write Feynman Rules *)

FO$WriteVertices[out0_,FeynmanRules0_,opts:options___]:=Block[{out=out0,FeynmanRules=FeynmanRules0,outfile,VRTX,verticesExpressions,InputParameterList,Qi,ParameterValue,Subst,Fermions,CouplingswIm,newCouplings,CouplingList,stringSubstitutions,fieldVertexSubs,exprsVertexSubs,ast,lineLength},

  (* Routine for creating the list of substitutions of the desired parameters *)

  InputParameterList=InternalParameters/.{options}/.Options[FO$WriteVertices];
  SubstitutionParameters={};
  Do[
    Qi=InputParameterList[[j]];
    ParameterValue=Value/.MR$ParameterRules[Qi];

    If[Head[(Value/.MR$ParameterRules[#]&)@Qi]=!=  List, Subst=Replace[Qi, {Qi-> {Qi-> Value/.MR$ParameterRules[Qi]}}], Subst=Replace[Qi,Qi-> Value/.MR$ParameterRules[Qi]]];

    SubstitutionParameters=AppendTo[SubstitutionParameters,Subst]
  ,{j,1,Length[InputParameterList]}];

  SubstitutionParameters=Flatten[SubstitutionParameters];



  (* Identify the fermions in the Feynman Rules  -  This is redundant as it is also computed in the Definitions routine  *)

  Fermions=Select[allFieldsInFR,FermionQ[#]&&Not[GhostFieldQ[#]===True]&];



  (* Organization of the vertices functions  *)

  VRTX=vrtx@@@(FeynmanRules0[[;;,1]]/.{
  {Field_,index_?(NumericQ[#]&)}:> Field [TagIndex[[index]],MomentumFO[[index]],LorentzIndexFO[[index]],IndicesSequence[Field,index]]
  });

  (* Substituted expressions of the Feynman Rules  *)


  (* verticesExpressions=FeynmanRules0[[;;,2]]//.SubstitutionsFRtoForm//.SubstitutionParameters//.{Complex[a_,b_]:>a+ im b}//Simplify; *)

(*  verticesExpressions = FO$GCVertices[FeynmanRules]//.SubstitutionsFRtoForm//.{Complex[a_,b_]:>a+ im b};*)
  verticesExpressions = FO$GCVertices[FeynmanRules]//.{Complex[a_,b_]:>a+ im b};

  (* Write the couplings file *)

   CouplingList = FO$WriteCouplings[out,opts];

  (* Factor out the imaginary from GCouplings *)

   CouplingList = StringReplace[CouplingList,{"=":>"->"}];
   CouplingList = ToExpression/@CouplingList;
   CouplingList = Simplify/@CouplingList; (* Ensures the correct identification of couplings with im *)

   CouplingswIm = GCouplingsWithIm[CouplingList];
   CouplingswIm = Replace[#, {x_ -> (x -> Im*x)}] & /@ CouplingswIm;

   verticesExpressions = verticesExpressions/.CouplingswIm;

  (* Re-write the couplings file *)

   newCouplings = Import[out<>"/Couplings.frm", "Text"];
   newCouplings = StringReplace[newCouplings, "im" -> "1"];

   Export[out<>"/Couplings.frm", newCouplings, "Text"];

  (* GCouplings that are equal for simplification of Feynman Rules *)
   
   equalGCouplings = GCouplingsEquivalencies[CouplingList];

   verticesExpressions = Simplify[verticesExpressions/.equalGCouplings];

  (* stringSubstitutions = {

    "iGluon10" :> "iGluon999"

  }; *)

  fieldVertexSubs = {
    ("PartTag"~~char:LetterCharacter):>"?"<>char,
    ("p"~~char:DigitCharacter):>"p"<>char<>"?",
    ("Lor"~~char:DigitCharacter):>"Lor"<>char<>"?",
    ("Spin"~~char:DigitCharacter):>"Spin"<>char<>"?",
    ("Gluon"~~char:DigitCharacter):>"Gluon"<>char<>"?",
    ("Colour"~~char:DigitCharacter):>"Colour"<>char<>"?",
    "["-> "(",
    "]"-> ")"
  };

  exprsVertexSubs={
      "Sqrt"-> "sqrt_",
      "Im"->"im",
      "FO$Metric"-> "d_",
      "FO$Delta"-> "d_",
      ("QuestionMark"~~char:LetterCharacter):>"?"<>char,
      "["-> "(",
      "]"-> ")",
      "FO$Gamma"-> "GammaM",
      "FO$Dot[" ~~Shortest[ p1___ ] ~~ ", " ~~ Shortest[p2___] ~~ "]" :> p1 <> "." <> p2
      (* "iGluon10" :> "iGluon999" *)
};


  outfile=out<>"/Vertices.frm";

  lineLength=100;
  wsp="                                                                                                                             ";   
  ast="*****************************************************************************************************************************";   
  info=(ToString[#1]<>" : "<>If[Head[#2]===List,StringDrop[StringJoin[(#<>", "&)/@#2],-2],#2]&)@@@M$Information;
  Print["Writing  vertices model in "<>outfile];
  OpenWrite[outfile];
    WriteString[outfile,StringTake[ast,lineLength+1]<>"\n"];
    WriteString[outfile,StringTake["* File automalically generated by FeynRules -"<>wsp,lineLength]<>"*\n"];
    WriteString[outfile,StringTake["* date and time of generation : "<>DateString[]<>wsp,lineLength]<>"*\n"];
    WriteString[outfile,StringTake["* FeynRules model information : "<>wsp,lineLength]<>"*\n"];
    WriteString[outfile, HeadFormat[M$Information,lineLength,"*","*"]<>"\n"];
    WriteString[outfile,StringTake[ast,lineLength+1]<>"\n"];
    WriteString[outfile,"\n"];
    WriteString[outfile,"\n"];
    WriteString[outfile,"*****  This file contains the expressions of the vertices according to the syntax from QGraf  *****\n"];
    WriteString[outfile,"\n"];
    WriteString[outfile,"\n"];
    WriteString[outfile,"#procedure FeynmanRules\n"];
    WriteString[outfile,".sort\n"];
    WriteString[outfile,"\n"];
    WriteString[outfile,"argument;\n"];
    WriteString[outfile,"id   ff?Field(x?partTagExt[n],?a) = ff(x,?a,ExtLor[n]);\n"];
    WriteString[outfile,"id   ff?Field(x?partTagInt[n],?a) = ff(x,?a,IntLor[n]);\n"];
    WriteString[outfile,"id   ff?Fermion(x?partTagExt[n],?a) = ff(x,?a,Spin[n]);\n"];
    WriteString[outfile,"id   ff?Fermion(x?partTagInt[n],?a) = ff(x,?a,IntSpin[n]);\n"];
    Do[
      WriteString[outfile,"id   ff?"<>ToString[FO$FieldsClassified[[j,1,1]]]<>"ed(x?partTagExt[n],?a) = ff(x,?a,"<>ToString[FO$FieldsClassified[[j,1,1]]]<>"[n]);\n"]
      WriteString[outfile,"id   ff?"<>ToString[FO$FieldsClassified[[j,1,1]]]<>"ed(x?partTagInt[n],?a) = ff(x,?a,Int"<>ToString[FO$FieldsClassified[[j,1,1]]]<>"[n]);\n"]
    ,{j,1,Length[FO$IndicesInFR]}];
    WriteString[outfile,"endargument;\n"];
    WriteString[outfile,"\n"];
    WriteString[outfile,"\n"];
    comma=False;
    Do[
      If[comma,WriteString[outfile,";\n\n"] ];

  	  comma=False;
      (* WriteString[outfile,"\nid \t"<>((StringReplace[#,("PartTag"~~char:LetterCharacter):>"?"<>char]&)@(StringReplace[#,("p"~~char:DigitCharacter):>"p"<>char<>"?"]&)@(StringReplace[#,("Lor"~~char:DigitCharacter):>"Lor"<>char<>"?"]&)@(StringReplace[#,("Spin"~~char:DigitCharacter):>"Spin"<>char<>"?"]&)@(StringReplace[#,("Gluon"~~char:DigitCharacter):>"Gluon"<>char<>"?"]&)@(StringReplace[#,("Colour"~~char:DigitCharacter):>"Colour"<>char<>"?"]&)@(StringReplace[#,"["-> "("]&)@(StringReplace[#,"]"-> ")"]&)@(ToString[InputForm[#] ])&)@VRTX[[j]]]; *)
      WriteString[outfile,"\nid \t"<>((StringReplace[#,fieldVertexSubs]&)@(ToString[InputForm[#] ])&)@VRTX[[j]]];

      (* WriteString[outfile,"\t = \t "<>((StringReplace[#,{"Sqrt"-> "sqrt_","Im"->"im"}]&)@(StringReplace[#,"FO$Metric"-> "d_"]&)@(StringReplace[#,"FO$Delta"-> "d_"]&)@(StringReplace[#,("QuestionMark"~~char:LetterCharacter):>"?"<>char]&)@(StringReplace[#,"["-> "("]&)@(StringReplace[#,"]"-> ")"]&)@(StringReplace[#,"FO$Gamma"-> "GammaM"]&)@(StringReplace[ #,"FO$Dot[" ~~Shortest[ p1___ ] ~~ ", " ~~ Shortest[p2___] ~~ "]" :> p1 <> "." <> p2]&)@(StringReplace[#,stringSubstitutions]&)@(ToString[InputForm[#] ])&)@verticesExpressions[[j]] ]; *)
      WriteString[outfile,"\t = \t "<>((StringReplace[#,exprsVertexSubs]&)@(ToString[InputForm[#] ])&)@verticesExpressions[[j]] ];
      WriteString[outfile,";\n"];
    ,{j,1,Length[verticesExpressions]}];
    WriteString[outfile,"\n"];
    WriteString[outfile,"\n"];
    WriteString[outfile,"#endprocedure FeynmanRules \n\n"];
  Close[outfile];

  (* the routines below add "repeat endrepeat" and gluon indices sums *)

  modifyVertices[outfile];
  modifyVerticesGluon[outfile];
  ChangeGluons
]