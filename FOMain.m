(* ::Module:: *)

(* ::Text:: *)
(*Authors : A. Vasquez 2025*)
(*    *)

FFO$Version = "v1.0\[Beta]1";

Print["Version: "<>FFO$Version<>", (28 Aug 2025)"];
Print["Authors:  C. Duhr, P. Mukherjee, A. Vasquez"];
Print["Please cite  "];
Print["https://github.com/ANATAR-hep/FeynRules-Model-Interface"];

(* ::Section:: *)
(*Some printout*)


<< "FOVertices.m";
<< "FOPolarizations.m";
<< "FOPropagators.m";
<< "FOCouplings.m";


QGCouplings::usage = "QGCouplings -> { \[Alpha]1, \[Alpha]2,... } is an option of WriteFFO and chooses the relevant couplings in the QGraf model, which allows diagram filtering. "

InternalParameters::usage = "InternalParameters is an option of WriteFFO that specifies the internal parameters that should be re-written in terms of input parameter."

DeleteByFields::usage = "DeleteByFields -> {u,d,...} is an option of WriteFFO. It discards the vertices with the particles listed. "

SelectByFields::usage = "SelectByFields -> {G,t,...} is an option of WriteFFO. It selects the vertices with the particles listed. "

(* Automatize this for an arbitrary representation  *)
MomentumFO=Array[Symbol["p"<>ToString[#]]&,10]
abcIndex=Array[Symbol[FromLetterNumber[#]]&,10]
TagIndex=Array[Symbol["PartTag"<>FromLetterNumber[#]]&,10]
LorentzIndexFO=Array[Symbol["Lor"<>ToString[#]]&,10]
iLorentzIndexFO=Array[Symbol["iLor"<>ToString[#]]&,10]
SpinIndexFO=Array[Symbol["Spin"<>ToString[#]]&,10]
iSpinIndexFO=Array[Symbol["iSpin"<>ToString[#]]&,10]
GluonIndexFO=Array[Symbol["Gluon"<>ToString[#]]&,10]
ColourIndexFO=Array[Symbol["Colour"<>ToString[#]]&,10]
GluonIntIndexFO=Array[Symbol["iGluon"<>ToString[#]]&,10]
ColourIntIndexFO=Array[Symbol["iColour"<>ToString[#]]&,10]


Unprotect[Colour, Spin];
Gluon /: IndicesType[Gluon] = GluonIndexFO;
Colour /: IndicesType[Colour] = ColourIndexFO;
Spin /: IndicesType[Spin] = SpinIndexFO;
Gluon /: IntIndicesType[Gluon] = GluonIntIndexFO;
Colour /: IntIndicesType[Colour] = ColourIntIndexFO;
Protect[Colour, Spin];

AddIndexGG[field_,NIndex_,IndexType_]:=If[ MemberQ[$IndList[field],Index[IndexType] ], IndicesType[IndexType][[NIndex]] ]

IndicesSequence[Field_,index_]:=Sequence@@((AddIndexGG[Field,index,#1]&)/@(Select[($IndList[Field]/.{Index[ind_]:> ind}),Not[#===Lorentz]&]))


HeadFormat[data_List,lineLength_,in_,end_] := Module[{lines, format, pad},
  
  (*format key->value into "Key: value"*)
  
  format[Rule[key_, value_] ] :=  StringJoin[ToString[key], ": ", StringReplace[ToString[value], "{" | "}" -> ""] ];
  
  pad[str_] := StringTake[in<>"   " <> StringPadRight[str, lineLength], lineLength]<>end;
  
  lines = format /@ data;
  
  StringRiffle[pad /@ lines, "\n"]
  ];

(*  Create lists containing the Structure Constants and the Representations of the model  *)

MR$StructureConstantsList={};
MR$RepresentationList={};
Do[

  If[!(Abelian/.M$GaugeGroups[[j,2]]),
  fj=StructureConstant/.M$GaugeGroups[[j,2]];
  MR$StructureConstantsList=AppendTo[MR$StructureConstantsList,fj]
  ];

  If[!(Abelian/.M$GaugeGroups[[j,2]]),
  fj=Representations/.M$GaugeGroups[[j,2]];
  MR$RepresentationList=AppendTo[MR$RepresentationList,fj]
  ];
  
,{j,1,Length[M$GaugeGroups]}];

MR$StructureConstantsList;
MR$RepresentationList=Flatten[MR$RepresentationList,{1,2}];
MR$GeneratorsList=MR$RepresentationList[[;;,1]];

StructureConstantQ[fj_]:=MemberQ[MR$StructureConstantsList,fj]
GeneratorQ[fj_]:=MemberQ[MR$GeneratorsList,fj]







(* Organization of the fields by their indices *)

ClassifyFieldsFO[allFieldsInFR_]:=Block[{currentField,currentIndex},
   
   FO$IndicesInFR={};
   FO$FieldsClassified={};
   Do[
   
      currentField=allFieldsInFR[[j]];
      currentIndex=(Select[#,(Not[#===Index[Lorentz]]&&Not[#===Index[Spin]])&]&)@$IndList[currentField]/.{Index[GG_]:> GG};
   
   
      If[currentIndex=!={},
   
         FO$FieldsClassified=AppendTo[FO$FieldsClassified,{currentIndex[[1]],currentField}];
         If[!MemberQ[FO$IndicesInFR,currentIndex[[1]]],(FO$IndicesInFR=AppendTo[FO$IndicesInFR,currentIndex[[1]]])]
   
      ],
   
   
      {j,1,Length[allFieldsInFR]}];
      FO$FieldsClassified=GatherBy[FO$FieldsClassified,First];
]

FO$WriteGeneral[out_,Fields_,antiFields_]:= Block[{outfile,Bosons,Fermions,allFieldsInFR,xFermions,yFermions,FO$Masses,parameters,extParamFiltered,parametersDefinitions,vectorBosons,scalarBosons,lineLength,ast,gaugeList},

   outfile=out<>"/General.m";

   gaugeList=M$GaugeGroups/.{(x_ == y___) :> {x, y}};
   gaugeList=Flatten/@gaugeList;
   gaugeList=Part[#, 1 ;; 2]&/@gaugeList;
   gaugeList=gaugeList/.{(Abelian -> True) :> Abelian, (Abelian -> False) :> NonAbelian};

   allFieldsInFR=Join[FieldsInFR,antiFieldsInFR];

   Fermions=Select[allFieldsInFR,FermionQ[#]&&Not[GhostFieldQ[#]===True]&];
   Bosons=Select[allFieldsInFR,BosonQ[#]&];
   vectorBosons=Select[allFieldsInFR,VectorFieldQ[#]&];
   scalarBosons=Select[allFieldsInFR,ScalarFieldQ[#]&];

   xFermions=Select[Fermions,(AntiFieldQ[#]===False)&];
   yFermions=Select[Fermions,(AntiFieldQ[#]===True)&];

   FO$Masses=If[(Mass[#] =!= 0), {#, Mass[#]}, Nothing]&/@allFieldsInFR;

   (* parameters=Select[(#1[[1]]&/@substitutionsParam2),Not[(MemberQ[M$Parameters[[;;,1]],#])]&];
   parameters=Join[M$Parameters[[;;,1]],parameters];
   parameters=StringReplace[ToString@FullForm@parameters,{"\\["~~Shortest[x___]~~"]":>x,"List"->"","["->"{","]"->"}", f_ ~~ "[" ~~ x_ ~~ ", " ~~ y_ ~~ "]" :> f <> x <> "x" <> y}]; *)

   parametersDefinitions={};
   parametersDefinitions=Flatten@Append[parametersDefinitions,Variables@Level[#,{-1}]]&/@IParamList[[;;,2]];
   parametersDefinitions=DeleteDuplicates@Flatten@parametersDefinitions;

   parameters=Select[IParamList[[;;,1]],Not[(MemberQ[M$Parameters[[;;,1]],#])]&];
   parameters=Join[M$Parameters[[;;,1]],parameters];
   extParamFiltered=Select[extParam[[;;,1]],Not[(MemberQ[parameters,#])]&];
   parameters=Join[extParamFiltered,parameters];
   parameters=Select[parameters,Not[(MemberQ[FO$Masses[[;; , 2]],#])]&];
   parameters=Select[parameters,Not[(MemberQ[FieldsInFR,#])]&];
   parameters=DeleteDuplicates@Join[parameters,parametersDefinitions];
   parameters=StringReplace[ToString@FullForm@parameters,{"\\["~~Shortest[x___]~~"]":>x,"List"->"","["->"{","]"->"}", f_ ~~ "[" ~~ x_ ~~ ", " ~~ y_ ~~ "]" :> f <> x <> "x" <> y}];
  
   extParam=If[Head[#[[2]]]=!=List,#[[1]]->#[[2]],Nothing]&/@extParam;
   extParam=extParam//.{Rule[a_, b_]:>{a, b}};
   extParam=StringReplace[ToString@FullForm@extParam,{"\\["~~Shortest[x___]~~"]":>x,"List"->"","["->"{","]"->"}"}];

  lineLength=100;
  wsp="                                                                                                                             ";   
  ast="*****************************************************************************************************************************";  

   info=(ToString[#1]<>" : "<>If[Head[#2]===List,StringDrop[StringJoin[(#<>", "&)/@#2],-2],#2]&)@@@M$Information;
   OpenWrite[outfile];
      WriteString[outfile,"("<>StringTake[ast,lineLength]<>")\n"];
      WriteString[outfile,StringTake["(* File automalically generated by FeynRules -"<>wsp,lineLength]<>"*)\n"];
      WriteString[outfile,StringTake["(* date and time of generation : "<>DateString[]<>wsp,lineLength]<>"*)\n"];
      WriteString[outfile,StringTake["(* FeynRules model information : "<>wsp,lineLength]<>"*)\n"];
      WriteString[outfile, HeadFormat[M$Information,lineLength,"(*","*)"]<>"\n"];
      WriteString[outfile,"("<>StringTake[ast,lineLength]<>")\n"];
      WriteString[outfile,"\n"];
      WriteString[outfile,"(*****  This file contains the fields that  appear in the Feynman Rules  *****)\n"];
      WriteString[outfile,"\n"];
      WriteString[outfile,"\n"];
      WriteString[outfile,"AN$ModelName = \""<>out<>"\";\n\n"];
      WriteString[outfile,"AN$GaugeGroups = "<>ToString@gaugeList<>";\n\n"];
      WriteString[outfile,"AN$Fields = "<>ToString@FieldsInFR<>";\n\n"];
      WriteString[outfile,"AN$antiFields = "<>ToString@antiFieldsInFR<>";\n\n"];
      WriteString[outfile,"AN$Fermions = "<>ToString@Fermions<>";\n\n"];
      WriteString[outfile,"AN$Bosons = "<>ToString@Bosons<>";\n\n"];
      WriteString[outfile,"AN$VectorFields = "<>ToString@vectorBosons<>";\n\n"];
      WriteString[outfile,"AN$ScalarFields = "<>ToString@scalarBosons<>";\n\n"];
      WriteString[outfile,"AN$xFermions = "<>ToString@xFermions<>";\n\n"];
      WriteString[outfile,"AN$yFermions = "<>ToString@yFermions<>";\n\n"];
      WriteString[outfile,"AN$NoCouplings = "<>ToString@FO$NoCouplings<>";\n\n"];
      WriteString[outfile,"AN$Parameters = "<>parameters<>";\n\n"];
      WriteString[outfile,"AN$Masses = "<>ToString@FO$Masses<>";\n\n"];
      WriteString[outfile,"AN$IndicesInFR = "<>ToString@FO$IndicesInFR<>";\n\n"];
      WriteString[outfile,"AN$FieldsClassified = "<>ToString@FO$FieldsClassified<>";\n\n"];
      WriteString[outfile,"AN$StructureConstantsList = "<>ToString@MR$StructureConstantsList<>";\n\n"];
      WriteString[outfile,"AN$GeneratorsList = "<>ToString@MR$GeneratorsList<>";\n\n"];
      WriteString[outfile,"AN$ExternalParameters = "<>extParam<>";\n\n"];
   Close[outfile];

]






(***************************)
(*   Writing QGraf model   *)
(***************************)

(* QGWrite uses the function FeynmanRules to find the vertices from the provided Lagrangians which is based in the function WriteUFO, because of this it can share the options of FeynmanRules and WriteUFO. *)





(* This writes the couplings orders *)

FO$WriteVertexQG[FRi_,couplings_,file_]:=Block[{currentCouplingPower,CouplingSpecification,StrCouplingSpecification,VertexHead,StrVertexHead},

  VertexHead=FRi[[1,;;,1]];
  StrVertexHead=StringJoin[ToString/@Riffle[VertexHead,", "]];

(*  currentCouplingPower=Exponent[(FRi[[2]]//.SubstitutionsFRtoForm//.SubstitutionParameters),couplings]; *)
  currentCouplingPower=Exponent[(FRi[[2]]//.SubstitutionsLorentzFRtoForm//.SubstitutionParameters),couplings];
  CouplingSpecification=MapThread[(#1<>"="<>#2)&,{ToString/@couplings,ToString/@currentCouplingPower}];
  StrCouplingSpecification=StringJoin[ToString/@Riffle[CouplingSpecification,", "]];

  WriteString[file,StringTake[wsp,8]<>"["<>StrVertexHead<>"; "<>StrCouplingSpecification<>" ]\n"]

];

(* QGParticleSelection deletes vertices containing the particles given as a list *)

QGParticleSelection[FeynmanRules_, Particles_List] := Block[{NewFR},

  NewFR=FeynmanRules;

  Do[
    NewFR =  Select[NewFR, (MemberQ[Part[#, 1, ;; , 1], Particles[[ii]]] === False) &];
    , {ii, 1, Length[Particles]}];

  Return[NewFR];

];


Options[QGWrite]:=Join[Options[FeynmanRules],Options[WriteUFO],Options[FO$WriteVertices] ];


QGWrite[lagrangians___,out0_, options:OptionsPattern[] ]:=Block[{lags={lagrangians},out=out0,outfile,vertices,ParticlesInVertices,Bosons,Fermions,PropagatorsBosons,PropagatorsFermions,FieldsInFR,particlesDeleted,lineLength,ast,selectbyFields,selFields},


  If[lags=!={} && Head[lags]===Plus,

	(*AppendTo[GenInt$LogFile,"   * Calling FeynmanRules for "<>ToString[Length[lags]]<>" Lagrangians."];*)

    (* Construct options for Feynman Rules ConservedQuantumNumbers is put to False. *)

	  opts=Table[Rule[Options[QGWrite][[iii,1]],OptionValue[Options[QGWrite][[iii,1]]] ],{iii,Length[Options[QGWrite] ]}];
	  opts=FilterRules[opts,First/@Options[FeynmanRules] ];
	  opts=FilterRules[opts,Except[ConservedQuantumNumbers|ScreenOutput] ];
	  opts=Join[opts,{ConservedQuantumNumbers->False,ScreenOutput->False}];

    (* Set FlavorExpand from Automatic to FR$AutoFlavorExpand,and from {} to False*)

	  (*opts=opts/.{Rule[FlavorExpand,Automatic]\[Rule]Rule[FlavorExpand,FR$AutoFlavorExpand]}/.{Rule[FlavorExpand,{}]\[Rule]Rule[FlavorExpand,False]};*)

	  vertices=Join@@(FeynmanRules[#,Sequence@@opts]&/@lags);
  ];

  If[Head[lagrangians]===List,vertices=lagrangians];

  (* vertices=NewFR; *)

  particlesDeleted=OptionValue[DeleteByFields];

  If[particlesDeleted=!=None, vertices=QGParticleSelection[vertices,particlesDeleted] ];

  selectbyFields[list_,fields_]:=Module[{allowed=Flatten@{fields} }, Select[list, SubsetQ[allowed,First/@First[#] ]&] ];

  selFields = OptionValue[SelectByFields];
  If[selFields=!=All,  vertices = selectbyFields[vertices,selFields] ];
 (* Print[selFields];
  Print[vertices]; *)
	ParticlesInVertices=vertices[[;;,1]];


  FO$Couplings=OptionValue[QGCouplings];
  InputParameterList=OptionValue[InternalParameters];

  TransformParameters[InputParameterList];

  (* The relevant fields in the Lagrangian are obtained from the vertices and then classified. Below they are arranged in the format used by QGraf. *)

	Print[" --- QGraf model ---"];

	FieldsInFR=(Select[#,!AntiFieldQ[#]&]&)@DeleteDuplicates@Flatten@ParticlesInVertices[[;;,;;,1]];
	Bosons=Select[FieldsInFR,BosonQ];
	Fermions=Select[FieldsInFR,FermionQ];
	PropagatorsBosons=Inner[List,Bosons,anti/@Bosons,List]//InputForm;
	PropagatorsFermions=Inner[List,Fermions,anti/@Fermions,List]//InputForm;

	If[Sort[FieldsInFR]=!=Sort[Join[Bosons,Fermions]],Print[Style["The union of the identified fermions and bosons does not match the list of fields entering the Feynman rules.",Red]]];

(* 	outfile=FO$PathQG<>"/"<>out<>".qgraf";  *)
   outfile = out<>"QG.qgraf";


  lineLength=79;
  wsp="                                                                                                                             ";   
  ast="*****************************************************************************************************************************";  

	  info=(ToString[#1]<>" : "<>If[Head[#2]===List,StringDrop[StringJoin[(#<>", "&)/@#2],-2],#2]&)@@@M$Information;
	Print["Writing  QGraf model in "<>outfile];
	OpenWrite[outfile];
    WriteString[outfile,"%"<>StringTake[ast,lineLength]<>"\n"];
    WriteString[outfile,StringTake["%* File automalically generated by FeynRules -"<>wsp,lineLength]<>"*\n"];
    WriteString[outfile,StringTake["%* date and time of generation : "<>DateString[]<>wsp,lineLength]<>"*\n"];
    WriteString[outfile,StringTake["%* FeynRules model information : "<>wsp,lineLength]<>"*\n"];
    WriteString[outfile, HeadFormat[M$Information[[1;;1]],lineLength,"%*","*"]<>"\n"];
    WriteString[outfile,"%"<>StringTake[ast,lineLength]<>"\n"];
	  WriteString[outfile,"\n\n"];
	  WriteString[outfile,"%       *****       Propagators       *****\n"];
	  WriteString[outfile,"\n\n"];
	  (WriteString[outfile,StringTake[wsp,8]<>"["<>((StringDelete[#,"}"] &)@(StringDelete[#,"{"] &)@(ToString[InputForm[#]]))<>", + ]\n"]&)/@PropagatorsBosons[[1]];
	  (WriteString[outfile,StringTake[wsp,8]<>"["<>((StringDelete[#,"}"] &)@(StringDelete[#,"{"] &)@(ToString[InputForm[#]]))<>", - ]\n"]&)/@PropagatorsFermions[[1]];
	  WriteString[outfile,"\n\n"];
	  WriteString[outfile,"%       *****       Vertices       ***** \n"];
	  WriteString[outfile,"\n\n"];
    (FO$WriteVertexQG[#,FO$Couplings,outfile]&)/@vertices;
	Close[outfile];
];





(*  Function that writes the amplitude  *)

Options[WriteFFO]:=Join[Options[QGWrite], Options[FO$WriteVertices] ];
  
WriteFFO[out0_,vertices0_, opts: OptionsPattern[] ]:=Block[{out=out0,outfile,vertices=vertices0,FieldsInFR,antiFieldsInFR,nLines,allFieldsInFR,partTag,SubstitutionsPropagators,SubstitutionsPolarization,currentField,currentIndex,currentPropagator,currentPolarization,parameters,extParam,notSubbed,intParam,substitutionsParam,selFields},

   Print["--- FORM FeynRules Output (FFO) v1.0\[Beta]1 ---"];
   Print[Style["Feynman Rules already computed",Orange,Bold] ];

  (* OPTIONS Select and delete by fields *)

   deletebyFields[list_, forb_] := Module[{forbidden = Flatten@{forb} },  Select[list, FreeQ[First/@First[#], Alternatives @@forbidden] &] ];
   selectbyFields[list_,fields_]:=Module[{allowed=Flatten@{fields} }, Select[list, SubsetQ[allowed,First/@First[#] ]&] ];

   forbFields = OptionValue[DeleteByFields];
   If[forbFields=!=None,  vertices = deletebyFields[vertices,forbFields] ];

   selFields = OptionValue[SelectByFields];
   If[selFields=!=All,  vertices = selectbyFields[vertices,selFields] ];

   (* GET FIELDS IN FEYNMAN RULES *)

   FieldsInFR=(Select[#,!AntiFieldQ[#]&]&)@DeleteDuplicates@Flatten@vertices[[;;,1]];
   antiFieldsInFR=anti/@Select[FieldsInFR,Not[#=== anti@#]&];
   allFieldsInFR=Join[FieldsInFR,antiFieldsInFR];
   ClassifyFieldsFO[allFieldsInFR];
   (* Create Directory *)

    parameters={#[[1]],Sequence@@#[[2]]}&/@M$Parameters;
    parameters=Replace[#,{{pp_,___,ParameterType->IE_,___,Value->FF_,___}:>{pp,IE,FF}}]&/@parameters;
    extParam = Cases[parameters, {x___, External, y___} -> {x, y}];

    notSubbed = {\[Alpha]EW, MW, sw2, cw, gs, ee, vev};

    intParam=Cases[parameters,{x___,Internal,y___}->{x,y}];
    substitutionsParam=Select[intParam,Not@MemberQ[notSubbed,#[[1]]]&];
    substitutionsParam=If[Head[#[[2]]]===List,Sequence@#[[2]],#[[1]]->#[[2]]]&/@substitutionsParam;
    substitutionsParam=Flatten@substitutionsParam;

    vertices=vertices//.{f___[Index[Generation, x_], Index[Generation, y_] ] :> f[x, y] };
    vertices=vertices//.substitutionsParam//.FR$RmDblExt;

   If[DirectoryQ[out],DeleteDirectory[out,DeleteContents->True] ];
   CreateDirectory[out];

  (* WriteVertices comes first as it calls the UFO routine *)

   FO$WriteVertices[out,vertices,opts]; (* COUPLINGS file is written in Vertices routine *)

   FO$WritePropagators[out,allFieldsInFR];

   FO$WritePolarizations[out,allFieldsInFR];

   FO$WriteGeneral[out,FieldsInFR,antiFieldsInFR];

   QGWrite[vertices,out<>"/"<>out, opts]; 

   Print["Done!"];

]