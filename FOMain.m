<< "FOVertices.m";
<< "FOPolarizations.m";
<< "FOPropagators.m";
<< "FOCouplings.m";


QGCouplings::usage = "QGCouplings -> { \[Alpha]1, \[Alpha]2,... } is an option of QGWrite and chooses the relevant couplings in the QGraf model, which allows diagram filtering. "

InternalParameters::usage = "InternalParameters is an option of QGWrite that specifies the internal parameters that should be re-written in terms of input parameter."

DeleteParticles::usage = "DeleteParticles -> {u,d,...} is an option of QGWrite. It discards the vertices with the particles listed. "

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

FO$WriteGeneral[out_,Fields_,antiFields_]:= Block[{outfile,Bosons,Fermions,allFieldsInFR,xFermions,yFermions,FO$Masses,parameters,extParamFiltered,parametersDefinitions,vectorBosons,scalarBosons},

   outfile=out<>"/General.m";

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
   extParam=StringReplace[ToString@FullForm@extParam,{"\\["~~Shortest[x___]~~"]":>x,"List"->"","["->"{","]"->"}"}];

   wsp="                                                                                                                                ";      
   info=(ToString[#1]<>" : "<>If[Head[#2]===List,StringDrop[StringJoin[(#<>", "&)/@#2],-2],#2]&)@@@M$Information;
   OpenWrite[outfile];
      WriteString[outfile,"(*****  This file contains the fields that  appear in the Feynman Rules  *****)\n"];
      WriteString[outfile,"\n"];
      WriteString[outfile,"\n"];
      WriteString[outfile,"FO$ModelName = \""<>out<>"\";\n\n"];
      WriteString[outfile,"FO$Fields = "<>ToString@FieldsInFR<>";\n\n"];
      WriteString[outfile,"FO$antiFields = "<>ToString@antiFieldsInFR<>";\n\n"];
      WriteString[outfile,"FO$Fermions = "<>ToString@Fermions<>";\n\n"];
      WriteString[outfile,"FO$Bosons = "<>ToString@Bosons<>";\n\n"];
      WriteString[outfile,"FO$VectorFields = "<>ToString@vectorBosons<>";\n\n"];
      WriteString[outfile,"FO$ScalarFields = "<>ToString@scalarBosons<>";\n\n"];
      WriteString[outfile,"FO$xFermions = "<>ToString@xFermions<>";\n\n"];
      WriteString[outfile,"FO$yFermions = "<>ToString@yFermions<>";\n\n"];
      WriteString[outfile,"FO$NoCouplings = "<>ToString@FO$NoCouplings<>";\n\n"];
      WriteString[outfile,"FO$Parameters = "<>parameters<>";\n\n"];
      WriteString[outfile,"FO$Masses = "<>ToString@FO$Masses<>";\n\n"];
      WriteString[outfile,"FO$IndicesInFR = "<>ToString@FO$IndicesInFR<>";\n\n"];
      WriteString[outfile,"FO$FieldsClassified = "<>ToString@FO$FieldsClassified<>";\n\n"];
      WriteString[outfile,"MR$StructureConstantsList = "<>ToString@MR$StructureConstantsList<>";\n\n"];
      WriteString[outfile,"MR$GeneratorsList = "<>ToString@MR$GeneratorsList<>";\n\n"];
      WriteString[outfile,"FO$ExternalParameters = "<>extParam<>";\n\n"];
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

QGParticleSelection[FeynmanRules_, Particles_List] := Block[{},

  NewFR=FeynmanRules;

  Do[
    NewFR =  Select[NewFR, (MemberQ[Part[#, 1, ;; , 1], Particles[[ii]]] === False) &];
    , {ii, 1, Length[Particles]}];

];




Options[QGWrite]:=Join[Options[FeynmanRules],Options[WriteUFO], {QGCouplings -> {},InternalParameters-> {},DeleteParticles->{}}];


QGWrite[lagrangians___,out0_, options:OptionsPattern[]]:=Block[{lags={lagrangians},out=out0,outfile,vertices,ParticlesInVertices,Bosons,Fermions,PropagatorsBosons,PropagatorsFermions,FieldsInFR,particlesDeleted},


  If[lags=!={} && Head[lags]===Plus,

	(*AppendTo[GenInt$LogFile,"   * Calling FeynmanRules for "<>ToString[Length[lags]]<>" Lagrangians."];*)

    (* Construct options for Feynman Rules ConservedQuantumNumbers is put to False. *)

	  opts=Table[Rule[Options[QGWrite][[iii,1]],OptionValue[Options[QGWrite][[iii,1]]]],{iii,Length[Options[QGWrite]]}];
	  opts=FilterRules[opts,First/@Options[FeynmanRules]];
	  opts=FilterRules[opts,Except[ConservedQuantumNumbers|ScreenOutput]];
	  opts=Join[opts,{ConservedQuantumNumbers->False,ScreenOutput->False}];

    (* Set FlavorExpand from Automatic to FR$AutoFlavorExpand,and from {} to False*)

	  (*opts=opts/.{Rule[FlavorExpand,Automatic]\[Rule]Rule[FlavorExpand,FR$AutoFlavorExpand]}/.{Rule[FlavorExpand,{}]\[Rule]Rule[FlavorExpand,False]};*)

	  vertices=Join@@(FeynmanRules[#,Sequence@@opts]&/@lags);
  ];

  If[Head[lagrangians]===List,vertices=lagrangians];


  particlesDeleted=OptionValue[DeleteParticles];
  QGParticleSelection[vertices,particlesDeleted];
  vertices=NewFR;

	ParticlesInVertices=vertices[[;;,1]];


  FO$Couplings=OptionValue[QGCouplings];
  InputParameterList=OptionValue[InternalParameters];

  TransformParameters[InputParameterList];

	Print["Feynman Rules computed"];

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

	wsp="                                                                                                                             ";      
	  info=(ToString[#1]<>" : "<>If[Head[#2]===List,StringDrop[StringJoin[(#<>", "&)/@#2],-2],#2]&)@@@M$Information;
	Print["Writing  QGraf model in "<>outfile];
	OpenWrite[outfile];
	  WriteString[outfile,"%*******************************************************************\n"];
	  WriteString[outfile,"%* Model automalically generated by FeynRules -"<>StringTake[wsp,20]<>"*\n"];
	  WriteString[outfile,"%* date and time of generation : "<>StringTake[DateString[]<>wsp,34]<>"*\n"];
	  WriteString[outfile,"%* FeynRules model information : "<>StringTake[wsp,34]<>"*\n"];
	  WriteString[outfile,"%*******************************************************************\n"];
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
	Print["done"];
];





(*  Function that writes the amplitude  *)

Options[WriteFFO]:=Options[QGWrite];
  
WriteFFO[out0_,vertices0_, opts: OptionsPattern[]]:=Block[{out=out0,outfile,vertices=vertices0,FieldsInFR,antiFieldsInFR,nLines,allFieldsInFR,partTag,SubstitutionsPropagators,SubstitutionsPolarization,currentField,currentIndex,currentPropagator,currentPolarization},

   FieldsInFR=(Select[#,!AntiFieldQ[#]&]&)@DeleteDuplicates@Flatten@vertices[[;;,1]];
   antiFieldsInFR=anti/@Select[FieldsInFR,Not[#=== anti@#]&];
   allFieldsInFR=Join[FieldsInFR,antiFieldsInFR];
   ClassifyFieldsFO[allFieldsInFR];
   (* Create Directory *)

   If[DirectoryQ[out],DeleteDirectory[out,DeleteContents->True]];
   CreateDirectory[out];

  (* FO$WriteGeneral[out,FieldsInFR,antiFieldsInFR]; *)

  (* WriteVertices has to be first as it calls the UFO routine *)

   FO$WriteVertices[out,vertices,opts];
Abort[];
   (* FO$WriteCouplings[out,opts] *)

   FO$WritePropagators[out,allFieldsInFR];

   FO$WritePolarizations[out,allFieldsInFR];

   FO$WriteGeneral[out,FieldsInFR,antiFieldsInFR];

   QGWrite[vertices,out<>"/"<>out, opts]; 

]