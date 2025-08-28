

(*****************************)
(*     Write Propagators     *)
(*****************************)

(* UNITARY GAUGE *)

FO$ScalarPropagator[field_]= I Den[p1,m]/.m-> Mass[field];

FO$GoldstonePropagator[field_,Xi_]= 0;

FO$VectorPropagatorMassless[field_,Xi_,ind_]:=-I Times@@FO$Delta@@@Transpose[{Through[ind[1]],Through[ind[2]]}](1 Den[p1,0] FO$Metric[Lor1,Lor2])/.{FF_[GG_[x_],GG_[y_]]:>FF[ ToExpression[ToString[GG]<>ToString[x]],ToExpression[ToString[GG]<>ToString[y]]]};

FO$VectorPropagatorMassive[field_,ind_]:=-I Den[p1,m] Times@@FO$Delta@@@Transpose[{Through[ind[1]],Through[ind[2]]}](FO$Metric[Lor1,Lor2]-(p1[Lor1]p1[Lor2])/m^2)/.m-> Mass[field];

FO$FermionPropagator[field_,ind_]:=(I(FO$Gamma[p1,Spin1,Spin2]+m FO$Gamma[1,Spin1,Spin2])) Den[p1,m] Times@@FO$Delta@@@Transpose[{Through[ind[1]],Through[ind[2]]}]/.{m-> Mass[field],FF_[GG_[x_],GG_[y_]]:>FF[ ToExpression[ToString[GG]<>ToString[x]],ToExpression[ToString[GG]<>ToString[y]]]};

FO$GhostPropagatorMassless[field_,ind_]:=I Den[p1,0] Times@@FO$Delta@@@Transpose[{Through[ind[1]],Through[ind[2]]}]/.{FF_[GG_[x_],GG_[y_]]:>FF[ ToExpression[ToString[GG]<>ToString[x]],ToExpression[ToString[GG]<>ToString[y]]]};

FO$GhostPropagatorMassive[field_,Xi_,ind_]:=0;

(* FEYNMAN GAUGE *)

FO$VectorPropagatorMassiveFeynman[field_,ind_]:=-I Den[p1,m] Times@@FO$Delta@@@Transpose[{Through[ind[1] ],Through[ind[2] ]}](FO$Metric[Lor1,Lor2])/.m->Mass[field];

FO$GoldstonePropagatorFeynman[field_,Xi_]= I Den[p1,m]/.m-> Mass[field];

FO$GhostPropagatorMassiveFeynman[field_,Xi_,ind_]:= -I Den[p1,m]/.m-> Mass[field];



FO$Propagator[field_]:=Which[

  IndicesField=(Select[#,(Not[#===Index[Lorentz]]&&Not[#===Index[Spin]])&]&)@$IndList[field]/.{Index[GG_]:> GG};

  ScalarFieldQ[field]===True&& Not[GoldstoneQ[field]=== True]&& Not[GhostFieldQ[field]=== True],FO$ScalarPropagator[field],
  ScalarFieldQ[field]===True&& GoldstoneQ[field]=== True,FO$GoldstonePropagator[field,Xi],
  VectorFieldQ[field]===True&&Mass[field]=== 0, FO$VectorPropagatorMassless[field,Xi,IndicesField],
  VectorFieldQ[field]===True&&Not[Mass[field]===0], FO$VectorPropagatorMassive[field,IndicesField],
  FermionQ[field]===True&& Not[GhostFieldQ[field]=== True], FO$FermionPropagator[field,IndicesField],
  GhostFieldQ[field]===True&&Mass[field]=== 0, FO$GhostPropagatorMassless[field,IndicesField],
  GhostFieldQ[field]===True&&Not[Mass[field]===0], FO$GhostPropagatorMassive[field,Xi,IndicesField]

]

FO$PropagatorFeynman[field_]:=Which[

  IndicesField=(Select[#,(Not[#===Index[Lorentz]]&&Not[#===Index[Spin]])&]&)@$IndList[field]/.{Index[GG_]:> GG};

  VectorFieldQ[field]===True&&Mass[field]=== 0, FO$VectorPropagatorMassless[field,Xi,IndicesField],
  ScalarFieldQ[field]===True&& GoldstoneQ[field]=== True,FO$GoldstonePropagatorFeynman[field,Xi],
  VectorFieldQ[field]===True&&Not[Mass[field]===0], FO$VectorPropagatorMassiveFeynman[field,IndicesField],
  GhostFieldQ[field]===True&&Not[Mass[field]===0], FO$GhostPropagatorMassiveFeynman[field,Xi,IndicesField]

]


(*  Write propagators *)

FO$WritePropagators[out0_,allFieldsInFR_]:=Block[{out=out0,outfile,SubstitutionsPropagators,currentField,currentPropagator,SubstitutionsFieldsString,SubstitutionsPropatorsString,lineLength,ast,goldstoneBosons,GaugedFields,nonGaugedFields,EWFields,noPropFields},
  
  (* Discriminate fields *)
 
  goldstoneBosons = Select[allFieldsInFR, GoldstoneQ];

  GaugedFields = {A, ghWm, W, ghWp, Z, ghZ, G};
  GaugedFields = Join[goldstoneBosons,GaugedFields];

  noPropFields = {ghWmbar,ghWpbar,ghZbar,ghAbar,ghGbar};

  nonGaugedFields = Complement[allFieldsInFR,GaugedFields];
  nonGaugedFields = Complement[nonGaugedFields,noPropFields];
  nonGaugedFields = DeleteCases[nonGaugedFields, _?(AntiFieldQ)];
  
  EWFields = {G0, GP, A, ghWm, W, ghWp, Z, ghZ};

  (* Substitutions to provide the indices that the fields should have at the propagator level  *)
  
  SubstitutionsPropagators={
  
    Field_/;(Length@$IndList[Field]<2&&ScalarFieldQ[Field]&&Not[GhostFieldQ[Field]=== True]):> Field [TagIndex[[1]],MomentumFO[[1]],LorentzIndexFO[[1]],LorentzIndexFO[[2]]],
  
    Field_/;(Length@$IndList[Field]<1&&GhostFieldQ[Field]):> Field [TagIndex[[1]],MomentumFO[[1]],LorentzIndexFO[[1]],LorentzIndexFO[[2]]],
  
    Field_/;(Length@$IndList[Field]<2&&VectorFieldQ[Field]):> Field [TagIndex[[1]],MomentumFO[[1]],LorentzIndexFO[[1]],LorentzIndexFO[[2]]],
       
    Field_/;(Length@$IndList[Field]<2&&FermionQ[Field]&&AddIndexGG[Field,1,Gluon]=== Null ):> Field [TagIndex[[1]],MomentumFO[[1]],SpinIndexFO[[1]],SpinIndexFO[[2]]],
  
    Field_/;(AddIndexGG[Field,1,Gluon]=!= Null):> Field [TagIndex[[1]],MomentumFO[[1]],LorentzIndexFO[[1]],LorentzIndexFO[[2]],AddIndexGG[Field,1,Gluon],AddIndexGG[Field,2,Gluon] ],
  
    Field_/;(AddIndexGG[Field,1,Colour]=!= Null):> Field [TagIndex[[1]],MomentumFO[[1]],SpinIndexFO[[1]],SpinIndexFO[[2]],AddIndexGG[Field,1,Colour],AddIndexGG[Field,2,Colour] ]
  };

  (* String substitutions *)
  
  SubstitutionsFieldsString = {

    ("PartTag"~~char:LetterCharacter):>"?"<>char,
    ("p"~~char:DigitCharacter):>"p"<>char<>"?",
    ("Lor"~~char:DigitCharacter):>"Lor"<>char<>"?",
    ("Spin"~~char:DigitCharacter):>"Spin"<>char<>"?",
    ("Gluon"~~char:DigitCharacter):>"Gluon"<>char<>"?",
    ("Colour"~~char:DigitCharacter):>"Colour"<>char<>"?",
    "["-> "(",
    "]"-> ")"
  };
  
  SubstitutionsPropatorsString = {

    "FO$Metric"-> "d_",
    "FO$Gamma"-> "GammaM",
    "FO$Delta"-> "d_",
    ("QuestionMark"~~char:LetterCharacter):>"?"<>char,
    "["-> "(",
    "]"-> ")"
  };
  
  outfile=out<>"/Propagators"<>".frm";
  Print["Writing  propagators in "<>outfile];

  lineLength=100;
  wsp="                                                                                                                             ";   
  ast="*****************************************************************************************************************************";  
  OpenWrite[outfile];
    WriteString[outfile,StringTake[ast,lineLength+1]<>"\n"];
    WriteString[outfile,StringTake["* File automalically generated by FeynRules -"<>wsp,lineLength]<>"*\n"];
    WriteString[outfile,StringTake["* date and time of generation : "<>DateString[]<>wsp,lineLength]<>"*\n"];
    WriteString[outfile,StringTake["* FeynRules model information : "<>wsp,lineLength]<>"*\n"];
    WriteString[outfile, HeadFormat[M$Information,lineLength,"*","*"]<>"\n"];
    WriteString[outfile,StringTake[ast,lineLength+1]<>"\n"];
    WriteString[outfile,"\n"];
    WriteString[outfile,"\n"];
    WriteString[outfile,"*****  This file contains the Propagators   *****\n"];
    WriteString[outfile,"\n"];
    WriteString[outfile,"\n"];
    WriteString[outfile,"#procedure Propagators\n"];
    WriteString[outfile,".sort\n\n"];
    WriteString[outfile,"*********************************************\n"];
    WriteString[outfile,"*    Propagators non dependent on Gauge     *\n"];
    WriteString[outfile,"*********************************************\n\n"];
    WriteString[outfile,"id   pro(ff?Field(x?,y?partTagInt[n],?a)) = pro(ff(x,y,?a,IntLor[n]));\n"];
    WriteString[outfile,"id   pro(ff?Fermion(x?,y?partTagInt[n],p1?,iLor1?,Spin1?,iLor2?)) = pro(ff(x,y,p1,Spin1,IntSpin[n]));\n"];
    WriteString[outfile,"id   pro(ff?Fermion(x?,y?partTagInt[n],p1?,iLor1?,Spin1?,Colour1?,iLor2?)) = pro(ff(x,y,p1,Spin1,IntSpin[n],Colour1,IntColour[n]));\n"];
    WriteString[outfile,"id   pro(ff?Gluoned(x?,y?partTagInt[n],p1?,iLor1?,Gluon1?,iLor2?)) = pro(ff(x,y,p1,iLor1,iLor2,Gluon1,IntGluon[n]));\n"];
    WriteString[outfile,"id   pro(ff?Gluoned(x?,y?partTagInt[n],0,iLor1?,Gluon1?,iLor2?)) = pro(ff(x,y,0,iLor1,iLor2,Gluon1,IntGluon[n]));\n"];
    WriteString[outfile,"\n"];
    WriteString[outfile,"\n"];
    WriteString[outfile,"*****  Propagators  *****\n"];
    WriteString[outfile,"\n"];
    Do[

      currentField=(StringReplace[#, SubstitutionsFieldsString ]&)@(ToString[InputForm[#] ]&)@(nonGaugedFields[[j]]/.SubstitutionsPropagators);
      currentPropagator=(StringReplace[#, SubstitutionsPropatorsString ]&)@(ToString[InputForm[#] ])&@(FO$Propagator[nonGaugedFields[[j]] ]//.{Complex[a_,b_]:>a+ im b});

      WriteString[outfile,"id   pro("<>currentField<>") = "<>currentPropagator<>";\n"]

    ,{j,1,Length[nonGaugedFields]} ];
    WriteString[outfile,"\n"];
    WriteString[outfile,"#endprocedure Propagators\n\n"];
    WriteString[outfile,"*********************************************\n"];
    WriteString[outfile,"*       Propagators QCD RXi Gauges          *\n"];
    WriteString[outfile,"*********************************************\n\n"];
    WriteString[outfile,"*****  Propagators QCD Feynman Gauge  *****\n\n"];
    WriteString[outfile,"#procedure PropagatorsQCDFeynmanGauge\n"];
    WriteString[outfile,".sort\n\n"];
    WriteString[outfile,"id   pro(G(?a, p1?, Lor1?, Lor2?, Gluon1?, Gluon2?)) = -im*Den(p1, 0)*d_(Gluon1, Gluon2)*(d_(Lor1, Lor2));\n\n"];
    WriteString[outfile,"#endprocedure PropagatorsQCDFeynmanGauge\n\n"];
    WriteString[outfile,"*****  Propagators RXi Gauge  *****\n\n"];
    WriteString[outfile,"#procedure PropagatorsQCDRXiGauge\n"];
    WriteString[outfile,".sort\n\n"];
    WriteString[outfile,"id   pro(G(?a, p1?, Lor1?, Lor2?, Gluon1?, Gluon2?)) = -im*Den(p1, 0)*d_(Gluon1, Gluon2)*(d_(Lor1, Lor2) - (1 - RXi)*(p1(Lor1)*p1(Lor2))*Den(p1,0));\n\n"];
    WriteString[outfile,"#endprocedure PropagatorsQCDRXiGauge\n\n"];
    WriteString[outfile,"*********************************************\n"];
    WriteString[outfile,"*          Propagators EW Gauges            *\n"];
    WriteString[outfile,"*********************************************\n\n"];
    WriteString[outfile,"*****  Propagators EW Feynman Gauges  *****\n\n"];
    WriteString[outfile,"#procedure PropagatorsEWFeynmanGauge\n"];
    WriteString[outfile,".sort\n\n"];
    Do[

      currentField=(StringReplace[#, SubstitutionsFieldsString ]&)@(ToString[InputForm[#] ]&)@(EWFields[[j]]/.SubstitutionsPropagators);
      currentPropagator=(StringReplace[#, SubstitutionsPropatorsString ]&)@(ToString[InputForm[#] ])&@(FO$PropagatorFeynman[EWFields[[j]] ]//.{Complex[a_,b_]:>a+ im b});

      WriteString[outfile,"id   pro("<>currentField<>") = "<>currentPropagator<>";\n"]

    ,{j,1,Length[EWFields]} ];
    WriteString[outfile,"\n"];
    WriteString[outfile,"#endprocedure PropagatorsEWUnitaryGauge\n\n"];
    WriteString[outfile,"*****  Propagators EW Unitary Gauges  *****\n\n"];
    WriteString[outfile,"#procedure PropagatorsEWUnitaryGauge\n"];
    WriteString[outfile,".sort\n\n"];
    Do[

      currentField=(StringReplace[#, SubstitutionsFieldsString ]&)@(ToString[InputForm[#] ]&)@(EWFields[[j]]/.SubstitutionsPropagators);
      currentPropagator=(StringReplace[#, SubstitutionsPropatorsString ]&)@(ToString[InputForm[#] ])&@(FO$Propagator[EWFields[[j]] ]//.{Complex[a_,b_]:>a+ im b});

      WriteString[outfile,"id   pro("<>currentField<>") = "<>currentPropagator<>";\n"]

    ,{j,1,Length[EWFields]} ];
    WriteString[outfile,"\n"];
    WriteString[outfile,"#endprocedure PropagatorsEWUnitaryGauge\n\n"];
  Close[outfile];
]