Dir="res/";

Group {
  
  bottom = Region[1];
  right = Region[2];
  top = Region[3];
  left = Region[4];
  
  D = Region[100];
  Bnd_D = Region[50];
  DomainDummy = Region[175];
  Extended_D = Region[{D, left, right}]; // Extending support for post-proc
}

Function {
  DefineConstant[ Freq, Val_Rint, Val_Rext ];
  DefineFunction[ ks0, js0, nxe];
  
  I[] = Complex[0, 1];
  Freq = 1e8;
  omega[] = 2*Pi*Freq;
  mu[] = 4*Pi*1e-7;
  eps[] = 8.85e-12;

  I0 = 1;
  a = 0.1;
  h0[] = I0/a * Vector[0,0,1];
  c02[] = mu[] * eps[]; 
  
  // Need for Gauged (electric field)
  Flag_CotreeGauge = 0;
}

Constraint {
  { Name e0 ; Type Assign;
    Case {
      { Region Region[{bottom, right, top}] ; Value 0. ; }
    }
  }
  { Name gauge ; Type Assign;
    Case {
     { Region D ; SubRegion Region [ {bottom, right, top} ] ; Value 0. ; }
    }
  }
  { Name xi0 ; Type Assign;
    Case {
      { Region Region[{bottom, right, top}] ; Value 0. ; }
    }
  }
  { Name h0 ; Type Assign;
    Case {
	  { Region Region[{right}] ; Value Norm[h0[]] ; }
    }
  }
  { Name h_source ;
    Case {
    }
  }

}

FunctionSpace {
  { Name H_e; Type Form1;
    BasisFunction {
      { Name se; NameOfCoef ve; Function BF_Edge; Support Extended_D; Entity EdgesOf[All]; }
    }
    Constraint {
      { NameOfCoef ve; EntityType EdgesOf; NameOfConstraint e0; }
      If(Flag_CotreeGauge)
        { NameOfCoef ve  ; EntityType EdgesOfTreeIn ; EntitySubType StartingOn ;
          NameOfConstraint gauge ; }
      EndIf
    }
  }
  { Name H_xi; Type Form0;
    BasisFunction {
      { Name sn; NameOfCoef vn; Function BF_Node; Support D; Entity NodesOf[All]; }
    }
    Constraint {
      { NameOfCoef vn; EntityType NodesOf; NameOfConstraint xi0; }
    }
  }
  { Name H_h; Type Form1P;
    BasisFunction {
      { Name sn; NameOfCoef hn; Function BF_PerpendicularEdge; Support Extended_D; Entity NodesOf[All]; }
	  { Name sn; NameOfCoef hn2; Function BF_PerpendicularEdge; Support Bnd_D; Entity NodesOf[All]; }
    }
    Constraint {
      { NameOfCoef hn; EntityType NodesOf; NameOfConstraint h0; }
	  { NameOfCoef hn2; EntityType NodesOf; NameOfConstraint h_source; }
    }
  }
}

Jacobian {
  { Name V;
    Case {
      { Region All; Jacobian Vol; }
    }
  }
  { Name S;
    Case {
      { Region All; Jacobian Sur; }
    }
  }
  { Name Jac1D;
    Case {
      { Region All; Jacobian Lin; }
    }
  }
}

Integration {
  { Name I;
    Case {
      { Type Gauss;
        Case {
		  { GeoElement Line; NumberOfPoints  4; }
          { GeoElement Triangle; NumberOfPoints  4; }
	}
      }
    }
  }
}

Formulation{
  // Electric field formulation
  { Name E_Full_wave ; Type FemEquation;
    Quantity {
       { Name hs; Type Local; NameOfSpace H_h; }
	   { Name e; Type Local; NameOfSpace H_e; }
    }
    Equation {
      Galerkin { [ Dof{d e} , {d e} ];
        In D; Integration I; Jacobian V;  }
      Galerkin { DtDtDof [ c02[] * Dof{e} , {e} ];
        In D; Integration I; Jacobian V;  }
	  
	  // Galerkin { [ Normal[] /\ Dof{d e} , {e} ] ;
	  //   In Bnd_D; Integration I; Jacobian S ; }

	  Galerkin { DtDof[ mu[] * Dof{d hs} , {e} ];
        In D; Integration I; Jacobian V;  }
	  Galerkin { [ ( Normal[] /\ Dof{d e} ) , {e} ];
        In Bnd_D; Integration I; Jacobian S;  }  
    }
  }
}

Resolution {
  { Name Analysis;
    System {
	  { Name A; NameOfFormulation E_Full_wave; Type Complex; Frequency Freq; }
    }
    Operation {
	  CreateDir["res/"];
	  Generate[A]; Solve[A]; SaveSolution[A]; PostOperation[FW_e_plots];
    }
  }
}

PostProcessing {
  
  { Name e; NameOfFormulation E_Full_wave;
    PostQuantity {
      /*{ Name hs; Value {
          Term { [ {hs} ]; In D; Jacobian V; } } }*/
	  { Name e; Value {
          Term { [ {e} ]; In D; Jacobian V; } } }
      { Name curle ;  Value{
          Term { [ {d e} ]; In D; Jacobian V; } } }
	  { Name h_from_e ;  Value{
          Local{ [ I[]*(1/mu[])*{d e}/(2*Pi*Freq) ] ; In D; Jacobian V; } } }

      { Name Voltage; Value {
          Integral { [ CompY[ {e} ] ]; In left; Jacobian S; Integration I; } } }
	  { Name L_int; Value {
        Integral { [ 1 ]; In left; Integration I; Jacobian Jac1D; } } }
	  /*{ Name Current; Value {
          Integral { [ CompZ[ {h} ] ]; In left; Integration I; Jacobian S; } } }
	  { Name Z_input; Value {
        Term { Type Global; [ $Input_Volt/(a*$H_int/$L_int) ]; In DomainDummy; } } }*/
    }
  }
}

PostOperation {
  { Name FW_e_plots ; NameOfPostProcessing e;
    Operation {
      //Print[ hs, OnElementsOf D, File StrCat[Dir,"hs.plots"] ] ;
	  Print[ e, OnElementsOf D, File StrCat[Dir,"e.plots"] ] ;
      Print[ curle, OnElementsOf D, File StrCat[Dir,"curle.pos"] ] ;
	  Print[ h_from_e, OnElementsOf D, File StrCat[Dir,"h_from_e.pos"] ] ;
	  //Print[ h, OnLine {{0,0,0}{2,0,0}} {40}, Format SimpleTable, File StrCat[Dir,"Cut_h.dat"] ];
	  Print[ Voltage[left], OnGlobal, Format Table, StoreInVariable $Input_Volt, File StrCat[Dir,"Voltage.txt"] ] ;
	  //Print[ Current[left], OnGlobal, Format Table, StoreInVariable $H_int, File StrCat[Dir,"H_int.txt"] ] ;
	  Print[ L_int[left], OnGlobal, Format Table, StoreInVariable $L_int, File StrCat[Dir,"L_int.txt"] ] ;
	  //Print[ Z_input, OnRegion DomainDummy, Format FrequencyTable, File >> StrCat[Dir,"Z.txt"] ] ;
    }
  }
}

DefineConstant[
  C_ = {"-solve", Name "GetDP/9ComputeCommand", Visible 0}
];
