# Stage 1 (GAP): character tables of GL2(Z/m), m = 8, 9, 5.
# Output tables_<m>.txt:  REP a b c d size   (one line per conjugacy class, GAP class order)
#                         L <conductor>
#                         CHAR v_1 ... v_k    (each value as comma-separated coefficients w.r.t. 1,z,...,z^(L-1), z=E(L))
for m in [8,9,5] do
  R:=Integers mod m;; G:=GL(2,R);; iso:=IsomorphismPermGroup(G);; P:=Image(iso);;
  T:=CharacterTable(P);; irr:=Irr(T);; cls:=ConjugacyClasses(T);;
  L:=Lcm(List(Flat(List(irr,ValuesOfClassFunction)),Conductor));;
  out:=OutputTextFile(Concatenation("tables_",String(m),".txt"),false);; SetPrintFormattingStatus(out,false);;
  for c in cls do r:=PreImagesRepresentative(iso,Representative(c));
    AppendTo(out,"REP ",Int(r[1][1])," ",Int(r[1][2])," ",Int(r[2][1])," ",Int(r[2][2])," ",Size(c),"\n"); od;
  AppendTo(out,"L ",L,"\n");
  for chi in irr do
    AppendTo(out,"CHAR ",JoinStringsWithSeparator(List(ValuesOfClassFunction(chi),v->JoinStringsWithSeparator(List(CoeffsCyc(v,L),String),","))," "),"\n");
  od;
  CloseStream(out); Print("m=",m,": ",Length(irr)," irreducibles, L=",L,"\n");
od;
QUIT;
