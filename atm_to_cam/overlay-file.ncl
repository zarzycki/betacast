begin

maxLev = 100.

print("Base_file is: "+base_file)
print("top_file is: "+top_file)
print("maxLev is: "+maxLev)

f1 = addfile(base_file,"w")
f2 = addfile(top_file,"r")

T_base=f1->T
Q_base=f1->Q
U_base=f1->U
V_base=f1->V
PS_base=f1->PS

T_rgm=f2->T
Q_rgm=f2->Q
U_rgm=f2->U
V_rgm=f2->V
PS_rgm=f2->PS

; This only overlays data in the troposphere because regional models generally low-topped
T_rgm(:,{0:maxLev},:) = T_rgm@_FillValue
Q_rgm(:,{0:maxLev},:) = Q_rgm@_FillValue
U_rgm(:,{0:maxLev},:) = U_rgm@_FillValue
V_rgm(:,{0:maxLev},:) = V_rgm@_FillValue

T_base = where(.not.ismissing(T_rgm),T_rgm,T_base)
Q_base = where(.not.ismissing(Q_rgm),Q_rgm,Q_base)
U_base = where(.not.ismissing(U_rgm),U_rgm,U_base)
V_base = where(.not.ismissing(V_rgm),V_rgm,V_base)
PS_base = where(.not.ismissing(PS_rgm),PS_rgm,PS_base)

print("Number of missing values in top_file is: "+num(ismissing(PS_rgm)))

f1->T = T_base
f1->Q = Q_base
f1->U = U_base
f1->V = V_base
f1->PS = PS_base

end
