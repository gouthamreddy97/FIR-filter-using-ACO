function Fil_Coef=Fil_Shape_coe(FIR,Rescale_Mag)
Length_L=length(Rescale_Mag);
if mod(Length_L,2)==0
    Rescale_Mag=Rescale_Mag((Length_L+2)/2:Length_L);
else
    Rescale_Mag=Rescale_Mag((Length_L+1)/2:Length_L);
end

FIR=FIR/FIR(1);
FIR=FIR(2:length(FIR));
C=(rand(size(FIR))-0.5);
Fil_Coef=zeros(size(FIR));
for i=1:length(Rescale_Mag)
    y0=sum(FIR.*C)+Rescale_Mag(i);
    yr=round(y0);    
    Fil_Coef(i)=yr;
    C=[yr-y0 C(1:length(C)-1)];
end
if mod(Length_L,2)==0 
    Fil_Coef=[fliplr(Fil_Coef) Fil_Coef];
else
    Fil_Coef=[fliplr(Fil_Coef) Fil_Coef(2:length(Fil_Coef))];
end      
