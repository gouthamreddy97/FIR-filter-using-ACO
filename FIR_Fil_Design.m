function FIR_F=FIR_Fil_Design(Fil_Order,Discre_Word_Length,Fp_value,Fst_value,No_Itera,No_Ants);
clc
clear all
close all
close all
delta_T=0.3;
if ~exist('Fil_Order'  ,'var');  
Fil_Order_Prompt = 'Enter the N value: ';
Fil_Order = input(Fil_Order_Prompt);
end
if ~exist('Discre_Word_Length' ,'var');
Word_Length_Prompt = 'Enter the p value: ';
Discre_Word_Length = input(Word_Length_Prompt);
end
if ~exist('Fp_value' ,'var');
Fp_Prompt = 'Enter the fp value: ';
Fp_value = input(Fp_Prompt);
end
if ~exist('Fst_value','var'); 
Fst_Prompt = 'Enter the fst value: ';
Fst_value = input(Fst_Prompt);   
end
if ~exist('q'  ,'var');  womega=16;       end
No_Itera_Prompt = 'Enter the requird number of Iterations: ';
No_Itera =input(No_Itera_Prompt);
if ~exist('No_Ants'  ,'var')
  No_Ants_Prompt = 'Enter the requird number of Ants: ';
No_Ants=input(No_Ants_Prompt);
end
Fil_Inputs  = randn(1,No_Ants); 
Fst_value=0.6;
Passband_mag  = pass_Mag_Resp(Fil_Order,Fp_value);    
n  = 0.1.*Fil_Inputs;
fd=fdesign.lowpass('n,fp,fst,ast',Fil_Order,0.3,Fst_value,No_Itera);
Fil_Design_Weights=fdesign.lowpass('n,fp,fst,ast',Fil_Order,delta_T,Fst_value,No_Itera);
Filter_Mag_Res=Filter_Design(Fil_Design_Weights, 'equiripple');
Fil_Mag_wei=Filter_Mag_Res.Numerator;
Filter_Shape_Pass=1E-1;
Filter_Shape_Stop=1E-2;
FIR=FIR_Bands(Discre_Word_Length,(Fp_value+Fst_value)/2,Filter_Shape_Pass,Filter_Shape_Stop,Fp_value,Fst_value,1E4);
Freq_Norm=linspace(0,2*(No_Ants-1)/No_Ants,No_Ants);
Path=find(Freq_Norm>Fst_value & Freq_Norm<1);
Scale=2^(womega-1);
Rescale_Mag=Scale*Fil_Mag_wei/max(abs(Fil_Mag_wei))*(Scale-1)/Scale;
Passband_Mag_ACO_Pro  = filter(Passband_mag,1,Fil_Inputs);
delay = 1;              
Pass_Ha = adaptfilt.dlms(Fil_Order,(Discre_Word_Length*10^-4),1,delay);
[Passband_Mag_GA,Passband_Mag_ACO_Conv] = filter(Pass_Ha,Fil_Inputs,Passband_Mag_ACO_Pro);
Passband_Mag_OPT = filter(Passband_mag,1,Fil_Inputs)+n;
k=size(Fil_Inputs);
Sum=sum(Rescale_Mag);
Rounded_Mag=round(Rescale_Mag);
Nearest_Path=Find_Path(Rounded_Mag,No_Ants,Path,Sum);
Next_Nearest_Path=Nearest_Path;
for i=1:100
    Fil_Coef=Fil_Shape_coe(FIR,Rescale_Mag);
    Final_Nearest_Path=Find_Path(Fil_Coef,No_Ants,Path,Sum);
    if Final_Nearest_Path<Next_Nearest_Path 
        Next_Nearest_Path=Final_Nearest_Path;
        FIR_F=Fil_Coef;
    end
end
if nargout==0
  h_bb=freqz(FIR/FIR(1),1,pi*Freq_Norm);
  GA_Coe=round(Rescale_Mag);
  GA=freqz(GA_Coe/Sum,1.2,pi*Freq_Norm);
  ACO_Conv=freqz(FIR_F/(Sum-1),1,pi*Freq_Norm);
  Fil_Mag_Wei=Filter_Mag_Res.Numerator;
  Optimal_Fun_Coe=Scale*Fil_Mag_Wei/max(abs(Fil_Mag_Wei))*(Scale-1)/Scale;
  Optimal_Mag=freqz(Optimal_Fun_Coe/Sum,1,pi*Freq_Norm);
  ACO_Pro_Fun_Coe=Scale*Fil_Mag_Wei/max(abs(Fil_Mag_Wei))*(Scale-1)/Scale;
  ACO_Pro_Coe=Scale*ACO_Pro_Fun_Coe/max(abs(ACO_Pro_Fun_Coe))*(Scale-1)/Scale;
  ACO_Prop=freqz(ACO_Pro_Coe/(Sum-1),0.9,pi*Freq_Norm);
  figure();
  plot(Freq_Norm,(20*log10(abs(GA))),'g');
  hold on
  plot(Freq_Norm,(20*log10(abs([ACO_Conv]))),'m');
  hold on
  plot(Freq_Norm,(20*log10(abs([ACO_Prop]))),'r','linewidth',2);
  hold on
  plot(Freq_Norm,(20*log10(abs([Optimal_Mag]))),'.k');
  hold on
  axis([0 1 -200 50]);
  grid on;
  title('Magnitude Response');
  xlabel('Normalized Frequency');
  ylabel('Magnitude (dB)');
  legend( 'GA', 'ACO(Conv)', 'ACO(Prop)','Optimal');
  grid on
  figure
 plot(1:k(1,2),Passband_Mag_ACO_Pro,'r');
 hold on
 plot(1:k(1,2),Passband_Mag_GA,'b');
 hold on
 plot(1:k(1,2),Passband_Mag_ACO_Conv,'g');
 hold on
 plot(1:k(1,2),Passband_Mag_OPT,'c');
title('Passband Magnitude Response');
legend('ACO(Prop)','GA','ACO(Conc)','Optimal');
xlabel('Normalized Frequency'); ylabel('Magnitude (dB)');
grid on

end


