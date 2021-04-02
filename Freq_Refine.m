function Max_Ref_Freq = Freq_Refine(best_Cosine_Coe,Val,Max_Ref_Freq)
best_Cosine_Coe = best_Cosine_Coe(:);
Val = Val(:)';
Max_Ref_Freq_w = Max_Ref_Freq(:);
m = length(best_Cosine_Coe)-1;
for k = 1:5
   H1 = -sin(Max_Ref_Freq_w*Val) * (Val'.*best_Cosine_Coe);
   H2 = -cos(Max_Ref_Freq_w*Val) * ((Val.^2)'.*best_Cosine_Coe);
   Max_Ref_Freq_w = Max_Ref_Freq_w - H1./H2;
end
Max_Ref_Freq(:) = Max_Ref_Freq_w;
