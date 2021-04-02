function Nearest_Path=Find_Path(Rounded_Mag,No_Ants,Path,Sum)
FIR_Co=fft([Rounded_Mag zeros(1,No_Ants-length(Rounded_Mag))]);
Near_freq=FIR_Co(Path);
Nearest_Path=10*log10(max(Near_freq.*conj(Near_freq))/(Sum*Sum));