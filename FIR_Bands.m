function [H_Val,best_Cosine_Coe] = FIR_Bands(Discre_Word_Length,W_O,Filter_Shape_Pass,Filter_Shape_Stop,varargin)
Low_Pass = 0;
High_Pass = 1;
Type_Filter = Low_Pass;
Edge_Val = 0;
Filter_Weights = 0;
Discre_Word_Length = signal.internal.sigcasttofloat(Discre_Word_Length,'double','FIR_Bands','N','allownumeric');
W_O = signal.internal.sigcasttofloat(W_O,'double','FIR_Bands','WO',...
  'allownumeric');
Filter_Shape_Pass = signal.internal.sigcasttofloat(Filter_Shape_Pass,'double','FIR_Bands','DEVP',...
  'allownumeric');
Filter_Shape_Stop = signal.internal.sigcasttofloat(Filter_Shape_Stop,'double','FIR_Bands','DEVS',...
  'allownumeric');
if Type_Filter == Low_Pass
    Boundary_Up = [1+Filter_Shape_Pass  Filter_Shape_Stop];
    Boundary_Low = [1-Filter_Shape_Pass -Filter_Shape_Stop];
end
[Discre_Word_Length,M1,M2,M_Obj] = Fil_order_Validatoin_Fun(Discre_Word_Length,1,Type_Filter);
Discre_Word_Length = Discre_Word_Length+1;  
if rem(Discre_Word_Length,2) == 1
    P_Parity = 1;
else
    P_Parity = 0;
end
W_O = W_O*pi;
if Filter_Weights
    Passband_Weights = Passband_Weights*pi;
    Stopband_Weights = Stopband_Weights*pi;
end
if Edge_Val, WT_Val = WT_Val*pi; end
Word_Length_Ceil = 2^ceil(log2(5*Discre_Word_Length));
Upper_Boundary = [Boundary_Up(1)*ones((round(W_O*Word_Length_Ceil/pi)),1); Boundary_Up(2)*ones(Word_Length_Ceil+1-(round(W_O*Word_Length_Ceil/pi)),1)];
Lower_Boundary = [Boundary_Low(1)*ones((round(W_O*Word_Length_Ceil/pi)),1); Boundary_Low(2)*ones(Word_Length_Ceil+1-(round(W_O*Word_Length_Ceil/pi)),1)];
if P_Parity  == 1
    median = (Discre_Word_Length-1)/2;
    if Filter_Weights
        Coe = 2*[Passband_Weights/sqrt(2); (sin(Passband_Weights*(1:median))./(1:median))']/pi;
    else
        Coe = 2*[W_O/sqrt(2); (sin(W_O*(1:median))./(1:median))']/pi;
    end
    if Type_Filter == High_Pass
        Coe = -Coe;
        Coe(1) = Coe(1)+sqrt(2);
    end
    Zeros = zeros(2*Word_Length_Ceil-Discre_Word_Length,1);
    Val = 0:median;
    T = 1 - 1/sqrt(2);
    Num_of_Param = median+1;	
else
    median = Discre_Word_Length/2;
    Val = (1:median)-1/2;
        if Filter_Weights
        Coe = (2*sin(Passband_Weights*Val)./(Val*pi))';
    else
        Coe = (2*sin(W_O*Val)./(Val*pi))';
    end
    Zeros = zeros(4*Word_Length_Ceil,1);
    T = 0;
    Num_of_Param = median;		
end
if Filter_Weights

else
    best_Cosine_Coe = Coe;
end
SN = 1e-8;              
while 1

    if P_Parity == 1
        H_Val = fft([best_Cosine_Coe(1)*sqrt(2); best_Cosine_Coe(2:median+1); Zeros; best_Cosine_Coe(median+1:-1:2)]);
        H_Val = real(H_Val(1:Word_Length_Ceil+1)/2);
    else
        Zeros(2:2:2*median) = best_Cosine_Coe;
        Zeros(4*Word_Length_Ceil-2*median+2:2:4*Word_Length_Ceil) = best_Cosine_Coe(median:-1:1);
        H_Val = fft(Zeros);
        H_Val = real(H_Val(1:Word_Length_Ceil+1)/2);
    end
    Maxima_H = local_max(H_Val);
    Minima_H = local_max(-H_Val);
    if P_Parity == 0
        Length_Maxima = length(Maxima_H);
        if Maxima_H(Length_Maxima) == Word_Length_Ceil+1, Maxima_H(Length_Maxima) = []; end
        Length_Minim = length(Minima_H);
        if Minima_H(Length_Minim) == Word_Length_Ceil+1, Minima_H(Length_Minim) = []; end
    end
    Max_Ref_Freq = (Maxima_H-1)*pi/Word_Length_Ceil;      
    Min_Ref_Freq = (Minima_H-1)*pi/Word_Length_Ceil;
    Max_Ref_Freq = Freq_Refine(best_Cosine_Coe,Val,Max_Ref_Freq);  
    Min_Ref_Freq = Freq_Refine(best_Cosine_Coe,Val,Min_Ref_Freq);
    if Edge_Val
        if Type_Filter == Low_Pass
            if edge_type == PASS
                Key_W = max(Max_Ref_Freq(Max_Ref_Freq<W_O));
                if WT_Val > Key_W
                    Minima_H = [Minima_H; 1];
                    Min_Ref_Freq = [Min_Ref_Freq; WT_Val];
                end
            else
                Key_W = min(Min_Ref_Freq(Min_Ref_Freq>W_O));
                if WT_Val < Key_W
                    Maxima_H = [Maxima_H; Word_Length_Ceil];
                    Max_Ref_Freq = [Max_Ref_Freq; WT_Val];
                end
            end
        else 
            if edge_type == PASS
                Key_W = min(Max_Ref_Freq(Max_Ref_Freq>W_O));
                if WT_Val < Key_W
                    Minima_H = [Minima_H; Word_Length_Ceil];
                    Min_Ref_Freq = [Min_Ref_Freq; WT_Val];
                end
            else
                Key_W = max(Min_Ref_Freq(Min_Ref_Freq<W_O));
            end
        end
    end
    H_Max_Refi_Freq = cos(Max_Ref_Freq*Val)*best_Cosine_Coe - T*best_Cosine_Coe(1);
    H_Min_Refi_Freq = cos(Min_Ref_Freq*Val)*best_Cosine_Coe - T*best_Cosine_Coe(1);
    V_1   = H_Max_Refi_Freq > Upper_Boundary(Maxima_H)-100*SN;
    V_2   = H_Min_Refi_Freq < Lower_Boundary(Minima_H)+100*SN;
    Maxima_H = Maxima_H(V_1); Minima_H = Minima_H(V_2);
    Max_Ref_Freq = Max_Ref_Freq(V_1); Min_Ref_Freq = Min_Ref_Freq(V_2);
    H_Max_Refi_Freq = H_Max_Refi_Freq(V_1); H_Min_Refi_Freq = H_Min_Refi_Freq(V_2);
    Length_Maxima   = length(Maxima_H);
    Length_Minim   = length(Minima_H);
    E_FIl  = max([H_Max_Refi_Freq-Upper_Boundary(Maxima_H); Lower_Boundary(Minima_H)-H_Min_Refi_Freq; 0]);
    if E_FIl < SN
        break
    end
    if P_Parity == 1
        Lagrange_Constant = [ones(Length_Maxima,median+1); -ones(Length_Minim,median+1)].*cos([Max_Ref_Freq; Min_Ref_Freq]*(0:median));
        Lagrange_Constant(:,1) = Lagrange_Constant(:,1)/sqrt(2);
        end
    Lagrange_Matrix = [Upper_Boundary(Maxima_H); -Lower_Boundary(Minima_H)];

    if Filter_Weights
        Lagrange_Multiplier = (Lagrange_Constant*Ri*Lagrange_Constant')\(Lagrange_Constant*Ri*Coe-Lagrange_Matrix);
    else
        Lagrange_Multiplier = (Lagrange_Constant*Lagrange_Constant')\(Lagrange_Constant*Coe-Lagrange_Matrix);
    end
    [Neg_Mu,K] = min(Lagrange_Multiplier);
    while Neg_Mu < 0
        Lagrange_Constant(K,:) = [];
        Lagrange_Matrix(K) = [];
        if Filter_Weights
            Lagrange_Multiplier = (Lagrange_Constant*Ri*Lagrange_Constant')\(Lagrange_Constant*Ri*Coe-Lagrange_Matrix);
        else
            Lagrange_Multiplier = (Lagrange_Constant*Lagrange_Constant')\(Lagrange_Constant*Coe-Lagrange_Matrix);
        end
        [Neg_Mu,K] = min(Lagrange_Multiplier);
    end

    if Filter_Weights
       best_Cosine_Coe = Ri*(Coe-Lagrange_Constant'*Lagrange_Multiplier);
    else
        best_Cosine_Coe = Coe-Lagrange_Constant'*Lagrange_Multiplier;
    end

end

if P_Parity == 1
    H_Val = [best_Cosine_Coe(median+1:-1:2)/2; best_Cosine_Coe(1)/sqrt(2); best_Cosine_Coe(2:median+1)/2]';
else
    H_Val = [best_Cosine_Coe(median:-1:1); best_Cosine_Coe]'/2;
end

if nargout > 1
   best_Cosine_Coe = 1;
end
