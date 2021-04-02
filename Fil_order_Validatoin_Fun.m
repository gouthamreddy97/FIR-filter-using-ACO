function [Discre_Word_Length,M1,M2,M_Obj] = Fil_order_Validatoin_Fun(Discre_Word_Length,Fend,Type_Filter,exception)
error(nargchk(3,4,nargin,'struct'));

if nargin == 3,
    exception = false;
end

M1 = '';
M2 = '';
M_Obj = [];
Ord_odd = false; 

if isempty(Discre_Word_Length) || length(Discre_Word_Length) > 1 || ~isnumeric(Discre_Word_Length) || ~isreal(Discre_Word_Length) || Discre_Word_Length~=round(Discre_Word_Length) || Discre_Word_Length<=0,
        return
end

if rem(Discre_Word_Length,2) == 1,
    Ord_odd = true; 
end
 
if (Type_Filter(end) ~= 0) && Fend == 1 && Ord_odd && ~exception,
    M_Obj = message('signal:firchk:NeedZeroGain');
    M2 = getString(M_Obj);
    Discre_Word_Length = Discre_Word_Length+1;
end
    

