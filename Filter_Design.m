function varargout = Filter_Design(Fil_Design_Weights, type, varargin)
Filter_Obj_Arg = {};
Filter_Obj_Arg_Design = {};
type = lower(type);
switch type
       case {'fir', 'iir'}
        Design_Types = designmethods(Fil_Design_Weights, type,Filter_Obj_Arg{:});
        if strcmpi(type, 'fir') && any(strcmpi(Design_Types, 'equiripple'))
            Design_Types = 'equiripple';
        end
        Filter_Evaluate_Quo = feval(Design_Types, Fil_Design_Weights,Filter_Obj_Arg{:});
        varargout = {Filter_Evaluate_Quo};
    otherwise
        if any(strcmpi(type, designmethods(Fil_Design_Weights))) || ...
                any(strcmpi(type, hiddenmethods(Fil_Design_Weights)))
            Filter_Evaluate_Quo = feval(type, Fil_Design_Weights, varargin{:},Filter_Obj_Arg_Design{:});
            varargout = {Filter_Evaluate_Quo};
        end
end
end