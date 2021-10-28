%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% scaleDir.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [step,par]=scaleDir(point,step,par)
% choose components of sensible sign and scale
%
% for details of input and output structures see lmbopt.m
%

function [step,par]=scaleDir(point,step,par)

sc    = zeros(point.n,1);
Jz    = find(point.x==0);
if ~isempty(Jz)
    sc(Jz)       = min(1,point.upp(Jz)-point.low(Jz));
    Jzz          = find(point.g(Jz)<0);
    Je           = setdiff(Jz,Jzz);
    step.p(Jzz)  = sc(Jzz);
    step.p(Je)   = -sc(Je);
end

Jnz   = setdiff(1:point.n,Jz);
if ~isempty(Jnz)
     sc(Jnz)     = abs(point.x(Jnz)); 
     Jl          = find(point.x(Jnz)==point.low(Jnz));
     step.p(Jl)  = sc(Jl);
     Ju          = find(point.x(Jnz)==point.upp(Jnz));
     step.p(Ju)  = -sc(Ju);
     Jlu         = union(Jl,Ju);
     Jg          = setdiff(Jnz,Jlu);
     Jg          = find(point.g(Jg)<0);
     step.p(Jg)  = sc(Jg);
     Je1         = setdiff(Jnz,union(Jlu,Jg));
     step.p(Je1) = -sc(Je1);
end


step.p=step.p(:); step.pI = step.p(point.I); par.CG = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   