%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% quasiNewtonDir.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [point,step,par]=quasiNewtonDir(point,step,par,tune)
% comput quasi Newton direction
%
% for details of input and output structures see lmbopt.m

function [point,step,par]=quasiNewtonDir(point,step,par,tune)

if point.nh>=point.m
    J=[1 point.m];  
else
    J=[point.m-point.nh+1 point.m]; 
end

Lam  = sqrt(sum(point.Y(point.I,J).^2'))./...
    sqrt(sum(point.S(point.I,J).^2'));
 
Ok   = (isnan(Lam) | Lam==0 | Lam>=tune.deltaD |Lam<=1/tune.deltaD);
        
Lam(Ok) = 1;  Lam = Lam(:); 


SS=point.S(point.I,:);
on = ones(size(SS,2),1);
U = point.Y(point.I,:)- Lam(:,on).*SS;
M = (point.Y(point.I,:)'*(point.Y(point.I,:)./Lam(:,on)))-point.H;
warning off
z     = M\(U'*(point.gI./Lam));
warning on 


step.pI = (U*z-point.gI)./Lam;
    
if any(isnan(step.pI))|any(isinf(step.pI))|all(step.pI==0)
   step.lamBB  = max(Lam);
   step.pI     = -step.lamBB*point.gI;
   step.gTp    = step.pI'*point.gI;   
end
 

point.Lam=Lam;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   