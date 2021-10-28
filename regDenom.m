
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% regDenom.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function par=regDenom(point,step,par,tune)
% construct regularize denomiator

function par=regDenom(point,step,par,tune)

par.denom=par.gamma-step.q'*step.r;
% regularize denominator
e2f=abs(point.fnew-point.f)+step.agood*(abs(point.g)'*abs(step.p));
dcor=tune.DeltaH*(e2f/step.agood^2+abs(step.q)'*abs(step.r));
if par.denom>=0, par.denom=par.denom+dcor; 
else, par.denom=par.denom-dcor; 
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   