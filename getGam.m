
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% getGam.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [point,step,par,info]=getGam(fun,point,step,par,tune,info)
% compute gamma

function [point,step,par,info]=getGam(fun,point,step,par,tune,info)
                     

point.df =   point.deltaf;                 

[step]=goodStep(point,step,tune);

if par.nstuck>=tune.nsmin, 
  % increase trial step to avoid remianing stuck
  step.agood=2*par.nstuck*step.agood;
end;


point.xnew = max(point.low,min(point.x+step.agood*step.p,point.upp));

% evaluate function
[point.fnew,point.res] = fun(point.xnew);
if isnan(point.fnew), point.fnew=inf; end; 
info.nf=info.nf+1;
info.nf2g=info.nf+2*info.ng;

% check stopping test
sec        = (cputime-info.initTime);
info.done  = (sec>=info.secmax)|(info.nf2g>=info.nf2gmax);
info.sec   = sec;

d2f=abs(point.fnew-point.f-step.agood*step.gTp)+eps;


par.gamma = 2*d2f/step.agood^2;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   