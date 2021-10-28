%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% getSuccess.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [point,step,par,info]=getSuccess(fun,point,step,par,tune,info)
% determine whether iteration is successful or not
%
% for details of input and output structures see lmbopt.m
%

function [point,step,par,info]=getSuccess(fun,point,step,par,tune,info)

par.quad=0;
estuck=0;
ok = par.cosine<0 & any(step.p~=0) & (par.CG>0 | par.nstuck>=tune.nsmin);
if ok
    % data are consistent with a definite quadratic function
    % evaluate function
    [point.fnew,point.res] = fun(point.xnew);
    if isnan(point.fnew), point.fnew=inf; end; 
    info.nf=info.nf+1;
    
    info.nf2g=info.nf+2*info.ng;
    
    % check stopping test
    sec        = (cputime-info.initTime);
    info.done  = (sec>=info.secmax)|(info.nf2g>=info.nf2gmax);
    info.sec   = sec;
  
    
      
    mu=(point.fnew-point.f)/step.gTp;
    if point.fnew==point.f | mu*abs(mu-1)>=tune.betaCG, 
      par.quad=par.CG;
    end;
    % here fnew=f is accepted since the function may be already flat
    % although a nontrivial step was taken
    % note that mu=0 may hold if fnew>f due to underflow!

    estuck=( par.nstuck>=tune.nsmin & point.fnew<=point.f+point.deltaf);
    % here a robust increase is counted as success if stuck enough
end;

par.success=(par.quad | estuck);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   