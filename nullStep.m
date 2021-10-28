
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% nullStep.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [point,step,par,info]=nullStep(fun,point,step,par,tune,info);
% try to prevent producing the null step
%
% for details of input and output structures see lmbopt.m

function [point,step,par,info]=nullStep(fun,point,step,par,tune,info)
%
norms = norm(step.s);
par.flags = (norms ==0);
if par.nnull>2&par.flags
    if info.eff == 4
        x1 = max(point.low,min(point.x0*(1-tune.del),point.upp));
    else
        x1 = max(point.low,min(point.x*(1-tune.del),point.upp));
    end
    ind       = (x1==0);
    x1(ind)   = tune.del;
    step.s    = x1-point.x;
    point.x   = x1;
    norms     = norm(step.s);
    par.flags = (norms==0);
    % evaluate function
    [point.fnew,point.res] = fun(point.x);
    if isnan(point.fnew), point.fnew=inf; end; 
    info.nf=info.nf+1;
    info.nf2g=info.nf+2*info.ng;
    
    if info.prt>=1
         disp('')
        if par.flags 
            disp('nullStep: the null step has not been removed')
        else
            disp('nullStep: the null step has been removed') 
        end
         disp('')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
