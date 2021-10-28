%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% AvoidZigzagDir.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [step]=AvoidZigzagDir(point,step,par,tune)
%
% comput a direction to avoid ziggzaging
%
% for details of input and output structures see lmbopt.m
%


function [step]=AvoidZigzagDir(point,step,par,tune)

step.pI = step.p(point.I);

Bet = abs(point.gI./step.pI);

I = find(~isinf(Bet)&~isnan(Bet));
if ~isempty(I)
    Bet = Bet(I);
    betak = tune.theta*max(Bet);
else
    betak = tune.theta;
end
lam        = (1+betak*step.gTp)/point.omega;
step.pI    = betak*step.pI -lam*point.gI;

if any(isnan(step.pI))|any(isinf(step.pI))|all(step.pI==0)
    if ~isnan(step.lamBB)
       if par.nlocal==tune.nwait & point.nh>=1
          step.lamBB  = max(point.Lam);
       end
       step.pI     = -step.lamBB*point.gI;
       step.gTp    = step.pI'*point.gI;
    else
       step.pI = -point.gI; step.gTp=point.gI'*step.pI;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   