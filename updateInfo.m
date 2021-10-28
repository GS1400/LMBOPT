
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% updateInfo.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [point,par]=updateInfo(point,par,tune,info)
% update deltaf and the best point
%
% for details of input and output structures see lmbopt.m
%
function [point,par]=updateInfo(point,par,tune,info)

% update fbest
if point.fnew<point.fbest, 
    par.nstuck=0; point.fbest=point.fnew;
else % possible cycle
     % but often the gradient still decreases nonmonotonically
     par.nstuck=par.nstuck+1;
end;  

% update f
dec = point.fnew<point.f;
if dec % improvement
    par.monotone=1;
    point.deltaf=point.f-point.fnew;flagmod = mod(info.ng, tune.mf);
    if (flagmod == 0),point.Df(tune.mf) = point.deltaf;
    else, point.Df(flagmod) = point.deltaf;end
    elseif point.fnew==point.f 
    % stalled
    par.monotone=0;
else % no descent
    par.monotone=0;
    point.deltaf = max(tune.Deltaf*point.deltaf,...
                   tune.Deltam*(abs(point.f)+abs(point.fnew)));
    flagmod = mod(info.ng, tune.mf);
    if (flagmod == 0),point.Df(tune.mf) = point.deltaf;
    else, point.Df(flagmod) = point.deltaf;
    end
end;


point.f=point.fnew;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
