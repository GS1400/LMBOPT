%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% InitialInfo.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [point]=InitialInfo(point,tune)
% Initialize best point and \Deltaf
%
% for details of input and output structures see lmbopt.m
%

function [point]=InitialInfo(point,tune)

if isnan(point.f), point.f=inf; end;
point.fbest=point.f;
if isfinite(point.f) & point.f~=0, 
    point.deltaf=tune.facf*abs(point.f);
else,  point.deltaf=1;
end;

point.Df(1:tune.mf)  = -Inf; point.Df(tune.mf)    = point.deltaf;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   