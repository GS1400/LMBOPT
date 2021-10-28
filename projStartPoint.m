%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% projStartPoint.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [point]=projStartPoint(point,tune)
% project starting point (to the interior when Deltax>0)
%
% for details of input and output structures see lmbopt.m
%


function [point]=projStartPoint(point,tune)

if tune.Deltax==0,
    
  point.x=max(point.low,min(point.x,point.upp));
else
  ind=find(point.x<=point.low);
  point.x(ind)=point.low(ind)+tune.Deltax*min(tune.Deltau*...
              abs(point.low(ind)),point.upp(ind)-point.low(ind));
  ind=find(point.x>=point.upp);
  point.x(ind)=point.upp(ind)-tune.Deltax*min(tune.Deltau*...
              abs(point.upp(ind)),point.upp(ind)-point.low(ind));
end;
ind=find(~isfinite(point.x));
% force starting point finite
point.x(ind)=max(point.low(ind),min(0,point.upp(ind)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
