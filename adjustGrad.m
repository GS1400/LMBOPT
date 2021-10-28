

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% adjustg.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function point=adjustGrad(point,tune);
% adjust meaningless gradient components (NaN,inf,-inf)
%

function point=adjustGrad(point,tune)


% the three disp lines should be deleted later
ind=isnan(point.g);
if any(ind), 
  % force freeability
  ind2=(point.x-point.low>point.upp-point.x);
  point.g(ind&ind2)=tune.Deltag;
  point.g(ind&~ind2)=-tune.Deltag;
end;
ind=(point.g==inf);
if any(ind), 
  point.g(ind)=tune.Deltag; 
end;
ind=(point.g==-inf);
if any(ind), 
  point.g(ind)=-tune.Deltag; 
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   




