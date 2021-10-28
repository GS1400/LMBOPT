%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% findFreePos.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [point,par,info]=findFreePos(point,par,info)
% update local and determine activity
%
% for details of input and output structures see lmbopt.m

function [point,par,info]=findFreePos(point,par,info)
                
rho        = (1/max(1,info.ng-1));
par.freeBounds = (par.freeBounds | point.ngIneg2<rho*point.ngred2);
% find appropriate free variables
if par.freeBounds
    % identify I = {i | l(i)<x(i)<u(i) or gred(i)=~0} 
    
    point.actl=find(point.x==point.low); 
    point.actu=find(point.x==point.upp);
    point.Ipos=union(find(point.gred),...
                setdiff(1:point.n,union(point.actl,point.actu)));
    
    point.nIpos=length(point.Ipos);
    if info.ng==1 | point.nIpos>point.nIneg,
      % freeing step 
      point.Ineg=point.Ipos; par.nlocal=0;
      if info.prt>=1, 
          disp(' '); disp('findFreePos: freeing step'); disp(' ');
      end
    end;
end; % if freeBounds

% update free index set
point.I     = point.Ineg;
point.nI    = length(point.I); 
info.nI     = point.nI; 
point.gI    = point.g(point.I); 
point.omega = point.gI'*point.gI; 
  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   