%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% findFreeNeg.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [point,par,info]=findFreeNeg(point,par,tune,info)
% update local and determine activity
%
% for details of input and output structures see lmbopt.m

function [point,par,info]=findFreeNeg(point,par,tune,info)

% check for new activities
% identify the set I={i| low(i)<x(i)<upp(i)} (in paper:I_{-})

if par.sub
   point.actl=find(point.x==point.low); 
   point.actu=find(point.x==point.upp);    
end
point.Ineg=setdiff(1:point.n,union(point.actl,point.actu));

point.nIneg=length(point.Ineg);
point.gIneg=point.g(point.Ineg); 
% normgIneg = norm(gIneg,inf);
point.ngIneg2 = norm(point.gIneg,2)^2;
fixed = (point.nIneg<point.nI);


% decide localness and freeing
if fixed,
    % activity changed
    par.nlocal=0;
    if info.prt>=1, 
       disp(' ')
       disp(['findFreeNeg: activity has been changed',...
              ' (nlocal=',num2str(par.nlocal),')']);
       disp(' ')   
    end
elseif par.nstuck>0 
    % avoid cycling
    if par.nlocal>tune.nwait+min(point.m,point.nh), 
      par.nlocal=tune.nwait; 
       if info.prt>=1, 
        disp(' ')   
        disp(['findFreeNeg: nlocal (=',num2str(par.nlocal),...
             ') has been restarted']); 
        disp(' ') 
      end
    else               
      if par.sub, par.nlocal=par.nlocal+1; end
       if info.prt>=1,
         disp(' ')  
         disp(['findFreeNeg: nlocal (=',num2str(par.nlocal),...
             ') has been updated']); 
         disp(' ')
       end
    end;
elseif ~par.quad    
    % restart      
    if par.nlocal>=tune.nwait, 
      par.nlocal=tune.nwait; 
       if info.prt>=1, 
          disp(['findFreeNeg: nlocal (=',num2str(par.nlocal),...
             ') has been restarted']); 
      end
    else        
      if par.sub, par.nlocal=par.nlocal+1; end 
      if info.prt>=1, 
         disp(' ') 
         disp(['findFreeNeg: nlocal (=',num2str(par.nlocal),...
             ') has been updated']); 
         disp(' ')
       end
    end;
elseif ~fixed, % local
      if par.sub, par.nlocal=par.nlocal+1; end
      if info.prt>=1, 
         disp(' ') 
         disp(['findFreeNeg: nlocal (=',num2str(par.nlocal),...
             ') has been updated']); 
         disp(' ')
      end
end;
par.freeBounds = ...
    (~par.monotone | point.nIneg>point.nI | par.nlocal>=tune.nlf);

 if info.prt>=1,
    disp(' ')
    % display some other information
    if ~par.monotone,
        disp('findFreeNeg: no descent');
    elseif point.nIneg==0,
        disp('findFreeNeg: corner')
    elseif fixed,
        disp('findFreeNeg: fixed')
    end;
    disp(' ')
 end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
