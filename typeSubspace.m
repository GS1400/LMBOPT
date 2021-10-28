%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% typeSubspace.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [point,par]=typeSubspace(point,par,tune,info)
% determine the type of subspace
%
% for details of input and output structures see lmbopt.m
%

function [point,par]=typeSubspace(point,par,tune,info)

  mhat = min(point.m,point.nh);
  if par.nlocal<tune.nwait,
    % ordinary subspace step: full subspace direction 
    % may be contaminated by now active components
    % may lead to premature freeing if used directly
    point.m0 = min(info.ng-1,mhat); par.hist =1:point.m0; par.CG = 0;
     if info.prt>=1, 
        disp(' ')
        disp('typeSubspace.m: the ordinary subspace step' ); 
        disp(' ')
    end
  elseif par.nlocal==tune.nwait
    % restart: quasi newton direction
    point.m0   = 0;
    par.hist = [];
    par.CG   = 1;
    
    % permute subspace basis so that oldest column is first
   

     perm        = [point.ch+1:point.m,1:point.ch];
     point.S     = point.S(:,perm);
     point.Y     = point.Y(:,perm);
     point.H     = point.H(perm,perm);
     point.ch    = 0;

    
    if info.prt>=1, 
        disp(' ')
        disp('typeSubspace.m: the subspace basis has been permuted ');
        disp(' ')
    end
    
  elseif par.nlocal<tune.nwait+mhat,
    % preserve conjugacy by restricting the subspace 
    point.m0    = par.nlocal-tune.nwait;
    par.hist  = 1:point.m0;
    par.CG    = 2;
    if info.prt>=1, 
       disp(' ')
       disp(['typeSubspace.m: the conjugacy has been preserved by', ...
              ' restricting the subspace']); 
       disp(' ')
    end
  else 
    % full subspace step preserves conjugacy
    point.m0   = mhat;
    par.hist = 1:mhat;
    par.CG   = 3;
    if info.prt>=1, 
       disp(' ')
       disp(['typeSubspace.m: the conjugacy has been preserved',...
             ' by the full subspace step']); 
       disp(' ')
    end
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   