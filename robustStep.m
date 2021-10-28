%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% robustStep %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [point,step]=robustStep(point,step,tune)
% find a robust step size
%
% for details of input and output structures see lmbopt.m


function [point,step]=robustStep(point,step,tune)
 

% return robust step 

[dbest,irob]=min(point.dlist);
if dbest<0 % improvement
  step.anew=step.alist(irob);
  point.fnew=point.f+point.dlist(irob);
  point.res=point.rlist(:,irob);
  return;
else
    % treat failed line search (no improvement)
    % find point with smallest robust change
    % check for robust change
    ind=find(point.dlist>0 & isfinite(point.dlist));
    [dnew,inew]=min(point.dlist(ind)); % empty if ind is empty
    irob=ind(inew);
    if dnew<=point.df
      % accept point with robust change (never the case if dnew is empty)
    elseif isempty(dnew) | any(point.dlist<=0),
      % function almost flat; take step with largest d
      ind=find(isfinite(point.dlist));
      [dnew,inew]=max(point.dlist(ind));
      irob=ind(inew);
      if dnew>0 % never the case if dnew is empty
      else % function is flat; take largest step
        [~,inew]=max(step.alist(ind));
        irob=ind(inew);
      end
    else
      if dnew<=tune.Deltar*point.df
        % accept point with nonrobust change (at most twice too large)
        % this factor could be changed; >50 is unstable
      else % take largest almost flat step 
        ind=find(point.dlist<=point.df);
        [~,inew]=max(step.alist(ind));
        irob=ind(inew);
      end
    end
    step.anew=step.alist(irob); point.fnew=point.f+point.dlist(irob);
    point.res=point.rlist(:,irob);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   



 
