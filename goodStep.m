
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% goodstep.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [step]=goodStep(point,step,tune);
%
% robust initial step for line search

function [step]=goodStep(point,step,tune);

                           
deltaa=tune.Deltaalpha;
exact = tune.exact;
q=tune.q;

if step.gTp>=0,
  step.p,step.gTp
  error('no descent direction')
end;


% first breakpoint
ind=(step.p<0 & point.x>point.low);
alow=min((point.low(ind)-point.x(ind))./step.p(ind));
if isempty(alow), alow=inf; end;
ind=(step.p>0 & point.x<point.upp);
aupp=min((point.upp(ind)-point.x(ind))./step.p(ind));
if isempty(aupp), aupp=inf; end;
step.abreak=min(alow,aupp);
step.abreak=step.abreak*(1+tune.Deltab) ;% allow for roundoff


% define minimal step size
if any(point.x==0 & step.p~=0),
   ratio = abs(point.f/step.gTp);
  %if ratio>=1
    step.amin=deltaa*ratio;
 % else
 %   step.amin=ratio;  
 % end
else 
  ind=(step.p~=0);
  if any(ind)
    ratio1 = abs(point.f/step.gTp);
    ratio2 = abs(point.x(ind)./step.p(ind));
   % if ratio1>=1 | ratio2>=1
       step.amin=deltaa*max(ratio1,min(ratio2));
   % else
   %    step.amin=max(ratio1,min(ratio2)); 
   % end
  else  
    % zero direction
    step.amin=1;step.step.agood=1;
    return;
  end;
end;

step.amin = min(step.amin,1);



% define target step size
atarget=max(step.amin,point.df/abs(step.gTp));
% ensure a free first step
if exact, atarget=min(atarget,step.abreak); end;
if atarget==step.abreak, 
else                
end;

% initial step size
if atarget<=step.abreak/q, 
  step.agood=atarget;
else                  
  step.agood=max(step.amin,step.abreak); 
  if step.agood==step.amin, end;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   



