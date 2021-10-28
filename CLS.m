

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLS.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [point,step,par,info]=CLS(fun,point,step,par,tune,info)
% find a step size alpha with mu(alpha)|mu(alpha)-1| >= beta
%
% for details of input and output structures see lmbopt.m


function [point,step,par,info]=CLS(fun,point,step,par,tune,info)

% full search direction
step.p=zeros(point.n,1); 
step.p(point.I)=step.pI;

if info.ng==1, point.df = point.deltaf;
else
  absf = abs(point.f);
  if mod(info.ng,tune.mdf)==0, point.df=point.deltaf*(absf+1);
  else, point.df=max(point.Df);
  end
end


if isfield(tune,'exact'), exact=tune.exact;
else, exact=0;
end


% check consistency of input
if ~(step.gTp<0)
 step.gTp;
 error('gTp <0 is required')
end



% compute robust initial step
if par.sub | info.ng~=1
  % compute robust initial step
  [step]=goodStep(point,step,tune);
else
   step.agood=1; step.amin =1; 
end

if step.agood>=1,tune.q=max(tune.qmin,tune.q/tune.Deltaq);end

% use robust initial step
step.anew=step.agood;


% define bracket [alow,aupp] for step sizes 
alow=0;
aupp=inf;

% initialize lists
alist(1)=0;  % step sizes 
dlist(1)=0;  % differences in function values
info.eff=NaN;
descent=0;
rob=0;
itc=0;



while 1
    
  % evaluate function at current step
  [point.fnew,res]=fun(max(point.low,min(point.x+step.anew*...
             step.p,point.upp)));
         
  if isnan(point.fnew), point.fnew=inf; end;
   
  itc=itc+1;
  
  info.nf = info.nf+1; info.nf2g = info.nf2g+1;
  
  % check stopping test
  sec        = (cputime-info.initTime);
  info.done  = (sec>=info.secmax)|(info.nf2g>=info.nf2gmax);
  info.sec   = sec;
  
  if info.done, break; end
  
 
  dnew=point.fnew-point.f;
  alist(itc+1)=step.anew;
  dlist(itc+1)=dnew;
  rlist(:,itc+1)=res;
  
  % quit line search?
  mu=dnew/(step.anew*step.gTp);
  if mu*abs(mu-1)>=tune.beta | info.eff==1 % line search efficient
    info.eff=1;
    if info.prt>=1, 
        disp(' ') 
        disp('CLS: line search efficient'); 
        disp(' ') 
    end
    if ~exact, 
      break; 
    end;
  elseif itc>1 & ~descent & rob>0 % robust nonmonotone step accepted
    info.eff=2;
     if info.prt>=1, 
         disp(' ') 
         disp('CLS: robust nonmonotone step accepted'); 
         disp(' ') 
     end
    break;
  elseif itc>=tune.lmax % limit on function values reached
    if descent % descent
      info.eff=3;
      if info.prt>=1, 
          disp(' ') 
          disp(['CLS: limit on function values reached',...
                ' with improvement on f']);
          disp(' ')   
      end
      break;
    else % no descent
      info.eff=4;
      if info.prt>=1, 
        disp(' ')   
        disp(['CLS: limit on function values reached',...
            ' without improvement on f']);
        disp(' ') 
      end
      break;
    end;
  end;

  % update bracket
  if descent % update bracket for descent
    if mu>0.5
      alow=step.anew;
    else
      aupp=step.anew;
    end
  elseif dnew<0 % first descent step
    descent=1;
    % create bracket for descent
    % lower part
    ind=find(dlist<0);
    alow=max(alist(ind));
    rob=-1;
    % upper part
    ind=find(dlist>=0 & alist>alow);
    if isempty(ind), aupp=inf;
    else, aupp=min(alist(ind));
    end
  else % no descent; update robust bracket
    if dnew<=point.df
      alow=step.anew;
      rob=dnew;
    else
      aupp=step.anew;
    end
  end

  % determine next step
  if itc==1 % first step
    if mu<1 % secant step for Goldstein quotient
      step.anew=step.anew*0.5/(1-mu); 
      if step.anew==0,step.anew=step.amin; end;
    else % extrapolation
      step.anew=max(step.amin,step.anew*tune.q);
    end
    exact=0; 
  else % robust bisection 
    if aupp==inf % extrapolation
      step.anew=step.anew*tune.q;
    elseif alow==0 % contraction
      step.anew=step.anew/tune.q;
    else  % geometric mean bisection
       alow0 = max(alow,step.amin); 
       step.anew  = alow0*sqrt(aupp/alow0);
    end;
  end; % of if itc

end; % of while 

point.dlist = dlist;
point.rlist = rlist;
step.alist = alist;

% return robust step 
[point,step]=robustStep(point,step,tune);

point.dlist = [];
step.alist  = [];
point.rlist = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   



 
