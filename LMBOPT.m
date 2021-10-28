%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LMBOPT.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x,f,info]=LMBOPT(fun,x,low,upp,tune,st)
% limited memory method for unconstrained and bound-constrained 
% optimization
%
% solves the smooth optimization problem
%     min  f(x)
%     s.t. low<=x<=upp
% lower and upper bounds may be infinite.
%
% fun      % function and gradient handle for f(.)
% x        % starting point (must be specified)
% st       % structure with stop and print criteria
%          % (indefinite run if no stopping criterion is given)
%  .secmax       %   stop if sec>=secmax     (default: inf)
%  .nfmax        %   stop if nf>=nfmax       (default: inf)
%  .ngmax        %   stop if ng>=ngmax       (default: inf)
%  .nf2gmax      %   stop if nf2g>=nf2gmax   (default: inf)
%  .prt          %   printlevel (default: -1)
%                %   -1: nothing, 0: little, >=1: more and more
% tune     % optional structure containing tuning parameters
%          %   for details see below
%
% x        % best point found 
% f        % function value at best point found 
% info     % structure with performance info
%          %   for details see below
% 
% 
function [x,f,info]=LMBOPT(fun,x,low,upp,tune,st);

 
% measure to get performace information 
% even when program is interrupted by funcall

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% check inputs %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% feasible region

info.error='';

n=length(low);
if length(upp)~=n, 
  nlow=n,nupp=length(upp),
  info.error='dimension mismatch'; 
end;
if n==0,
  % no variables
  x=zeros(0);f=inf;
  info.error='no variables';
  return;
end;
if ~all(low<=upp), 
  % feasible domain empty
  x=NaN;f=NaN;
  info.error='feasible domain empty';
  return;
end;

if all(low<=-1e20) & all(upp>=1e20)
   par.probcase = 0; % unconstrained case
else
   par.probcase = 1; % bound constrained case
end


% default parameters
resfile='results.txt';  % result file

% initialize structure containing all tuning parameters 
%
% tune % structure containing all tuning parameters 
%      % all parameters have a default that can be overwritten 
%      % by specifying it as an input
%  .nsmin      % how many stucks before taking special action?, >0
%  .nwait      % number of local steps before CG is started
%  .rfac       % restart after rfac*I local steps
%  .nlf        % number of local steps before freeing is allowed
%  .Deltam     % tiny factor for regularizing Deltaf if not monotone 
%  .Deltapg    % tiny factor for regularizing gTp
%  .DeltaH     % tiny regularization factor for subspace Hessian, <<1
%  .Deltaalpha % tiny factor for initial step
%  .lmax       % iteration limit in efficient line search
%  .beta       % threshold for determining efficiency (0<beta<=0.2)
%  .betaCG     % threshold for efficiency of CG (0<betaCG<0.25)
%  .scale      % initial scaling step?
%  .Deltaangle % regularization angle
%  .Deltareg   % factor for regularization
%  .Deltaw     % for generating w>0 in enforceAngle.m
%  .facf       % relative accuracy of f in first step
%  .Deltax     % tiny factor for interior move, <<1
%  .m          % subspace dimension
%  .mf         % length Df
%  .epsilon    % absolute accuracy requested for reduced gradient
%  .typeH      % choose update formula (0 or 1)
%  .nnullmax   % iteration limit in null steps
%  .del        % factor for null step
%  .Deltar     % factor for finding almost flat step
%  .Deltag     % factor for adjusting gradient
%  .Deltab     % tiny factor for breakpoint
%  .Deltau     % factor for upp
%  .theta      % parameter for adjusting direction
%  .stuckmax   % iteration limit in number of stuck
%  .zetamin    % lower bound for zeta
%  .zetamax    % upper bound for zeta
%  .deltaD     % parameter to adjust the scaling matrix D
%  .q          % parameter for extrapolation in CLS
%  .qmin       % lower bound for q when agood >= 1
%  .Deltaq     % parameter for reducing q when agood >= 1
%  .Deltaf     % parameter for expanding deltaf


if ~exist('tune'), tune=[]; end;

% how many stucks before taking special action?, >0
if ~isfield(tune,'nsmin'), tune.nsmin = 1; end;
% number of local steps before CG is started
if ~isfield(tune,'nwait'), tune.nwait = 1; end;
% restart after rfac*I local steps
if ~isfield(tune,'rfac'), tune.rfac = 2.5; end;
% number of local steps before freeing is allowed
if ~isfield(tune,'nlf'), tune.nlf = 2; end;
% tiny factor for regularizing Deltaf if not monotone 
if ~isfield(tune,'Deltam'), tune.Deltam = 1e-13; end;
% tiny factor for regularizing gTp
if ~isfield(tune,'Deltapg'), tune.Deltapg = eps; end;
% tiny regularization factor for subspace Hessian, <<1
if ~isfield(tune,'DeltaH'), tune.DeltaH = eps; end;
% tiny factor for initial step
if ~isfield(tune,'Deltaalpha'), tune.Deltaalpha = 5*eps; end;
% iteration limit in efficient line search
if ~isfield(tune,'lmax'), tune.lmax = 4; end;
% threshold for determining efficiency (0<beta<=0.2)
if ~isfield(tune,'beta'), tune.beta = 0.02; end;
% threshold for efficiency of CG (0<betaCG<0.25)
if ~isfield(tune,'betaCG'), tune.betaCG = 0.001; end;
% regularization angle
if ~isfield(tune,'Deltaangle'), tune.Deltaangle = 1e-12; end;
% factor for regularization
if ~isfield(tune,'Deltareg'), tune.Deltareg = 1e-12; end;
% for generating w>0 in enforceAngle.m
if ~isfield(tune,'Deltaw'), tune.Deltaw = eps; end;
% relative accuracy of f in first step
if ~isfield(tune,'facf'), tune.facf = 1e-8; end;
% tiny factor for interior move, <<1
if ~isfield(tune,'Deltax'), tune.Deltax = 1e-20; end;
% subspace dimension
if ~isfield(tune,'m'), tune.m = 12; end;
% length Df
if ~isfield(tune,'mf'), tune.mf = 2; end;
% parameter for Df
if ~isfield(tune,'mdf'), tune.mdf = 20; end;

% choose update formula (0 or 1)
if ~isfield(tune,'typeH'), tune.typeH = 0; end;
% iteration limit in null steps
if ~isfield(tune,'nnullmax'), tune.nnullmax = 3; end;
% factor for null step
if ~isfield(tune,'del'), tune.del = 1e-10; end;
% factor for finding almost flat step
if ~isfield(tune,'Deltar'), tune.Deltar = 20; end;
% factor for adjusting gradient
if ~isfield(tune,'Deltag'), tune.Deltag = 100; end;
% tiny factor for breakpoint
if ~isfield(tune,'Deltab'), tune.Deltab = eps; end;
% factor for upp
if ~isfield(tune,'Deltau'), tune.Deltau = 1000; end;
% parameter for adjusting direction
if ~isfield(tune,'theta'), tune.theta = 1e-8; end;
% enforce exact line search on quadratic
if ~isfield(tune,'exact'), tune.exact = 0; end;
% gradient tolerance for skipping update (used in Powell condition)
if ~isfield(tune,'Deltapo'), tune.Deltapo = eps; end;
% iteration limit in number of stuck
if ~isfield(tune,'nstuckmax'), tune.nstuckmax = inf; end; 
% lower bound for zeta
if ~isfield(tune,'zetamin'), tune.zetamin = -1e+50; end; 
% upper bound for zeta
if ~isfield(tune,'zetamax'), tune.zetamax = -1e-50; end; 
% parameter to adjust the scaling matrix D
if ~isfield(tune,'deltaD'), tune.deltaD = 1e10; end; 
% parameter for extrapolation in CLS
if ~isfield(tune,'q'), tune.q = 25; end;
% lower bound for q when agood >= 1
if ~isfield(tune,'qmin'), tune.qmin = 2.5; end; 
% parameter for reducing q when agood >= 1
if ~isfield(tune,'Deltaq'), tune.Deltaq = 10; end; 
% parameter for expanding deltaf
if ~isfield(tune,'Deltaf'), tune.Deltaf = 2; end; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% initialize point %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

point.n    = n;   point.x    = x; 
point.low  = low; point.upp  = upp;
point.Lam = []; point.y= []; point.Ineg=NaN; 
point.m=max(1,min(tune.m,point.n));
point.nIneg=NaN; point.Y=zeros(point.n,point.m); 
point.S=zeros(point.n,point.m); point.nh = 0; 
point.H=zeros(point.m); point.ch=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% improve the initial point %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[point]=projStartPoint(point,tune);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% compute f and g at x %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[f,res] = fun(point.x); g = fun(res);
info.nf=1; info.ng=1;
point.f=f;point.g=g; point.res=res;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% initialize Info %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[point]=InitialInfo(point,tune);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% adjust meaningless gradient components %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

point=adjustGrad(point,tune);

point.ngIneg2 = norm(point.g,2)^2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% initialize step %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

step.p=[]; step.lamBB= NaN;
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% initialize par %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

par.nlocal=-1; par.freeBounds=1; par.nstuck=0; 
par.nnull=0; par.flags = 0; par.sub=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% initialize info %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

info.status = -1; 


if ~exist('st'), st=[]; end;


% print level
if isfield(st,'prt'), info.prt = st.prt; 
else,                 info.prt = 0; 
end
% stopping criteria
if isfield(st,'secmax'), info.secmax=st.secmax;
else, info.secmax=inf;
end
if isfield(st,'nfmax'), info.nfmax=st.nfmax;
else, info.nfmax=1000*n;
end

if isfield(st,'ngmax'), info.ngmax=st.ngmax;
else, info.ngmax=1000*n;
end

if isfield(st,'nf2gmax'), info.nf2gmax=st.nf2gmax;
else, info.nf2gmax=info.nfmax+2*info.ngmax;
end

if isfield(st,'time2print'), info.time2print=st.time2print;
else, info.time2print=cputime+1;
end



% absolute accuracy requested for reduced gradient
if isfield(st,'epsilon'), info.epsilon=st.epsilon;
else, info.epsilon = 1e-6; 
end;

prt=info.prt;

info.nf2g = 3; info.nI=NaN; 

info.initTime = cputime;
if prt>=1, itc = 0; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% append to current file %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid=fopen(resfile,'a');  
if fid<0, 
  resfile
  error('file not accessible - permissions?'); 
end;
fprintf(fid,'\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% main loop %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while 1
    
    if prt>=1
        itc = itc+1;
        disp('=============================================== ') 
        disp(['The ',num2str(itc),'th iteration: ']);
        disp('===============================================') 
    end
    
    

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%% get reduced gradient %%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  [point]=redGrad(point);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%% find I_+ %%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  [point,par,info]=findFreePos(point,par,info);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%% check stopping tests %%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  info.acc=norm(point.gred,Inf); 

  % check stopping test
  sec=(cputime-info.initTime);
  info.done=(sec>=info.secmax)|(info.nf2g>=info.nf2gmax);
  info.done=(info.done|info.acc<info.epsilon|...
             par.nstuck>=tune.nstuckmax);
  info.sec = sec;
  
  
  
  
  if info.prt==0 
      if cputime>=info.time2print
         info.time2print=cputime+1;
         disp('============================================')
         disp(['nf= ',num2str(info.nf),', ng= ',num2str(info.ng),...
             ', nf2g= ', num2str(info.nf2g),', time= ',num2str(sec)])
         disp(' ..........................................  ')
         disp(['||gred||_inf= ',num2str(info.acc)])
         disp(' .......................................... ')
         disp(['fbset= ',num2str(point.fbest),', f= ',num2str(point.f)]) 
         disp(' ..........................................  ')
         disp(['|Ineg|= ',num2str(point.nIneg),...
               ' , |Ipos|= ',num2str(point.nIpos)])
         disp(' ..........................................  ')
         disp(['The length of subspace, m0:', num2str(point.m0)])
         disp(' ..........................................  ')
         disp(['The number of stuck, nstuck:', num2str(par.nstuck)])
         disp(' ..........................................  ')
         disp(['The number of null step, nnull:', num2str(par.nnull)])
         disp(' ..........................................  ')
         disp('============================================')
      end;
  end

   if prt>=0
      disp(' ') 
      disp(['f= ',num2str(point.f),' improved at nf2g=',...
            num2str(info.nf2g)]) 
      disp(' ')   
  end    
  
  % check stopping tests 
  if info.done, break;  end;
  

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%% test for local restart %%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  if par.nlocal>max(tune.nwait,tune.rfac*point.n)
    % restart CG
    par.nlocal=tune.nwait; 
   
    if prt>=1
        disp(' ') 
        disp('local restart is done');
        disp(' ') 
    end
  end;
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%% determine type of subspace %%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  [point,par]=typeSubspace(point,par,tune,info);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%% construct search direction %%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  [point,step,par]=searchDir(point,step,par,tune,info);
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%% construct ConjGradDir direction %%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
                          
  [point,step,par,info]= ...
             ConjGradDir(fun,point,step,par,tune,info);
         
  if info.done 
      % nf2gmax or secmax is exceeded inside ConjGradDir
      point.x = point.xnew; point.f = point.fnew;
      break; 
  end
 
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%% scaled Cauchy point %%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
  
  if par.probcase % only for bound constrained case
     [point,step,par]=scaleCauchy(point,step,par,tune,info);
  else
      step.gTp=point.gI'*step.pI; 
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%% Is iteration successful ? %%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  [point,step,par,info]= ...
              getSuccess(fun,point,step,par,tune,info);
  
  if info.done
      
      % nf2gmax or secmax is exceeded inside getSuccess
      point.x = point.xnew; point.f = point.fnew;
      break; 
  end
  
  if par.success, % subspace step is considered successful
   % update best point 
    point.x = point.xnew; step.s  = step.p; par.flags=0;
     par.successold =1;
    if prt>=1
        disp(' ') 
        disp('getSuccess: the subspace step is successful');
        disp(' ') 
    end  
    
    
  else % function not close to quadratic
       % perform a line search along a regularized direction
    if prt>=1
        disp(' ') 
        disp('getSuccess: the subspace step is unsuccessful');
        disp(' ') 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%% enforce angle condition %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    par.firstAngle=0;
    [step]=enforceAngle(point,step,par,tune,info);
    
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% run Curve line search %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
     if ~par.probcase, par.sub=1; end
    
    [point,step,par,info]=CLS(fun,point,step,par,tune,info);
    
    if info.done
        % nf2gmax or secmax is exceeded inside CLS
        point.x = point.xnew; 
        point.f = point.fnew;
        break; 
    end
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% update best point %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % update best point
    point.x0 = point.x;
    step.s   = point.x; 
    point.x  = max(point.low,min(point.x+step.anew*...
               step.p,point.upp)); 
    step.s  = point.x-step.s;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% remove null steps %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [point,step,par,info]= ...
                 nullStep(fun,point,step,par,tune,info);
   
    if par.flags % null step; avoid cycling
       par.nnull=par.nnull+1;  
       if par.nnull>tune.nnullmax, break; end;
    end
    
    
    
  end; % of if success
  
    
  if ~par.flags
    % significant step; get new gradient
    par.nnull=0; g0=point.g; 
    point.g = fun(point.res);
    info.ng   = info.ng+1; 
    info.nf2g = info.nf2g+1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% adjust meaningless gradient components %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    point=adjustGrad(point,tune); 
    
    point.y=point.g-g0;
    
    
   end; 


   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%% update information %%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   [point,par]=updateInfo(point,par,tune,info);
  
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%% find I_-  %%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  [point,par,info]=findFreeNeg(point,par,tune,info);
 
  
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%% update subspace %%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
   [point]=updateSubspace(point,step,par,tune);
  

  
end; % of while 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% create output info %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% summarize line searches
% status detection

x=point.x; f=point.f; 

if info.acc<=info.epsilon,  info.status=0; 
elseif par.nstuck>=tune.nstuckmax, info.status=1;
elseif par.nnull>=tune.nnullmax, info.status=2;
elseif info.nf2g>=info.nf2gmax, info.status=3;
else, info.status=4;
end


fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

