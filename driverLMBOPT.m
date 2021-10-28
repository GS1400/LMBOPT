%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% driverLMBOPT.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file illustrates the use of LMBOPT
% by showing how to minimize the function f(x):=||Ax-b||_p^e for p=2,e=1
% in dimension n=100 from a random starting point.
% Here A =diag(logspace(5,0,m))*Q' with [Q,~]= qr(randn(n,m),0)
% The goal assumed in the example is the reduction of the initial
% objective function value by a factor of 0.01 with at most nfmax=200*n
% function evaluations, to return when one of the two limiting
% conditions is first satisfied.
clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create problem definition

rng(0);
m = 50; n = 100;
[Q,~] = qr(randn(n,m),0);
A  = diag(logspace(5,0,m))*Q';


% exactly 10^5,so Hessian has condition number 10^10
b = zeros(m,1);
% make 0 vector feasible (so we're sure this is a feasible problem) 
xLS = A\b; % e.g., xLS=0. This is not unique,

% as all constraints are inactive. 

low = -1e-3*n*ones(n,1); % for n=100, this is -0.1
upp = -low; % for n=100, this is +0.1 

% Define function values and gradient by subprograms
% provided in analogy to the following

prob.A = A; prob.b = b;
fun    = @(x)getfg(x,prob);

% To solve your own problem, simply replace in these lines
% the expression after @(x) by your expression or function call.
% Parameters in this expression other than x are passed by value,
% hence are fixed during minimization.

% problem definition complete
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% starting point
x = rand(n,1);
% start and stop info
nf2gmax = 1e6;
secmax  = inf;
% stop after secmax seconds 
epsilon = 1e-3;
% absolute accuracy requested for reduced gradient
% customization
prt=0; % print level
       % -1: nothing, 0: little, >=1: more and more 
       
Tuning = 1;   % 0: no tuning of LMBOPT (recommended at first)
              % 1: full tuning (for specialists only) 
fullinfo = 1; % 0: pass stop and print criteria inside LMBOPT
              % 1: pass stop and print criteria
if fullinfo    % define stop and print criteria
               % (indefinite run if no stopping criterion is given)
    st.secmax = secmax; % (default: inf)
    st.nf2gmax = nf2gmax; 
    st.prt = prt; % (default: -1)
    st.epsilon = epsilon;
    if prt>=0,st.time2print = cputime+1; end
else
    st = []; % budgets are chosen inside LMBOPT
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  solve the problem with LMBOPT
if Tuning % self-tuning and info
    % given are the defaults
    % only the deviating parameters need to be set!
    % iteration limit in efficient line search 
    tune.lmax = 3;
    % threshold for determining efficiency (0<beta<=0.2)
    tune.beta = 0.02;
    % threshold for efficiency of CG (0<betaCG<0.25)
    tune.betaCG = 0.001;
    % subspace dimension 
    tune.m = 12;
    % choose update formula (0 or 1)
    tune.typeH = 0;
    
else
   tune = []; % full tuning inside LMBOPT is used
end
if prt>=0
   disp(['cond(A)=',num2str(cond(A))]);
end
% call LMBOPT 
[x,f,info] = LMBOPT(fun,x,low,upp,tune,st);
if ~isempty(info.error),
   error = info.error
else
    if prt>=0 & isempty(info.error)
     % display output
     disp('============================================================')
     if n<=100, disp('best point found: xbest=');x',end
     % best point found
     disp('============================================================')
     disp(['function value at xbest: f(xbest)=', num2str(f)]);
     disp(['nf+2*ng:', num2str(info.nf2g)]); 
     secused =cputime-info.initTime;
     disp(['time used up to now: ', num2str(secused)]);
     disp(['reduced gradient at xbest: grad(xbest)=', ...
     num2str(info.acc)]);
     disp('============================================================')
     info
     % progress report by LMBOPT
     disp('============================================================')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
