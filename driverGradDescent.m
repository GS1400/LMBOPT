
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% driverGradDescent.m %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

rng(0);
m = 50;
n = 100;
[Q,~] = qr(randn(n,m),0);
A = diag(logspace(5,0,m))*Q';
cond(A) % exactly 10^5, so Hessian has condition number 10^10
b = zeros(m,1); % make 0 vector feasible (so we're sure this is a feasible problem)
xLS = A\b; % e.g., xLS = 0. This is not unique, as all constraints are inactive.
low = -1e-3*n*ones(n,1); % for n=100, this is -0.1
upp = 1e-3*n*ones(n,1); % for n=100, this is +0.1

prob.A=A; prob.b=b;
% getwhat % 1: get f, 2: get g, 3: get f and g
fun  = @(x) getfg(x,prob); 
x0 = ones(n,1); % all algorithms will use the same starting point


nf2gmax = inf;
secmax  = 60;

cas = 3; % cas: 1  % Gradient Descent
         %      2  % FISTA
         %      3  % LMBOPT

switch cas
    case 1 % Solve via Gradient Descent
        
    getwhat=3;  % compute f and g
    
    t = 1/norm(A)^2;
    
    x = x0;
    y = x;
    if ~isinf(nf2gmax)
       hist = zeros(nf2gmax,1);
    end
    proj = @(x) max( min( x , upp ), low );
    nf=0; ng=0;
    tic;
    initTime=cputime;
    while 1
        xOld= x;
        [fx,res] = fun(x);
        gx       = fun(res);
        x = proj( x - t*gx );
        hist(nf+1) = fx;
        nf=nf+1; ng=ng+1;
        nf2g = nf+2*ng;
        sec=(cputime-initTime);
        if (nf2g>=nf2gmax | sec>=secmax)
            disp('Basic gradient descent ended')
            break;
        end
        disp(['f= ',num2str(fx),' improved at nf2g=',num2str(nf2g)])
    end
    
    toc

    case 2 % Solve via FISTA.
    t = 1/norm(A)^2;    
    x = x0;
    y = x;
    
    if ~isinf(nf2gmax)
       hist = zeros(nf2gmax,1);
    end
    proj = @(x) max( min( x , upp ), low );
    nf=0; ng=0;
    tic
    initTime=cputime;
    while 1
        xOld= x;   

        [~,resy] = fun(y);
        gy       = fun(resy);
        x = proj( y - t*gy);
        k = nf+1;
        y = x + k/(k+3)*(x-xOld);
        [fx,~]    = fun(x);
        hist(k)   = fx;
        nf=nf+1; ng=ng+1;
        nf2g = nf+2*ng;

        sec=(cputime-initTime);
        if (nf2g>=nf2gmax | sec>=secmax)
            disp('FISTA ended')
            break;
        end
        disp(['f= ',num2str(fx),' improved at nf2g=',num2str(nf2g)])
    end
    toc
    
    case 3 % LMBOPT
     
    st.secmax  = secmax;       
    st.nf2gmax = nf2gmax;
    st.nfmax   = inf;
    st.ngmax   = inf;
    st.prt     = 0;          
    st.epsilon = 0;
    tune       = [];
    
    % call LMBOPT
    [x,f,info] = LMBOPT(fun,x0,low,upp,tune,st);
    disp('LMBOPT ended')
    info
    sec=info.sec
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



