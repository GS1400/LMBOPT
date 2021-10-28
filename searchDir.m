%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% searchDir.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [point,step,par]=searchDir(point,step,par,tune,info)
% compute the search direction
%
% for details of input and output structures see lmbopt.m


function [point,step,par]=searchDir(point,step,par,tune,info)

if info.ng==1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%% call scaleDir %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [step,par]=scaleDir(point,step,par);
    step.gTp=point.gI'*step.pI;
    if info.prt>=1
       disp('')
       disp('searchDir: pinit has computed by scaleDir.m')
       disp('')
    end
    
elseif par.nlocal==tune.nwait & point.nh>=1
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% call quasiNewtonDir %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [point,step,par]=quasiNewtonDir(point,step,par,tune);
   
    step.gTp=point.gI'*step.pI;
    
    
    if info.prt>=1
       disp('')
       disp('searchDir: pinit has computed by quasiNewtonDir.m')
       disp('')
    end
    
   
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% call AvoidZigzagDir %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    [step]=AvoidZigzagDir(point,step,par,tune);
    
    step.gTp=point.gI'*step.pI;
    
    
    if info.prt>=1 
        disp('')
        disp('searchDir: pinit has computed by AvoidZigzagDir.m')
        disp('')
    end
end  

% enforce current activities and angle condition
par.firstAngle=1;
[step]=enforceAngle(point,step,par,tune,info);



step.p=zeros(point.n,1); step.p(point.I)=step.pI;

step.pinit=step.p; step.pinitI=step.pI;






 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
