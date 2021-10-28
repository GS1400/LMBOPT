
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% updateSubspace.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [point]=updateSubspace(point,step,par,tune)
% update the subspace such as the memory matrices S, Y and H
%
% for details of input and output structures see lmbopt.m

function [point]=updateSubspace(point,step,par,tune)


 % is close to 0 when step was tiny
   
flagnull  = (par.nnull==0); 
if  flagnull 
   gTy        = point.g'*point.y;
   powell     = abs(gTy)/point.omega;
   flagpowell = powell>=tune.Deltapo; 
   flaglocal  = (par.nlocal<=tune.nwait);
   subOk      = flagnull & (~flaglocal | flagpowell);
   if subOk % update the memory matrices S, Y, and H
        if point.ch<point.m, point.ch=point.ch+1; else, point.ch=1; end;
        point.S(:,point.ch)=step.s; point.Y(:,point.ch)=point.y;
        if tune.typeH, point.H(point.ch,:)=step.s'*point.Y;
        else, point.H(point.ch,:)=point.y'*point.S; 
        end;
        point.H(:,point.ch)=point.H(point.ch,:)';
        point.nh=point.nh+1;
   end
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
