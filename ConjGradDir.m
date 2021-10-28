
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ConjGradDir.m %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [point,step,par,info]=ConjGradDir(fun,point,step,par,tune,info)
% construct the conjugate gradient direction

function [point,step,par,info]=ConjGradDir(fun,point,step,par,tune,info)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% get gamma %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[point,step,par,info]=getGam(fun,point,step,par,tune,info);


if ~info.done % if both nf2gmax and secmax are not exceeded
              % inside getGam, do the following

    % find subspace direction
    if point.m0>0,

        % find ingredients for projected subspace minimizer
        h=point.S(point.I,par.hist)'*point.gI;
        % do the following only if h is significant                   
        Hoh=point.H(par.hist,par.hist);
        step.q=(step.pI'*point.Y(point.I,par.hist))';
        rhs=[-h step.q];
        warning off
        sol=Hoh\rhs;
        warning on 

        ok = (any(isnan(sol))|any(isinf(sol)));
        if ~ok
           z=sol(:,1);
           step.r=sol(:,2);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%% get denom %%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [par]      = regDenom(point,step,par,tune);
            zeta       = (step.gTp+step.q'*z)/par.denom;
            okz        = isnan(zeta);
            if okz, zeta = tune.zetamin; end
            % finite only if z, zeta finite
            z          = z+zeta*step.r;
            par.cosine = -step.gTp*zeta+h'*z; 
            step.pI     = -step.pI*zeta+point.S(point.I,par.hist)*z;

            if info.prt>=1, 
                disp(' ') 
                disp(['ConjGradDir: the subspace direction has been',...
                 ' computed (m0=',num2str(point.m0),')']); 
                disp(' ') 
            end

        else
            % no subspace step possible
            % here doing instead an ''exact'' line search saves
            % function values
            zeta        = step.gTp/par.gamma;
            okz        = isnan(zeta);
            if okz, zeta = tune.zetamin; end
            zeta        = min(tune.zetamax,max(tune.zetamin,zeta));
            par.cosine  = -step.gTp*zeta; % is finite only if z is finite 
            step.pI      = -step.pI*zeta;
            if info.prt>=1, 
                disp(' ') 
                disp(['ConjGradDir: no subspace step possible;',...
                     ' m0=',num2str(point.m0)]);
                disp(' ')  
            end
        end

    else
        % no subspace step possible
        % here doing instead an ''exact'' line search saves function values

         zeta       = step.gTp/par.gamma;
         okz        = isnan(zeta);
         if okz, zeta = tune.zetamin; end
         zeta       = min(tune.zetamax,max(tune.zetamin,zeta));

         par.cosine  = -step.gTp*zeta; % is finite only if z is finite 
         step.pI      = -step.pI*zeta;

         if info.prt>=1, 
             disp(' ') 
             disp(['ConjGradDir: no subspace step possible;',...
                     ' m0=',num2str(point.m0)]); 
             disp(' ') 
          end

    end;
    
    step.p=zeros(point.n,1);step.p(point.I)=step.pI;

    point.xnew=max(point.low,min(point.x+step.p,point.upp));

    step.p=point.xnew-point.x;
    
   
   
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   