%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% enforceAngle.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [step]=enforceAngle(point,step,par,tune,info)
% enforce the angle condition
%
% for details of input and output structures see lmbopt.m
%


function [step]=enforceAngle(point,step,par,tune,info)

if ~par.firstAngle, 
    step.pI=step.p(point.I); step.gTp=point.gI'*step.pI;
end

if isnan(step.gTp)
  if ~isnan(step.lamBB)
       if par.nlocal==tune.nwait & point.nh>=1
          step.lamBB  = max(point.Lam);
       end
       step.pI     = -step.lamBB*point.gI;
       step.gTp    = step.pI'*point.gI;
   else
       step.pI = -point.gI; step.gTp=point.gI'*step.pI;
   end
else
    if par.firstAngle

        if step.gTp>0,
          % move away from maximizer or saddle point
           step.pI=-step.pI;
           step.gTp=point.gI'*step.pI; 
        end

        if abs(step.gTp) <= tune.Deltapg*(abs(step.pI)'*abs(point.gI))
           if ~isnan(step.lamBB)
              if par.nlocal==tune.nwait & point.nh>=1
                 step.lamBB  = max(point.Lam);
              end
               step.pI     = -step.lamBB*point.gI;
               step.gTp    = step.pI'*point.gI;
           else
               step.pI = -point.gI; step.gTp=point.gI'*step.pI;
           end 
            
         end

    else
        if step.gTp>=0,
          % move away from maximizer or saddle point
          step.pI=step.pinitI; step.gTp = step.pI'*point.gI; 
          step.p=step.pinit;
        end
    end

end


ok = (par.firstAngle | ~par.firstAngle & step.gTp <0);


if ok
    sigma1  = point.gI'*point.gI;
    sigma2  = step.pI'*step.pI;
    sigma12 = sigma1*sigma2;
    sigmanew = step.gTp/sqrt(sigma12);

    if ~isnan(sigmanew)& isfinite(sigmanew)

        if sigmanew<=-tune.Deltaangle & step.gTp<0

           if ~isnan(step.gTp)
              if abs(step.gTp) <= tune.Deltapg*...
                         (abs(step.pI)'*abs(point.gI))
                 step.pI = -point.gI;  step.gTp = step.pI'*point.gI;
              else
                   if info.prt>=1, 
                      disp(' ') 
                      disp('The angle condition was satisfied');
                      disp(' ')  
                   end
              end
           else
                step.pI = -point.gI; 
                step.gTp=point.gI'*step.pI;
           end

        else
          w = (sigma12*max(tune.Deltaw,1-sigmanew^2))/...
               (1-tune.Deltaangle^2);
          t = (step.gTp+tune.Deltaangle*sqrt(w))/sigma1;
          if w>0 & isfinite(t) & t>=0, 
              step.pI=step.pI-t*point.gI; 
              step.gTp = step.pI'*point.gI;
               if ~isnan(step.gTp)
                  if abs(step.gTp) <= ...
                             tune.Deltapg*(abs(step.pI)'*abs(point.gI))
                     step.pI = -point.gI;  
                     step.gTp = step.pI'*point.gI;
                  else
                      if info.prt>=1 
                        disp(' ') 
                        disp('The angle condition was enforced');
                        disp(' ')  
                      end
                      
                  end
               else
                   if ~isnan(step.lamBB)
                       if par.nlocal==tune.nwait & point.nh>=1
                           step.lamBB  = max(point.Lam);
                       end
                       step.pI     = -step.lamBB*point.gI;
                       step.gTp    = step.pI'*point.gI;
                   else
                    step.pI = -point.gI; step.gTp=point.gI'*step.pI;
                   end
               end
          else 
            % this includes the case where pI=0  

            if ~isnan(step.lamBB)
                if par.nlocal==tune.nwait & point.nh>=1
                   step.lamBB  = max(point.Lam);
                end
                step.pI     = -step.lamBB*point.gI;
                step.gTp    = step.pI'*point.gI;
            else
                if ~par.firstAngle
                   step.pI    = step.pinitI;step.p=step.pinit;
                   step.gTp   = step.pI'*point.gI;
                   if ~isnan(step.gTp)
                     if abs(step.gTp) <= ...
                          tune.Deltapg*(abs(step.pI)'*abs(point.gI))
                        step.pI = -point.gI;  
                        step.gTp = step.pI'*point.gI;
                     end
                   else
                      step.pI = -point.gI; 
                      step.gTp=point.gI'*step.pI;
                   end
                end
            end
          end
        end

    else
        if ~isnan(step.lamBB)
            if par.nlocal==tune.nwait & point.nh>=1
               step.lamBB  = max(point.Lam);
            end
            step.pI     = -step.lamBB*point.gI;
            step.gTp    = step.pI'*point.gI;
        else
            if ~par.firstAngle
               step.pI     = step.pinitI;step.p=step.pinit;
               step.gTp   = step.pI'*point.gI;
               if ~isnan(step.gTp)
                 if abs(step.gTp) <= ...
                        tune.Deltapg*(abs(step.pI)'*abs(point.gI))
                  step.pI = -point.gI;  step.gTp = step.pI'*point.gI;
                 end
               else
                  step.pI = -point.gI; 
                  step.gTp=point.gI'*step.pI;
               end
            end
        end

    end
end % of if sigmanew
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   