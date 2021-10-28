
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% scaleCauchy.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [point,step,par]=scaleCauchy(point,step,par,tune,info)
% compute a scaled Cauchy point

function [point,step,par]=scaleCauchy(point,step,par,tune,info)

if par.freeBounds
    
     point.actlc=find(point.xnew==point.low); 
     point.actuc=find(point.xnew==point.upp);

     ok = (length(point.actl)~=length(point.actlc));
     ok = (ok|length(point.actu)~=length(point.actuc));

     if ok, par.sub=0;
     else
        par.sub = all(point.actl==point.actlc);
        par.sub = (par.sub & all(point.actu==point.actuc));
     end  

    if ~par.sub
        
        if info.prt>=1
          disp(' ')
          disp('scaleCauchy changes the activity'); 
          disp(' ')
        end

       if point.nh>=1
           D  = sqrt(sum(point.Y.^2'))./sqrt(sum(point.S.^2'));
       else
           D = abs(step.p./point.g);
      end

        Ok   = (isnan(D)|D>=tune.deltaD |D<=1/tune.deltaD);

        D(Ok) = 1;  D = D(:); 

        alpl = point.low-point.x;
        alpu = point.upp-point.x;
        step.p = zeros(size(point.x));
        ind0 = find(~D);

        indgp = find(point.g(ind0)>0);

        if ~isempty(indgp)
            step.p(indgp)=alpl(indgp);
        end

        indgn = find(point.g(ind0)<0);

         if ~isempty(indgn)
            step.p(indgn)=alpu(indgn);
        end

        ind = find(D);
        dfl = alpl(ind).*(point.g(ind)+alpl(ind).*D(ind));
        dfu = alpu(ind).*(point.g(ind)+alpu(ind).*D(ind));
        j = find(D(ind)<0);
        if ~isempty(j)
          [~,ii] = min([dfl(j) dfu(j)],[],2);
          iii = find(ii==1);
          if ~isempty(iii)
            step.p(ind(j(iii))) = alpl(ind(j(iii)));
          end
          iii = find(ii==2);
          if ~isempty(iii)
            step.p(ind(j(iii))) = alpu(ind(j(iii)));
          end
        end
        j = find(D(ind)>0);
        if ~isempty(j)
          alp = -0.5*point.g(ind(j))./D(ind(j));
          step.p(ind(j)) = max(min(alp,alpu(ind(j))),alpl(ind(j)));
        end
        point.xnew = point.x+step.p; step.pI = step.p(point.I);
    else
        par.sub=1;
    end
else
     par.sub=1;
end
step.gTp=point.gI'*step.pI;



 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   