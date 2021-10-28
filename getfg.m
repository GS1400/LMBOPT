%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% getfg.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% input
%      x        % point
%      prob     % ingredients of test problem
%
% output   
%       f     % function value
%       res   % gradient vector
%
% [f,res] = getfg(x,prob) compute f and res at x
% g = grtfg(res) computes the gradient at x
% 

function [f,res] = getfg(x,prob)

A=prob.A; b=prob.b; 

if nargout==2
   res = A*x-b; f = 0.5*(res'*res); 
elseif nargout==1 % x contains res, f returns g
    f = A'*x;
else
    error('at least one output needs')
end
     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%