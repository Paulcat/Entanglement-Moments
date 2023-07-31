function [sol,res] = is_sos(P,varargin)
%IS_SOS Check whether polynomial is sum-of-squares
%   [SOL,RES] = IS_SOS(P,'Name1',Value1,...) warper for Yalmip solvesos(P).
%   RES is the residual difference P - v'*G*v in case a sos decomposition
%   was found (if not, RES is the largest coefficient of P).
%
%   Supported options
%      - 'verbose': 1 | 0
%      - 'solver' : any solvers available in Yalmip
%
%   Uses Yalmip
%
%   See also lower_bound


warning('YALMIP changes its output depending on the verbose parameter!!!');


% default options
defaults = {...
	'verbose', 1, ...
	'solver', 'mosek' };

[verbose, solver] = process_options(varargin,defaults{:});


opt = sdpsettings( ... % options
	'solver', solver, ...
	'verbose', verbose);


F = sos(P); % sos constraint


[sol,~,~,res] = solvesos(F,[],opt);

end

