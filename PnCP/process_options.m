function varargout = process_options(params,varargin)
% A simplified version of the parseparams code in octave
% author:
% Jean-Francois Cardoso
%
% original authors:
% Alexander Barth
% Aida Alvera

% default options
names  = varargin(1:2:end);
defvalues = varargin(2:2:end);

% specified options
pnames  = params(1:2:end);
values = params(2:2:end);
if (length(pnames)~=length(values)) || ~iscellstr(pnames)
	error('options must be given as name-value pairs');
end

% match
varargout = defvalues;
for i = 1:length(pnames)
	pname = pnames{i};
	id = find(strcmp(names,pname));

	% set values
	if id==0
		error('unknown option: %s', pname);
	end
	varargout{id} = values{i};
end

end