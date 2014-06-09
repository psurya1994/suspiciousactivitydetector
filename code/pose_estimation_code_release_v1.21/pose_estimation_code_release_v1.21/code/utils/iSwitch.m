function v = iSwitch(varargin)

assert(rem(length(varargin),2)==1);

for i = 1:2:length(varargin)-1
	if varargin{i}
		v = varargin{i+1};
		return;
	end;
end; clear i;

v = varargin{end};
return;
