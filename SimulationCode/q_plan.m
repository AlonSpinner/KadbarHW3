function q=q_plan(prof,t,varargin)
%---input
%prof - motion profile: 'const','trapz','poly'
%t - scalar
%varargin input:
%elbows - [t2_elbow,d3_eblow]: accepts -1 or 1 values

%---output
%q - 3XN matrix of [t1,t2,d3]' joint locations in [radians/m]

elbows  = [1,1]; 
for ind	= 1:2:length(varargin)
	comm	= lower(varargin{ind});
	switch comm
		case 'elbows'
			elbows	= lower(varargin{ind+1});
	end
end
x=x_plan(prof,t);
q=inverse_kin(x,elbows);
end