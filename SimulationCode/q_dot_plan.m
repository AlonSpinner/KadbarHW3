function q_dot=q_dot_plan(prof,t,varargin)
%---input
%prof - motion profile: 'const','trapz','poly'
%t - scalar
%varargin input:
%method - 'jacobian'/'numeric' for velocity computation
%elbows - [t2_elbow,d3_eblow]: accepts -1 or 1 values

%---output
%q_dot - 3XN matrix of [t1,t2,d3]' joint locations in [radians/m]

method='jacobian';
elbows  = [1,1]; 
for ind	= 1:2:length(varargin)
	comm	= lower(varargin{ind});
	switch comm
        case 'method'
			method	= lower(varargin{ind+1});
		case 'elbows'
			elbows	= lower(varargin{ind+1});
	end
end

load('Parameters.mat');
dt=prm.dt;

x=x_plan(prof,t);
q=inverse_kin(x,elbows);
switch method
    case 'numeric'    
        q_dot=gradient(q,dt);
    case 'jacobian'
        v=v_plan(prof,t);
        N=size(v,2);
        [JL,~]=jacobian_mat(q);
        q_dot=zeros(3,N);
        for i=1:N
            q_dot(:,i)=linsolve((JL(:,:,i)),v(:,i));
        end
end
end