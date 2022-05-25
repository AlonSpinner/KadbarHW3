function q_dot2=q_dot2_plan(prof,t,varargin)
%---input
%prof - motion profile: 'const','trapz','poly'
%t - scalar
%varargin input:
%method - 'jacobian'/'numeric' for velocity computation
%elbows - [t2_elbow,d3_eblow]: accepts -1 or 1 values

%---output
%q_dot2 - 3XN matrix of [t1,t2,d3]' joint locations in [radians/m]

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
q_dot=q_dot_plan(prof,t,varargin{:});
switch method
    case 'numeric'  
        q_dot2=gradient(q_dot,dt);
    case 'jacobian'
        [JL,~]=jacobian_mat(q);
        [JL_dot,~]=jacobian_dot_mat(q,q_dot);
        %compute a
        a=a_plan(prof,t);
        %compute qdot2=inv(JL)*(a-JLdot*qdot)
        q_dot2=zeros(3,size(a,2));
        for i=1:size(a,2)
            q_dot2(:,i)=linsolve(JL(:,:,i),a(:, i)-JL_dot(:,:,i)*q_dot(:,i));
        end
end
end