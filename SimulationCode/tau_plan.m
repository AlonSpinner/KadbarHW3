function tau=tau_plan(prof,t,varargin)
%---input
%prof - motion profile: 'const','trapz','poly'
%t - vector of time values
%varargin input:
%elbows - [t2_elbow,d3_eblow]: accepts -1 or 1 values
%Mchoice - [1 if M=0.5Kg, 2 if M=0.6Kg]

N=length(t);

elbows  = [1,1]; 
Mchoice = 1;
Fe=zeros(6,N);
for ind	= 1:2:length(varargin)
	comm	= lower(varargin{ind});
	switch comm
		case 'elbows'
            elbows	= lower(varargin{ind+1});
        case 'mchoice'
            Mchoice	= lower(varargin{ind+1});
        case 'Fe'
            Fe=lower(varargin{ind+1});
	end
end

q=q_plan(prof,t,'elbows',elbows);
q_dot=q_dot_plan(prof,t,'elbows',elbows,'method','jacobian');
q_dot2=q_dot2_plan(prof,t,'elbows',elbows,'method','jacobian');

[H,C,G]=dynamics_mat(q,q_dot,'Mchoice',Mchoice);
[~,~,J]=jacobian_mat(q);

tau=zeros(3,N);
for i=1:N
    tau(:,i)=H(:,:,i)*q_dot2(:,i)+C(:,:,i)*q_dot(:,i)+G(:,i)-J(:,:,i)'*Fe(:,i);
end