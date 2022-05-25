function Xdot=state_eq(t,X,varargin) 
%Inputs:
%t - time
% X state [q;qdot]

%varargin inputs:
%Mchoice - [1 -> 0.5kg, 2->0.6kg]
%sample_time - time vector of samples for tauVec and FeVec
%tauVec - motor inputs 3xN
%FeVec - external forces on tool 6xN

%initalization
tauVec = zeros(3,1);
FeVec   = zeros(6,1);

%varargin inputs
Mchoice=1;
sample_time = zeros(1,1);
tau = zeros(3,1); tauFlag=0;
Fe = zeros(6,1); FeFlag=0;
for ind	= 1:2:length(varargin)
	comm	= lower(varargin{ind});
	switch comm
        case 'mchoice'
            Mchoice=varargin{ind+1};
		case 'sample_time'
			sample_time=varargin{ind+1};
		case 'tauvec'
			tauVec=varargin{ind+1};
            tauFlag=1;
		case 'Fevec'
			FeVec=varargin{ind+1};
            FeFlag=1;
	end
end

q = X(1:3);
qdot = X(4:6);

if  tauFlag
	for i = 1:3
		tau(i)=interp1(sample_time,tauVec(i,:),t);
	end
end
if  FeFlag
	for i=1:6
		Fe(i)=interp1(sample_time,FeVec(i,:),t);
	end
end

[~,~,J]=jacobian_mat(q);
[H,C,G]=dynamics_mat(q,qdot,'Mchoice',Mchoice);

Xdot = zeros(6,1);
Xdot(1:3) = qdot;
Xdot(4:6) = H\(tau+J'*Fe-G-C*qdot); %qdot2
end

