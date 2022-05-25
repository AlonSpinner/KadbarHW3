%% Main
clear
clc

load('Parameters.mat');
T=prm.T;
dt=prm.dt;
% sample_time=0:dt:T;
sample_time=0:dt:T;
integrate_time=0:dt/10:T;
elbows=[1,1];
prof='poly';
%% Poly xyz
x=x_plan(prof,sample_time);
v=v_plan(prof,sample_time);
a=a_plan(prof,sample_time);

fig=figure('color',[1,1,1],'position',[300,100,1100,600]);
ax11=subplot(1,3,1,'parent',fig); 
hold(ax11,'on'); grid(ax11,'on'); xlabel(ax11,'[s]'); ylabel(ax11,'[m]'); title('Position')
ax12=subplot(1,3,2,'parent',fig);
hold(ax12,'on'); grid(ax12,'on'); xlabel(ax12,'[s]'); ylabel(ax12,'[m/s]'); title('Velocity')
ax13=subplot(1,3,3,'parent',fig); 
hold(ax13,'on'); grid(ax13,'on'); xlabel(ax13,'[s]'); ylabel(ax13,'[m/s^2]'); title('Acceleration')

plot(ax11,sample_time,x','linewidth',2);
plot(ax12,sample_time,v','linewidth',2);
plot(ax13,sample_time,a','linewidth',2);

legend(ax11,'x','y','z');
sgtitle(fig,'Polynomial Trajectory - Tool Kinematics');
%% Q2 - Calcualte tau for polynomial profile
tauVec=tau_plan(prof,sample_time,'elbows',elbows,'Mchoice',1);
%% Plot Q2
fig=figure('color',[1,1,1],'position',[300,100,1100,600]);
Titles={'\theta_1','\theta_2','d_3'};
Ylabels={'[Nm]','[Nm]','[N]'};
for i=1:3
subplot(1,3,i,'parent',fig); 
hold('on'); grid('on'); xlabel('[s]'); title(Titles{i}); ylabel(Ylabels{i});
plot(sample_time,tauVec(i,:)','linewidth',2);
end
%% Solve Q3
[Rz,tauVecNewton]=Newton_Euler(prof,sample_time,'elbows',elbows,'Mchoice',1);
%% Check Q3 tau vs Q2 tau
fig=figure('color',[1,1,1],'position',[300,100,1100,600]);
Titles={'\theta_1','\theta_2','d_3'};
Ylabels={'[Nm]','[Nm]','[N]'};
for i=1:3
subplot(1,3,i,'parent',fig); 
hold('on'); grid('on'); xlabel('[s]'); title(Titles{i}); ylabel(Ylabels{i});
plot(sample_time,tauVec(i,:)','linewidth',2);
plot(sample_time,tauVecNewton(i,:)','linewidth',2);
legend('Lagrange','Newton-Euler');
end
%% Plot Rz solution - Q3
fig=figure('color',[1,1,1],'position',[300,100,1100,600]);
plot(sample_time,Rz); 
grid('on'); xlabel('[s]'); ylabel('N');
title('Reaction Force from ground in the z direction');
%% Solve Q4 - ODE45
q_planned=q_plan(prof,sample_time,'elbows',elbows);
qdot_planned=q_dot_plan(prof,sample_time,'elbows',elbows,'method','jacobian');
q0=q_planned(:,1);
qdot0=qdot_planned(:,1);

Mchoice=1;

[tsim,q_actual,qdot_actual]=RunSimulation_ODE45(tauVec,q0,qdot0,sample_time,integrate_time,Mchoice);
%% Plot Q4 - ODE45
%Plot Joints
fig=figure('color',[1,1,1],'position',[300,100,1100,600]);
Titles={'\theta_1','\theta_2','d_3'};
Ylabels={'[degree]','[degree]','[m]'};
for i=1:3
subplot(1,3,i,'parent',fig); 
hold('on'); grid('on'); xlabel('[s]'); title(Titles{i}); ylabel(Ylabels{i});
if i==1 || i==2
    plot(sample_time,rad2deg(q_planned(i,:)'),'linewidth',2);
    plot(tsim,rad2deg(q_actual(i,:)'),'linewidth',2);
else
    plot(sample_time,q_planned(i,:)','linewidth',2);
    plot(tsim,q_actual(i,:)','linewidth',2);
end
legend('Planned','Actual');
end

%Plot WYZ Error
xyz_actual=forward_kin(q_actual);
q_planned=q_plan(prof,integrate_time);
xyz_planned=forward_kin(q_planned);

arclen_planned=arclength(xyz_planned');
Err_percent=vecnorm(xyz_actual-xyz_planned,2,1)/arclen_planned*100;

fig=figure('color',[1,1,1],'position',[300,100,1100,600]);
title('Error: |Actual-Planned|/Planned ArchLength');
ylabel('Error [%]'); xlabel('time [s]'); grid('on'); hold('on');
plot(integrate_time,Err_percent)
%% Solve Q4 - Simulink
[tsim,tau_simulink,q_actual,qdot_actual,xyz_actual,q_planned,xyz_planned]=RunSimulation_Simulink();
%% Plot Q4 - Simulink
%Plot Joints
fig=figure('color',[1,1,1],'position',[300,100,1100,600]);
Titles={'\theta_1','\theta_2','d_3'};
Ylabels={'[degree]','[degree]','[m]'};
for i=1:3
subplot(1,3,i,'parent',fig); 
hold('on'); grid('on'); xlabel('[s]'); title(Titles{i}); ylabel(Ylabels{i});
if i==1 || i==2
    plot(tsim,rad2deg(q_planned(i,:)'),'linewidth',2);
    plot(tsim,rad2deg(q_actual(i,:)'),'linewidth',2,'linestyle','--');
else
    plot(tsim,q_planned(i,:)','linewidth',2);
    plot(tsim,q_actual(i,:)','linewidth',2,'linestyle','--');
end
legend('Planned','Actual');
end

%Plot WYZ Error
arclen_planned=arclength(xyz_planned');
Err_percent=vecnorm(xyz_actual-xyz_planned,2,1)/arclen_planned*100;

fig=figure('color',[1,1,1],'position',[300,100,1100,600]);
title('Error: |Actual-Planned|/Planned ArchLength');
ylabel('Error [%]'); xlabel('time [s]'); grid('on'); hold('on');
plot(tsim,Err_percent,'linewidth',2);
%% Compare tau_simulink to tau_Q2
fig=figure('color',[1,1,1],'position',[300,100,1100,600]);
Titles={'\theta_1','\theta_2','d_3'};
Ylabels={'[Nm]','[Nm]','[N]'};
for i=1:3
subplot(1,3,i,'parent',fig); 
hold('on'); grid('on'); xlabel('[s]'); title(Titles{i}); ylabel(Ylabels{i});
plot(sample_time,tauVec(i,:)','linewidth',2');
plot(tsim,tau_simulink(i,:),'linewidth',2,'linestyle','--');
legend('MATLAB','SIMULINK');
end
sgtitle('Torque/Force in joints: Simulation Comparison');
%% Q5 
%first thing:
%Change Z0 in Initalize/t2xyz to 0.26 from 0.25.
%THIS DOES NOT HAPPEN AUTOMATICLY. DONT FORGET TO CHANGE IT BACK

[tsim,tau_simulink,q_actual,qdot_actual,xyz_actual,q_planned,xyz_planned]=RunSimulation_Simulink();
%% Plot Q5 - Simulink
%Plot Joints
fig=figure('color',[1,1,1],'position',[300,100,1100,600]);
Titles={'\theta_1','\theta_2','d_3'};
Ylabels={'[degree]','[degree]','[m]'};
for i=1:3
subplot(1,3,i,'parent',fig); 
hold('on'); grid('on'); xlabel('[s]'); title(Titles{i}); ylabel(Ylabels{i});
if i==1 || i==2
    plot(tsim,rad2deg(q_planned(i,:)'),'linewidth',2);
    plot(tsim,rad2deg(q_actual(i,:)'),'linewidth',2,'linestyle','--');
else
    plot(tsim,q_planned(i,:)','linewidth',2);
    plot(tsim,q_actual(i,:)','linewidth',2,'linestyle','--');
end
legend('Planned','Actual');
end

%Plot WYZ Error
arclen_planned=arclength(xyz_planned');
Err_percent=vecnorm(xyz_actual-xyz_planned,2,1)/arclen_planned*100;

fig=figure('color',[1,1,1],'position',[300,100,1100,600]);
title('Error: |Actual-Planned|/Planned ArchLength');
ylabel('Error [%]'); xlabel('time [s]'); grid('on'); hold('on');
plot(tsim,Err_percent,'linewidth',2);
%% Q6
%first thing:
%Change Mchoice_actual to 2 in he initalization function
%DONT FORGET TO CHANGE IT BACK

[tsim,tau_simulink,q_actual,qdot_actual,xyz_actual,q_planned,xyz_planned]=RunSimulation_Simulink();
%% Plot Q6 - Simulink
%Plot Joints
fig=figure('color',[1,1,1],'position',[300,100,1100,600]);
Titles={'\theta_1','\theta_2','d_3'};
Ylabels={'[degree]','[degree]','[m]'};
for i=1:3
subplot(1,3,i,'parent',fig); 
hold('on'); grid('on'); xlabel('[s]'); title(Titles{i}); ylabel(Ylabels{i});
if i==1 || i==2
    plot(tsim,rad2deg(q_planned(i,:)'),'linewidth',2);
    plot(tsim,rad2deg(q_actual(i,:)'),'linewidth',2,'linestyle','--');
else
    plot(tsim,q_planned(i,:)','linewidth',2);
    plot(tsim,q_actual(i,:)','linewidth',2,'linestyle','--');
end
legend('Planned','Actual');
end

%Plot WYZ Error
arclen_planned=arclength(xyz_planned');
Err_percent=vecnorm(xyz_actual-xyz_planned,2,1)/arclen_planned*100;

fig=figure('color',[1,1,1],'position',[300,100,1100,600]);
title('Error: |Actual-Planned|/Planned ArchLength');
ylabel('Error [%]'); xlabel('time [s]'); grid('on'); hold('on');
plot(tsim,Err_percent,'linewidth',2);
%% Functions
function  arclen=arclength(data)
seglen = sqrt(sum(diff(data,[],1).^2,2));
arclen = sum(seglen);
end