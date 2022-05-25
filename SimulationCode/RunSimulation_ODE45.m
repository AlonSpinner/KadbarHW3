function [tsim,q_actual,qdot_actual]=RunSimulation_ODE45(tauVec,q0,qdot0,sample_time,integrate_time,Mchoice)
X0=[q0;qdot0];
[tsim,y] = ode45(@(t,y) state_eq(t,y,'tauVec',tauVec,'sample_time',sample_time,'Mchoice',Mchoice),...
    integrate_time,X0);
y=y'; %lets have observations in columns instead of rows
q_actual=y(1:3,:);
qdot_actual=y(4:6,:);
end