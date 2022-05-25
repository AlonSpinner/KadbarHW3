function plotrobot(q,varargin)
load('Parameters.mat');

pausetime=0;
xlimits=[-0.1,0.4];
ylimits=[-0.2,0.2];
zlimits=[-0.3,0.6];
x=[];
for ind	= 1:2:length(varargin)
    comm	= lower(varargin{ind});
    switch comm
        case 'pausetime'
            pausetime	= lower(varargin{ind+1});
        case 'xlimits'
            xlimits	= lower(varargin{ind+1});
        case 'ylimits'
            ylimits	= lower(varargin{ind+1});
        case 'zlimits'
            zlimits	= lower(varargin{ind+1});
        case 'x'
            x	= lower(varargin{ind+1});
    end
end

fig=figure('color',[1,1,1]);
ax=axes(fig);
grid(ax,'on'); hold(ax,'on');
xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');
view(ax,3); axis(ax,'equal');
ax.XLim=xlimits;
ax.YLim=ylimits;
ax.ZLim=zlimits;

if ~isempty(x) %plot track if provded
    plot3(ax,x(1,:),x(2,:),x(3,:),'linewidth',1,'color',[0,0.5,0]);
    plot3(ax,x(1,1),x(2,1),x(3,1),'marker','o','color',[0,0,0],'markersize',6);
    plot3(ax,x(1,end),x(2,end),x(3,end),'marker','*','color',[0,0,0],'markersize',6);
end

P0=[0,0,0]';
P1=[0,0,prm.H]';
A2t0=prm.f_A2t0(q(:,1));
P2=A2t0(1:3,4);
A3t0=prm.f_A3t0(q(:,1));
P3=A3t0(1:3,4);

LineColors=lines(3);
JointColors=lines(4);
L1=plot3(ax,[P0(1),P1(1)],[P0(2),P1(2)],[P0(3),P1(3)],'color',LineColors(1,:),'linewidth',2);
L2=plot3(ax,[P1(1),P2(1)],[P1(2),P2(2)],[P1(3),P2(3)],'color',LineColors(2,:),'linewidth',2);
L3=plot3(ax,[P2(1),P3(1)],[P2(2),P3(2)],[P2(3),P3(3)],'color',LineColors(3,:),'linewidth',2);
J0=scatter3(ax,P0(1),P0(2),P0(3),100,JointColors(1,:),'filled','markeredgecolor','k');
J1=scatter3(ax,P1(1),P1(2),P1(3),100,JointColors(2,:),'filled','markeredgecolor','k');
J2=scatter3(ax,P2(1),P2(2),P2(3),100,JointColors(3,:),'filled','markeredgecolor','k');
J3=scatter3(ax,P3(1),P3(2),P3(3),100,JointColors(4,:),'filled','markeredgecolor','k');


Tool=plot3(ax,P3(1),P3(2),P3(3),'linewidth',2,'linestyle','--');

for i=2:size(q,2)
    %Solve for O2
    A2t0=prm.f_A2t0(q(:,i));
    P2=A2t0(1:3,4);
    L2.XData=[0,P2(1)];
    L2.YData=[0,P2(2)];
    L2.ZData=[prm.H,P2(3)];
    J2.XData=P2(1);
    J2.YData=P2(2);
    J2.ZData=P2(3);
    %Solve for O3
    A3t0=prm.f_A3t0(q(:,i));
    P3=A3t0(1:3,4);
    L3.XData=[P2(1),P3(1)];
    L3.YData=[P2(2),P3(2)];
    L3.ZData=[P2(3),P3(3)];
    J3.XData=P3(1);
    J3.YData=P3(2);
    J3.ZData=P3(3);
    
    Tool.XData=[Tool.XData,P3(1)];
    Tool.YData=[Tool.YData,P3(2)];
    Tool.ZData=[Tool.ZData,P3(3)];
    
    drawnow;
    pause(pausetime);
end

end