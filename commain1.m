
tic;
x0=[0.1;0.2;0.25;-0.05;0;1;0;1;0;0;
    0;0;0;0;0;0;0;0;0;0;
    0;0;0;0;0;0;0;0];
Tf=[0:0.001:40];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x]=ode5('commodel1',Tf,x0,[]);
L1=length(Tf);
for i=1:L1
    [dx,x_11(i),x_12(i),x_21(i),x_22(i),x_31(i),x_32(i),x_41(i),x_42(i),s_11(i),s_21(i),s_31(i),s_41(i),u_1(i),u_2(i),u_3(i),u_4(i)]=comcontrol1(Tf(i),x(i,:)');
end

y1=0.05*ones(1,L1);
y2=-0.05*ones(1,L1);
y3=0.5*ones(1,L1);
y4=-0.5*ones(1,L1);

figure(2)
subplot(2,1,2);
plot(Tf,s_11,'r',Tf,s_21,'b',Tf,s_31,'-.r',Tf,s_41,'-.b',Tf,y1,'k',Tf,y2,'k',Tf,y3,'-.k',Tf,y4,'-.k','linewidth',1);
legend('$h_{1,1}$','$h_{2,1}$','$h_{3,1}$','$h_{4,1}$','0.05','-0.05','0.5','-0.5');
title('(b) Under the controller designed in [58]')
xlabel(' Time(sec)');
set(gca,'FontSize',10,'Fontname', 'Times New Roman');

figure(7)
subplot(2,1,2);
plot(Tf,u_1,'r',Tf,u_2,'b',Tf,u_3,'-.r',Tf,u_4,'-.b','linewidth',1);
xlabel(' Time(sec)');
legend('$u_1$','$u_2$','$u_3$','$u_4$')
title('(b) The control input in [58]')
set(gca,'FontSize',10,'Fontname', 'Times New Roman');
