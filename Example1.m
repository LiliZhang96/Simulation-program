function tys
l=[-16;-64];

beta=[45;45];
xi=[0.01;0.01];
rho=[1;1];

lags=[1,1,1,1]*2.5;

sol=dde23(@tys,lags,@tyshist,[0,40],[],xi,l,rho,beta);

figure(1)
plot(sol.x,sol.y(5,:),'r',sol.x,sol.y(18,:),'b',sol.x,sol.y(28,:),'-.r',sol.x,sol.y(38,:),'-.b',sol.x,sol.yp(9,:),'-.k','linewidth',1);
xlabel(' Time(sec)');
legend('$y_{1}$', '$y_{2}$', '$y_{3}$', '$y_{4}$','$y_{d_1}$')
set(gca,'FontSize',10,'Fontname', 'Times New Roman');


figure(2)
subplot(2,1,1);
plot(sol.x,sol.yp(11,:),'r',sol.x,sol.yp(23,:),'b',sol.x,sol.yp(33,:),'-.r',sol.x,sol.yp(43,:),'-.b',sol.x,sol.yp(12,:),'-k',sol.x,sol.yp(13,:),'-k','linewidth',1);
xlabel(' Time(sec)');
legend('$h_{1,1}$','$h_{2,1}$','$h_{3,1}$','$h_{4,1}$','$-D(t)$','$D(t)$')
title('(a) Under the proposed controller')
set(gca,'FontSize',10,'Fontname', 'Times New Roman');
hold on;

figure(3)
plot(sol.x,sol.y(3,:),'r',sol.x,sol.y(16,:),'b',sol.x,sol.y(26,:),'-.r',sol.x,sol.y(36,:),'-.b','linewidth',1);
xlabel(' Time(sec)');
legend('$e_{1,1}$', '$e_{2,1}$', '$e_{3,1}$', '$e_{4,1}$')

figure(4)
plot(sol.x,sol.y(4,:),'r',sol.x,sol.y(17,:),'b',sol.x,sol.y(27,:),'-.r',sol.x,sol.y(37,:),'-.b','linewidth',1);
xlabel(' Time(sec)');
legend('$e_{1,2}$', '$e_{2,2}$', '$e_{3,2}$', '$e_{4,2}$')
set(gca,'FontSize',10,'Fontname', 'Times New Roman');


figure(5)
plot(sol.x,sol.y(7,:),'r',sol.x,sol.y(20,:),'b',sol.x,sol.y(30,:),'-.r',sol.x,sol.y(40,:),'-.b','linewidth',1);
xlabel(' Time(sec)');
legend('$\hat{\varphi}_{1,1}$','$\hat{\varphi}_{2,1}$','$\hat{\varphi}_{3,1}$','$\hat{\varphi}_{4,1}$')

figure(6)
plot(sol.x,sol.y(8,:),'r',sol.x,sol.y(21,:),'b',sol.x,sol.y(31,:),'-.r',sol.x,sol.y(41,:),'-.b','linewidth',1);
xlabel(' Time(sec)');
legend('$\hat{\varphi}_{1,2}$','$\hat{\varphi}_{2,2}$','$\hat{\varphi}_{3,2}$','$\hat{\varphi}_{4,2}$')
set(gca,'FontSize',10,'Fontname', 'Times New Roman');


figure(7)
plot(sol.x,sol.yp(10,:),'r',sol.x,sol.yp(22,:),'b',sol.x,sol.yp(32,:),'-.r',sol.x,sol.yp(42,:),'-.b','linewidth',1);
xlabel(' Time(sec)');
legend('$u_1$','$u_2$','$u_3$','$u_4$')
set(gca,'FontSize',10,'Fontname', 'Times New Roman');

function dydt=tys(t,y,Z,xi,l,rho,beta)
ylag1 = Z(:,1);
ylag2 = Z(:,2);
ylag3 = Z(:,3);
ylag4 = Z(:,4);
%
xi11=xi(1);
xi12=xi(2);
xi21=xi(1);
xi22=xi(2);
xi31=xi(1);
xi32=xi(2);
xi41=xi(1);
xi42=xi(2);

l1=l(1);
l2=l(2);
rho11=rho(1);
rho12=rho(2);
beta11=beta(1);
beta12=beta(2);

rho21=rho(1);
rho22=rho(2);
beta21=beta(1);
beta22=beta(2);

rho31=rho(1);
rho32=rho(2);
beta31=beta(1);
beta32=beta(2);

rho41=rho(1);
rho42=rho(2);
beta41=beta(1);
beta42=beta(2);

yd=sin(t);

r1=1;T=5;q=0.05;
DL1=(r1*(sin((pi/(2*T))*(T-t)))^3)*(t>=0&t<=T)+0*(t>T);
D1=DL1+q;
dotD=(-((3*r1*pi)/(2*T))*cos((pi/(2*T))*(T-t))*(sin((pi/(2*T))*(T-t)))^2)*(t>=0&t<=T)+0*(t>T);
center=[-pi/2;0;pi/2];

h=[y(5);yd];
S11=FuzzySet(h,center);
g11=1+norm(S11);

h=[y(5);yd;y(7);D1];
S12=FuzzySet(h,center);
g12=1+norm(S12);

h=[y(18);y(5);y(6)];
S21=FuzzySet(h,center);
g21=1+norm(S21);

h=[y(15);y(18);y(20);y(5);y(2);D1];
S22=FuzzySet(h,center);
g22=1+norm(S22);

h=[y(28);y(5);y(6)];
S31=FuzzySet(h,center);
g31=1+norm(S31);

h=[y(25);(28);y(30);y(5);y(2);D1];
S32=FuzzySet(h,center);
g32=1+norm(S32);

h=[y(38);y(18);y(19);y(28);y(29)];
S41=FuzzySet(h,center);
g41=1+norm(S41);

h=[y(35);y(38);y(40);y(15);y(18);y(25);y(28);D1];
S42=FuzzySet(h,center);
g42=1+norm(S42);

h11=y(5)-y(18);
f=-(2*D1*h11^3)/((D1*D1-h11*h11)^2);
F=dotD*f; 
gamma1=(D1*D1*(D1*D1+h11*h11))/((D1*D1-h11*h11)^2);
Phi1=(D1*D1*h11)/((D1-h11)*(D1+h11)); 
alph11=-F/gamma1-((beta11*Phi1)/gamma1)-0.25*Phi1*gamma1-y(7)*g11*tanh(Phi1*gamma1*g11/0.01);
law11=-xi11*y(7)+rho11*Phi1*gamma1*g11*tanh(Phi1*gamma1*g11/0.01);
h12=y(2)-alph11;
u1=-beta12*h12-y(8)*g12*tanh(h12*g12/0.01);
law12=-xi12*y(8)+rho12*h12*g12*tanh(h12*g12/0.01);

h21=y(18)-yd;
f=-(2*D1*h21^3)/((D1*D1-h21*h21)^2);
F=dotD*f; 
gamma2=(D1*D1*(D1*D1+h21*h21))/((D1*D1-h21*h21)^2);
Phi2=(D1*D1*h21)/((D1-h21)*(D1+h21)); 
alph21=-F/gamma2-((beta21*Phi2)/gamma2)-0.25*Phi2*gamma2-y(20)*g21*tanh(Phi2*gamma2*g21/0.01);
law21=-xi21*y(20)+rho21*Phi2*gamma2*g21*tanh(Phi2*gamma2*g21/0.01);
h22=y(15)-alph21;
u2=-beta22*h22-y(21)*g22*tanh(h22*g22/0.01);
law22=-xi22*y(21)+rho22*h22*g22*tanh(h22*g22/0.01);

h31=y(28)-y(18);
f=-(2*D1*h31^3)/((D1*D1-h31*h31)^2);
F=dotD*f; 
gamma3=(D1*D1*(D1*D1+h31*h31))/((D1*D1-h31*h31)^2);
Phi3=(D1*D1*h31)/((D1-h31)*(D1+h31)); 
alph31=-F/gamma3-((1.4*beta31*Phi3)/gamma3)-0.25*Phi3*gamma3-y(30)*g31*tanh(Phi3*gamma3*g31/0.01);
law31=-xi31*y(30)+rho31*Phi3*gamma3*g31*tanh(Phi3*gamma3*g31/0.01);
h32=y(25)-alph31;
u3=-1.4*beta32*h32-y(31)*g32*tanh(h32*g32/0.01);
law32=-xi32*y(31)+rho32*h32*g32*tanh(h32*g32/0.01);

h41=y(38)-y(28);
f=-(2*D1*h41^3)/((D1*D1-h41*h41)^2);
F=dotD*f; 
gamma4=(D1*D1*(D1*D1+h41*h41))/((D1*D1-h41*h41)^2);
Phi4=(D1*D1*h41)/((D1-h41)*(D1+h41)); 
alph41=-F/gamma4-((beta41*Phi4)/gamma4)-0.5*Phi4*gamma4-y(40)*g41*tanh(Phi4*gamma4*g41/0.01);
law41=-xi41*y(40)+rho41*Phi4*gamma4*g41*tanh(Phi4*gamma4*g41/0.01);
h42=y(35)-alph41;
u4=-beta42*h42-y(41)*g42*tanh(h42*g42/0.01);
law42=-xi42*y(41)+rho42*h42*g42*tanh(h42*g42/0.01);

obs11=y(2)-l1*(y(5)-y(1));
obs12=u1-l2*(y(5)-y(1));
obs21=y(15)-l1*(y(18)-y(14));
obs22=u2-l2*(y(18)-y(14));
obs31=y(25)-l1*(y(28)-y(24));
obs32=u3-l2*(y(28)-y(24));
obs41=y(35)-l1*(y(38)-y(34));
obs42=u4-l2*(y(38)-y(34));

obe11=y(4)+(l1+1)*y(3)+y(5)/(1+y(5)^4);
obe12=l2*y(3)+0.1*sin(y(5)-(y(2)+y(4)))*exp(-(y(5)^2+(y(2)+y(4))^4));
obe21=y(17)+(l1+1)*y(16)+y(18)/(1+y(18)^4);
obe22=l2*y(16)+0.15*sin(y(18)-(y(15)+y(17)))-exp(-(y(18)^2+(y(15)+y(17))^4));
obe31=y(27)+(l1+1)*y(26)+y(28)/(1+y(28)^4);
obe32=l2*y(26)+0.15*sin(y(28)-(y(25)+y(27)))-exp(-(y(28)^2+(y(25)+y(27))^4));
obe41=y(37)+(l1+1)*y(36)+y(38)/(1+y(38)^2);
obe42=l2*y(36)+0.2*sin(y(38)-(y(35)+y(37)))-exp(-(y(38)^2+(y(35)+y(37))^4));

g11=y(2)+y(4)+y(5)/(1+y(5)^4);
g12=y(2)+y(4);
g21=y(15)+y(17)+y(18)/(1+y(18)^4);
g22=y(15)+y(17);
g31=y(25)+y(27)+y(28)/(1+y(28)^4);
g32=y(25)+y(27);
g41=y(35)+y(37)+y(38)/(1+y(38)^2);
g42=y(35)+y(37);

  dydt = [obs11;obs12;obe11;obe12;g11;g12;law11;law12;yd;u1;h11;-D1;D1;obs21;obs22;obe21;obe22;g21;g22;law21;law22;u2;h21;obs31;obs32;obe31;obe32;g31;g32;law31;law32;u3;h31;obs41;obs42;obe41;obe42;g41;g42;law41;law42;u4;h41];
 
function S = tyshist(t,xi,l,rho,beta)

 S = [0.1;0.2;0;0;0.1;0.2;0;0;0;0;0;-1.05;1.05;0.25;-0.05;0;0;0.25;-0.05;0;0;0;0;0.1;-0.15;0;0;0.1;-0.15;0;0;0;0;-0.25;0.3;0;0;-0.25;0.3;0;0;0;0];
