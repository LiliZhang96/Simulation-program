function tys

lam=[1;5;5];

l=-[8.6485;336.2915;-290.6175];
d=[0.5;0.5;0]*2;
k=[30;15;15];
lags=[1,1,1,1]*2.5; 

sol=dde23(@tys,lags,@tyshist,[0,40],[],lam,l,d,k);

figure(1)
plot(sol.x,sol.y(7,:),'r',sol.x,sol.y(22,:),'-.r',sol.x,sol.y(34,:),'-b',sol.x,sol.y(46,:),'-.b',sol.x,sol.yp(11,:),'-.k','linewidth',1);
xlabel(' Time(sec)');
legend('$y_1$', '$y_2$','$y_3$','$y_4$','$y_d$')
set(gca,'FontSize',10,'Fontname', 'Times New Roman');

figure(2)
subplot(2,1,1)
plot(sol.x,sol.yp(13,:),'r',sol.x,sol.yp(27,:),'-.r',sol.x,sol.yp(39,:),'b',sol.x,sol.yp(51,:),'-.b',sol.x,sol.yp(14,:),'-.k',sol.x,sol.yp(15,:),'-.k','linewidth',1);
xlabel(' Time(sec)');
legend('$h_{1,1}$','$h_{2,1}$','$h_{3,1}$','$h_{4,1}$','$-D(t)$','$D(t)$')
set(gca,'FontSize',10,'Fontname', 'Times New Roman');

figure(3)
plot(sol.x,sol.y(4,:),'r',sol.x,sol.y(19,:),'-.r',sol.x,sol.y(31,:),'b',sol.x,sol.y(43,:),'-.b','linewidth',1);
xlabel(' Time(sec)');
legend('$e_{1,1}$','$e_{2,1}$','$e_{3,1}$','$e_{4,1}$')
set(gca,'FontSize',10,'Fontname', 'Times New Roman');

figure(4)
plot(sol.x,sol.y(5,:),'r',sol.x,sol.y(20,:),'-.r',sol.x,sol.y(32,:),'b',sol.x,sol.y(44,:),'-.b','linewidth',1);
xlabel(' Time(sec)');
legend('$e_{1,2}$','$e_{2,2}$','$e_{3,2}$','$e_{4,2}$')
set(gca,'FontSize',10,'Fontname', 'Times New Roman');

figure(5)
plot(sol.x,sol.y(6,:),'r',sol.x,sol.y(21,:),'-.r',sol.x,sol.y(33,:),'b',sol.x,sol.y(45,:),'-.b','linewidth',1);
xlabel(' Time(sec)');
legend('$e_{1,3}$','$e_{2,3}$','$e_{3,3}$','$e_{4,3}$')
set(gca,'FontSize',10,'Fontname', 'Times New Roman');


figure(6)
plot(sol.x,sol.y(10,:),'r',sol.x,sol.y(25,:),'-.r',sol.x,sol.y(37,:),'b',sol.x,sol.y(49,:),'-.b','linewidth',1);
xlabel(' Time(sec)');
legend('$\hat{\varphi}_{1,2}$','$\hat{\varphi}_{2,2}$','$\hat{\varphi}_{3,2}$','$\hat{\varphi}_{4,2}$')
set(gca,'FontSize',10,'Fontname', 'Times New Roman');

figure(7)
plot(sol.x,sol.y(52,:),'r',sol.x,sol.y(53,:),'-.r',sol.x,sol.y(54,:),'b',sol.x,sol.y(55,:),'-.b','linewidth',1);
xlabel(' Time(sec)');
legend('$\hat{\varphi}_{1,3}}$','$\hat{\varphi}_{2,3}$','$\hat{\varphi}_{3,3}$','$\hat{\varphi}_{4,3}$')
set(gca,'FontSize',10,'Fontname', 'Times New Roman');

figure(8)
plot(sol.x,sol.yp(12,:),'r',sol.x,sol.yp(26,:),'-.r',sol.x,sol.yp(38,:),'-b',sol.x,sol.yp(50,:),'-.b','linewidth',1);
xlabel(' Time(sec)');
legend('$u_1$','$u_2$','$u_3$','$u_4$')
set(gca,'FontSize',10,'Fontname', 'Times New Roman');

function dydt=tys(t,y,Z,lam,l,d,k)
ylag1 = Z(:,1);
ylag2 = Z(:,2);
ylag3 = Z(:,3);
ylag4 = Z(:,4);
%
lam1=lam(1);
lam2=lam(2);
lam3=lam(3);

l1=l(1);
l2=l(2);
l3=l(3);
d1=d(1);
d2=d(2);
d3=d(3);
k1=k(1);
k2=k(2);
k3=k(3);
yd=sin(t);ydd=cos(t);


r1=1;T=5;q=0.05;
DL1=(r1*(sin((pi/(2*T))*(T-t)))^4)*(t>=0&t<=T)+0*(t>T);
D1=DL1+q;
dotD=(-((2*r1*pi)/T)*cos((pi/(2*T))*(T-t))*(sin((pi/(2*T))*(T-t)))^3)*(t>=0&t<=T)+0*(t>T);

center=[-3;-2;-1;0;1;2;3];

h=[y(1);y(2);y(7);yd;ydd;D1];
S11=FuzzySet(h,center);

h=[y(1);y(2);y(3);y(7);y(8);y(10)];
S12=FuzzySet(h,center);

h=[y(16);y(17);y(22);D1];
S21=FuzzySet(h,center);

h=[y(16);y(17);y(18);y(22);y(23);y(25)];
S22=FuzzySet(h,center);

h=[y(28);y(29);y(34);D1];
S31=FuzzySet(h,center);

h=[y(28);y(29);y(30);y(34);y(35);y(37)];
S32=FuzzySet(h,center);

h=[y(30);y(41);y(46);D1];
S41=FuzzySet(h,center);

h=[y(40);y(41);y(42);y(46);y(47);y(49)];
S42=FuzzySet(h,center);

M=(0.001625/0.9)+((0.506*0.305*0.305)/(3*0.9))+((0.434*0.305*0.305)/(1*0.9))+((2*0.434*0.023*0.023)/(5*0.9));
N=((0.506*0.305*9.8)/(2*0.9))+((0.434*0.305*9.8)/0.9);
L=0.025;
B=0.01625/0.9;K=0.9;R=0.5;
z11=y(7)-yd;
f1=-(2*D1*z11^3)/((D1*D1-z11*z11)^2);
F1=dotD*f1;
omga11=(D1*D1*(D1*D1+z11*z11))/((D1*D1-z11*z11)^2);
elta11=(D1*D1*z11)/((D1-z11)*(D1+z11)); 
alph11=-F1/omga11+ydd-((k1*elta11)/omga11)-((elta11*omga11)/4)+d3;
alph12=-(k2+0.5)*(y(2)-alph11)-(((S11'*S11)*y(10)*(y(2)-alph11))/(2*d1*d1));
u1=-(0.1*k3+0.5)*(y(3)-alph12)-(((S12'*S12)*y(52)*(y(3)-alph12))/(2*d2*d2));

z21=y(22)-y(7);
f2=-(2*D1*z21^3)/((D1*D1-z21*z21)^2);
F2=dotD*f2;
omga21=(D1*D1*(D1*D1+z21*z21))/((D1*D1-z21*z21)^2);
elta21=(D1*D1*z21)/((D1-z21)*(D1+z21)); 
alph21=-F2/omga21+ydd-((k1*elta21)/omga21)-((elta21*omga21)/4)+d3;
alph22=-(k2+0.5)*(y(17)-alph21)-(((S21'*S21)*y(25)*(y(17)-alph21))/(2*d1*d1));
u2=-(0.1*k3+0.5)*(y(18)-alph22)-(((S22'*S22)*y(53)*(y(18)-alph22))/(2*d2*d2));

z31=y(34)-y(7);
f3=-(2*D1*z31^3)/((D1*D1-z31*z31)^2);
F3=dotD*f3;
omga31=(D1*D1*(D1*D1+z31*z31))/((D1*D1-z31*z31)^2);
elta31=(D1*D1*z31)/((D1-z31)*(D1+z31)); 
alph31=-F3/omga31+ydd-((k1*elta31)/omga31)-((elta31*omga31)/4)+d3;
alph32=-(k2+0.5)*(y(29)-alph31)-(((S31'*S31)*y(37)*(y(29)-alph31))/(2*d1*d1));
u3=-(0.1*k3+0.5)*(y(30)-alph32)-(((S32'*S32)*y(54)*(y(30)-alph32))/(2*d2*d2));

z41=y(46)-y(34)+y(46)-y(22);
f4=-(2*D1*z41^3)/((D1*D1-z41*z41)^2);
F4=dotD*f4;
omga41=(D1*D1*(D1*D1+z41*z41))/((D1*D1-z41*z41)^2);
elta41=(D1*D1*z41)/((D1-z41)*(D1+z41)); 
alph41=-F4/omga41+ydd-((k1*elta41)/omga41)-((elta41*omga41)/4)+d3;
alph42=-(k2+0.5)*(y(41)-alph41)-(((S41'*S41)*y(49)*(y(41)-alph41))/(2*d1*d1));
u4=-(0.1*k3+0.5)*(y(42)-alph42)-(((S42'*S42)*y(55)*(y(42)-alph42))/(2*d2*d2));


obs11=y(2)-l1*y(4);
obs12=y(3)/M-l2*y(4)-B*y(2)/M;
obs13=u1/L-l3*y(4)-K*y(2)/L-R*y(3)/L;

obs21=y(17)-l1*y(19);
obs22=y(18)/M-l2*y(19)-B*y(17)/M;
obs23=u2/L-l3*y(19)-K*y(17)/L-R*y(18)/L;

obs31=y(29)-l1*y(31);
obs32=y(30)/M-l2*y(31)-B*y(29)/M;
obs33=u3/L-l3*y(31)-K*y(29)/L-R*y(30)/L;

obs41=y(41)-l1*y(43);
obs42=y(42)/M-l2*y(43)-B*y(41)/M;
obs43=u4/L-l3*y(43)-K*y(41)/L-R*y(42)/L;

obe11=y(5)+l1*y(4);
obe12=y(6)/M+l2*y(4)-(N/M)*sin(y(7))-(B/M)*y(5)+(B/M)*cos((y(2)+y(5)));%*sin((y(3)+y(6)));%*sin((y(1)+y(4)));(y(1)+y(4))*(y(1)+y(4))*((y(2)+y(5))^3);
obe13=l3*y(4)-(K/L)*y(5)-(R/L)*y(6)+(R/L)*sin(y(7));

obe21=y(20)+l1*y(19);
obe22=y(21)/M+l2*y(19)-(N/M)*sin(y(22))-(B/M)*y(20)+(B/M)*cos((y(17)+y(20)));%*sin((y(3)+y(6)));%*sin((y(1)+y(4)));(y(1)+y(4))*(y(1)+y(4))*((y(2)+y(5))^3);
obe23=l3*y(19)-(K/L)*y(20)-(R/L)*y(21)+(R/L)*sin(y(22));

obe31=y(32)+l1*y(31);
obe32=y(33)/M+l2*y(31)-(N/M)*sin(y(34))-(B/M)*y(32)+(B/M)*cos((y(29)+y(32)));%*sin((y(3)+y(6)));%*sin((y(1)+y(4)));(y(1)+y(4))*(y(1)+y(4))*((y(2)+y(5))^3);
obe33=l3*y(31)-(K/L)*y(32)-(R/L)*y(33)+(R/L)*sin(y(34));

obe41=y(44)+l1*y(43);
obe42=y(45)/M+l2*y(43)-(N/M)*sin(y(46))-(B/M)*y(44)+(B/M)*cos((y(41)+y(44)));%*sin((y(3)+y(6)));%*sin((y(1)+y(4)));(y(1)+y(4))*(y(1)+y(4))*((y(2)+y(5))^3);
obe43=l3*y(43)-(K/L)*y(44)-(R/L)*y(45)+(R/L)*sin(y(46));

g11=y(2)+y(5);
g12=y(2)+y(5);
g13=y(3)+y(6);

g21=y(17)+y(20);
g22=y(17)+y(20);
g23=y(18)+y(21);

g31=y(29)+y(32);
g32=y(29)+y(32);
g33=y(30)+y(33);

g41=y(41)+y(44);
g42=y(41)+y(44);
g43=y(42)+y(45);

law11=((lam2*(y(2)-alph11)*(y(2)-alph11)*(S11'*S11))/(2*d1*d1*1))-lam1*y(10);
law12=((lam2*(y(17)-alph21)*(y(17)-alph21)*(S21'*S21))/(2*d1*d1*1))-lam1*y(25);
law13=((lam2*(y(29)-alph31)*(y(29)-alph31)*(S31'*S31))/(2*d1*d1*1))-lam1*y(37);
law14=((lam2*(y(41)-alph41)*(y(41)-alph41)*(S41'*S41))/(2*d1*d1*1))-lam1*y(49);

law21=((lam3*(y(3)-alph12)*(y(3)-alph12)*(S12'*S12))/(2*d2*d2*1))-lam1*y(52);
law22=((lam3*(y(18)-alph22)*(y(18)-alph22)*(S22'*S22))/(2*d2*d2*1))-lam1*y(53);
law23=((lam3*(y(30)-alph32)*(y(30)-alph32)*(S32'*S32))/(2*d2*d2*1))-lam1*y(54);
law24=((lam3*(y(42)-alph42)*(y(42)-alph42)*(S42'*S42))/(2*d2*d2*1))-lam1*y(55);



 dydt = [obs11;obs12;obs13;obe11;obe12;obe13;g11;g12;g13;law11;yd;u1;z11;-D1;D1;obs21;obs22;obs23;obe21;obe22;obe23;g21;g22;g23;law12;u2;z21;obs31;obs32;obs33;obe31;obe32;obe33;g31;g32;g33;law13;u3;z31;obs41;obs42;obs43;obe41;obe42;obe43;g41;g42;g43;law14;u4;z41;law21;law22;law23;law24];

function S = tyshist(t,lam,l,d,k)

 S = [0.1;0.1;0.1;0.05;0.05;0.05;0.15;0.15;0.15;0;0;0;0.15;-1.05;1.05;0.1;0.1;0.1;0.05;0.05;0.05;0.15;0.15;0.15;0;0;0;0.1;0.1;0.1;0.05;0.05;0.05;0.15;0.15;0.15;0;0;0;0.1;0.1;0.1;0.05;0.05;0.05;0.15;0.15;0.15;0;0;0;0;0;0;0];
