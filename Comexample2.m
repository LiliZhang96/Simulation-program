function tys

lam=[1;5;5];

l=-[8.6485;336.2915;-290.6175];
d=[0.5;0.5;0]*2;
k=[5;5;5];
Tf=[0:0.001:30];
L1=length(Tf);

lags=[1,1,1,1]*2.5; 


sol=dde23(@tys,lags,@tyshist,[0,30],[],lam,l,d,k);


figure(2)
subplot(2,1,2);
plot(sol.x,sol.yp(13,:),'r',sol.x,sol.yp(27,:),'-.r',sol.x,sol.yp(39,:),'b',sol.x,sol.yp(51,:),'-.b',sol.x,sol.yp(14,:),'-.k',sol.x,sol.yp(15,:),'-.k','linewidth',1);
xlabel(' Time(sec)');
title('(b) Under the controller designed in [59]')
legend('$h_{1,1}$','$h_{2,1}$','$h_{3,1}$','$h_{4,1}$','0.05','-0.05')
set(gca,'FontSize',10,'Fontname', 'Times New Roman');

function dydt=tys(t,y,Z,lam,l,d,k)
y1=0.05;
y2=-0.05;
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

center=[-3;-2;-1;0;1;2;3];

h=[y(1);y(2);y(7);yd;ydd];
S11=FuzzySet(h,center);

h=[y(1);y(2);y(3);y(7);y(8);y(10)];
S12=FuzzySet(h,center);

h=[y(16);y(17);y(22)];
S21=FuzzySet(h,center);

h=[y(16);y(17);y(18);y(22);y(23);y(25)];
S22=FuzzySet(h,center);

h=[y(28);y(29);y(34)];
S31=FuzzySet(h,center);

h=[y(28);y(29);y(30);y(34);y(35);y(37)];
S32=FuzzySet(h,center);

h=[y(30);y(41);y(46)];
S41=FuzzySet(h,center);

h=[y(40);y(41);y(42);y(46);y(47);y(49)];
S42=FuzzySet(h,center);

M=(0.001625/0.9)+((0.506*0.305*0.305)/(3*0.9))+((0.434*0.305*0.305)/(1*0.9))+((2*0.434*0.023*0.023)/(5*0.9));
N=((0.506*0.305*9.8)/(2*0.9))+((0.434*0.305*9.8)/0.9);
L=0.025;
B=0.01625/0.9;K=0.9;R=0.5;
z11=y(7)-yd;
alph11=-k1*z11+ydd+d3;
alph12=-(k2+0.5)*(y(2)-alph11)-(((S11'*S11)*y(10)*(y(2)-alph11))/(2*d1*d1));
u1=-(0.1*k3+0.5)*(y(3)-alph12)-(((S12'*S12)*y(52)*(y(3)-alph12))/(2*d2*d2));

z21=y(22)-y(7);
alph21=-k1*z21+ydd+d3;
alph22=-(k2+0.5)*(y(17)-alph21)-(((S21'*S21)*y(25)*(y(17)-alph21))/(2*d1*d1));
u2=-(0.1*k3+0.5)*(y(18)-alph22)-(((S22'*S22)*y(53)*(y(18)-alph22))/(2*d2*d2));

z31=y(34)-y(7);
alph31=-k1*z31+ydd+d3;
alph32=-(k2+0.5)*(y(29)-alph31)-(((S31'*S31)*y(37)*(y(29)-alph31))/(2*d1*d1));
u3=-(0.1*k3+0.5)*(y(30)-alph32)-(((S32'*S32)*y(54)*(y(30)-alph32))/(2*d2*d2));

z41=y(46)-y(34)+y(46)-y(22);
alph41=-k1*z41+ydd+d3;
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



 dydt = [obs11;obs12;obs13;obe11;obe12;obe13;g11;g12;g13;law11;yd;u1;z11;y1;y2;obs21;obs22;obs23;obe21;obe22;obe23;g21;g22;g23;law12;u2;z21;obs31;obs32;obs33;obe31;obe32;obe33;g31;g32;g33;law13;u3;z31;obs41;obs42;obs43;obe41;obe42;obe43;g41;g42;g43;law14;u4;z41;law21;law22;law23;law24];

function S = tyshist(t,lam,l,d,k)

 S = [0.1;0.1;0.1;0.05;0.05;0.05;0.15;0.15;0.15;0;0;0;0.15;0.05;0.05;0.1;0.1;0.1;0.05;0.05;0.05;0.15;0.15;0.15;0;0;0;0.1;0.1;0.1;0.05;0.05;0.05;0.15;0.15;0.15;0;0;0;0.1;0.1;0.1;0.05;0.05;0.05;0.15;0.15;0.15;0;0;0;0;0;0;0];
