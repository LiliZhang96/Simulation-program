function [dx,x_11,x_12,x_21,x_22,x_31,x_32,x_41,x_42,s_11,s_21,s_31,s_41,u1,u2,u3,u4]=comcontrol1(t,x)
x_11=x(1);
x_12=x(2);
x_21=x(3);
x_22=x(4);
x_31=x(5);
x_32=x(6);
x_41=x(7);
x_42=x(8);
theta_11=x(9);
theta_21=x(10);
theta_31=x(11);
theta_41=x(12);
theta_12=x(13);
theta_22=x(14);
theta_32=x(15);
theta_42=x(16);
theta_j12=x(17);
theta_j21=x(18);
theta_j32=x(19);
theta_j43=x(20);
alphabar_12=x(21);
alphabar_22=x(22);
alphabar_32=x(23);
alphabar_42=x(24);
u_1=x(25);
u_2=x(26);
u_3=x(27);
u_4=x(28);

a12=1;a32=1;a43=1;a20=1;d_1=1; d_2=0;d_3=1;d_4=1;

c_11=5;c_21=5;c_31=5;c_41=5;
c_12=5;c_22=5;c_32=5;c_42=5;
tau_1=0.5;tau_2=0.6;tau_3=0.4;tau_4=0.4;
k_11=4;k_12=10;k_21=10;k_22=30;k_31=10;k_32=10;k_41=2.85;k_42=13;
k1_b1=0.5; gamma_11=100;gamma_12=100;gamma_21=5;gamma_22=5;gamma_31=100;gamma_32=100;gamma_41=100;gamma_42=100;
delta_11=0.3;delta_21=0.3;delta_31=25;delta_41=25;
delta_12=0.3;delta_22=0.3;delta_32=10;delta_42=10;
eta_11=2;eta_21=2;eta_31=2;eta_42=2;

%%%%%%%%%%%%%%ding%%%%%%%%%%%%%%%
lambda_12=0.1;lambda_22=0.1;lambda_32=0.1;lambda_42=0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y_d=sin(t);
y_d_dot=cos(t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f11=x_11/(1+x_11^4);
f12=0.1*sin(x_11-x_12)*exp(-x_11^2-x_12^4);
f21=x_21/(1+x_21^4);
f22=0.15*sin(x_21-x_22)*exp(-x_21^2-x_22^4);
f31=x_31/(1+x_31^4);
f32=0.15*sin(x_31-x_32)*exp(-x_31^2-x_32^4);
f41=x_41/(1+x_41^4);
f42=0.2*sin(x_41-x_42)*exp(-x_41^2-x_42^4);




D111=exp(-0.5*(x_11+5)^2);
D112=exp(-0.5*(x_11+2)^2);
D113=exp(-0.5*(x_11+0)^2);
D114=exp(-0.5*(x_11-2)^2);
D115=exp(-0.5*(x_11-5)^2);
D11=D111+D112+D113+D114+D115;
R11=[D111/D11;D112/D11;D113/D11;D114/D11;D115/D11];
S11=norm(R11);

D121=exp(-0.5*(x_11+5)^2)*exp(-0.5*(x_12+5)^2);
D122=exp(-0.5*(x_11+2)^2)*exp(-0.5*(x_12+2)^2);
D123=exp(-0.5*(x_11+0)^2)*exp(-0.5*(x_12+0)^2);
D124=exp(-0.5*(x_11-2)^2)*exp(-0.5*(x_12-2)^2);
D125=exp(-0.5*(x_11-5)^2)*exp(-0.5*(x_12-5)^2);
D12=D121+D122+D123+D124+D125;
R12=[D121/D12;D122/D12;D123/D12;D124/D12;D125/D12];
S12=norm(R12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D211=exp(-0.5*(x_21+5)^2);
D212=exp(-0.5*(x_21+2)^2);
D213=exp(-0.5*(x_21+0)^2);
D214=exp(-0.5*(x_21-2)^2);
D215=exp(-0.5*(x_21-5)^2);
D21=D211+D212+D213+D214+D215;
R21=[D211/D21;D212/D21;D213/D21;D214/D21;D215/D21];
S21=norm(R21);

D221=exp(-0.5*(x_21+5)^2)*exp(-0.5*(x_22+5)^2);
D222=exp(-0.5*(x_21+2)^2)*exp(-0.5*(x_22+2)^2);
D223=exp(-0.5*(x_21+0)^2)*exp(-0.5*(x_22+0)^2);
D224=exp(-0.5*(x_21-2)^2)*exp(-0.5*(x_22-2)^2);
D225=exp(-0.5*(x_21-5)^2)*exp(-0.5*(x_22-5)^2);
D22=D221+D222+D223+D224+D225;
R22=[D221/D22;D222/D22;D223/D22;D224/D22;D225/D22];
S22=norm(R22);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


D311=exp(-0.5*(x_31+5)^2);
D312=exp(-0.5*(x_31+2)^2);
D313=exp(-0.5*(x_31+0)^2);
D314=exp(-0.5*(x_31-2)^2);
D315=exp(-0.5*(x_31-5)^2);
D31=D311+D312+D313+D314+D315;
R31=[D311/D31;D312/D31;D313/D31;D314/D31;D315/D31];
S31=norm(R31);

D321=exp(-0.5*(x_31+5)^2)*exp(-0.5*(x_32+5)^2);
D322=exp(-0.5*(x_31+2)^2)*exp(-0.5*(x_32+2)^2);
D323=exp(-0.5*(x_31+0)^2)*exp(-0.5*(x_32+0)^2);
D324=exp(-0.5*(x_31-2)^2)*exp(-0.5*(x_32-2)^2);
D325=exp(-0.5*(x_31-5)^2)*exp(-0.5*(x_32-5)^2);
D32=D321+D322+D323+D324+D325;
R32=[D321/D32;D322/D32;D323/D32;D324/D32;D325/D32];
S32=norm(R32);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D411=exp(-0.5*(x_41+5)^2);
D412=exp(-0.5*(x_41+2)^2);
D413=exp(-0.5*(x_41+0)^2);
D414=exp(-0.5*(x_41-2)^2);
D415=exp(-0.5*(x_41-5)^2);
D41=D411+D412+D413+D414+D415;
R41=[D411/D41;D412/D41;D413/D41;D414/D41;D415/D41];
S41=norm(R41);

D421=exp(-0.5*(x_41+5)^2)*exp(-0.5*(x_42+5)^2);
D422=exp(-0.5*(x_41+2)^2)*exp(-0.5*(x_42+2)^2);
D423=exp(-0.5*(x_41+0)^2)*exp(-0.5*(x_42+0)^2);
D424=exp(-0.5*(x_41-2)^2)*exp(-0.5*(x_42-2)^2);
D425=exp(-0.5*(x_41-5)^2)*exp(-0.5*(x_42-5)^2);
D42=D421+D422+D423+D424+D425;
R42=[D421/D42;D422/D42;D423/D42;D424/D42;D425/D42];
S42=norm(R42);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s_11=(x_11-x_21);
s_21=(x_21-y_d);
s_31=(x_31-x_21);
s_41=(x_41-x_31);
s_12=x_12-alphabar_12;
s_22=x_22-alphabar_22;
s_32=x_32-alphabar_32;
s_42=x_42-alphabar_42;

% alpha_12=-c_11*s_11-2*s_11/(k1_b1^2-s_11^2)+a12*(x_22+theta_j12*S21)-theta_11*S11;
% alpha_22=-c_21*s_21+y_d_dot-3*s_21/(k1_b1^2-s_21^2)-theta_21*S21;
% alpha_32=-c_31*s_31-2*s_31/(k1_b1^2-s_31^2)+a32*(x_32+theta_j32*S21)-theta_31*S31;
% alpha_42=-c_41*s_41-2*s_41/(k1_b1^2-s_41^2)+a43*(x_42+theta_j43*S31)-theta_41*S41;
% 
alpha_12=-c_11*s_11-2*s_11/(k1_b1^2-s_11^2)+a12*(x_22+f21)-f11;
alpha_22=-c_21*s_21+y_d_dot-3*s_21/(k1_b1^2-s_21^2)-f21;
alpha_32=-c_31*s_31-2*s_31/(k1_b1^2-s_31^2)+a32*(x_32+f21)-f31;
alpha_42=-c_41*s_41-2*s_41/(k1_b1^2-s_41^2)+a43*(x_42+f31)-f41;

% u1=-c_12*s_12-s_12-theta_12*S12+(alpha_12-alphabar_12)/lambda_12;
% u2=-c_22*s_22-s_22-theta_22*S22+(alpha_22-alphabar_22)/lambda_22;
% u3=-c_32*s_32-s_32-theta_32*S32+(alpha_32-alphabar_32)/lambda_32;
% u4=-c_42*s_42-s_42-theta_42*S42+(alpha_42-alphabar_42)/lambda_42;

u1=-c_12*s_12-s_12-f12+(alpha_12-alphabar_12)/lambda_12;
u2=-c_22*s_22-s_22-f22+(alpha_22-alphabar_22)/lambda_22;
u3=-c_32*s_32-s_32-f32+(alpha_32-alphabar_32)/lambda_32;
u4=-c_42*s_42-s_42-f42+(alpha_42-alphabar_42)/lambda_42;
% 


x_11_dot=x_12+f11;
x_12_dot=0.3*u1+f12; 
x_21_dot=x_22+f21;
x_22_dot=0.4*u2+f22; 
x_31_dot=x_32+f31;
x_32_dot=0.5*u3+f32; 
x_41_dot=x_42+f41;
x_42_dot=0.5*u4+f42; 

theta_dot11=gamma_11*(s_11/(k1_b1^2-s_11^2)*S11-delta_11*theta_11);
theta_dot21=gamma_21*(s_21/(k1_b1^2-s_21^2)*S21-delta_21*theta_21);
theta_dot31=gamma_31*(s_31/(k1_b1^2-s_31^2)*S31-delta_31*theta_31);
theta_dot41=gamma_41*(s_41/(k1_b1^2-s_41^2)*S41-delta_41*theta_41);

theta_dot12=gamma_12*(s_12*S12-delta_12*theta_12);
theta_dot22=gamma_22*(s_22*S22-delta_22*theta_22);
theta_dot32=gamma_32*(s_32*S32-delta_32*theta_32);
theta_dot42=gamma_42*(s_42*S42-delta_42*theta_42);

theta_dotj12=eta_11*(-s_11/(k1_b1^2-s_11^2)*S21-delta_11*theta_j12);
theta_dotj21=0;
theta_dotj32=eta_21*(-s_31/(k1_b1^2-s_31^2)*S21-delta_21*theta_j32);
theta_dotj43=eta_31*(-s_41/(k1_b1^2-s_41^2)*S31-delta_31*theta_j43);

alphabar_12dot=(alpha_12-alphabar_12)/lambda_12;
alphabar_22dot=(alpha_22-alphabar_22)/lambda_22;
alphabar_32dot=(alpha_32-alphabar_32)/lambda_32;
alphabar_42dot=(alpha_42-alphabar_42)/lambda_42;

theta_j21=0;
    tt=t
dx=[x_11_dot;x_12_dot;x_21_dot;x_22_dot;x_31_dot;x_32_dot;x_41_dot;x_42_dot;theta_dot11;theta_dot21;theta_dot31;theta_dot41;theta_dot12;theta_dot22;theta_dot32;theta_dot42;theta_dotj12;theta_dotj21;theta_dotj32;theta_dotj43;alphabar_12dot;alphabar_22dot;alphabar_32dot;alphabar_42dot;u1;u2;u3;u4];
