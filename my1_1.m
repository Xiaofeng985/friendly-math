set(0, 'defaultLegendInterpreter', 'latex');
set(0, 'defaultLegendFontName', 'Times New Roman');

dt=0.001;
t=0:dt:10;
t1=0:dt:9.999;
N=length(t);

%设计参数
r=3/5;
%控制器参数
k11=1;k12=1;k13=1;k21=2;k22=2;k23=2;
s11=2;s12=2;s13=2;s21=2;s22=2;s23=2;
%大规模参数
psi1=1;psi2=1;
%补偿信号参数
l11=5;l12=5;l13=5;l21=5;l22=5;l23=5;
%滤波器参数
r111=500;r112=500;r121=500;r122=500;
r211=500;r212=500;r221=500;r222=500;

%参考信号
x11c=@(t)(sin(t));
x21c=@(t)(0.5*sin(t)+0.5*sin(0.5*t));
%定义函数
H111=@(t)(t^4);H112=@(t)(t^4);H211=@(t)(t^2);

%系统初值
x11_0=0.3;x12_0=0.1;x13_0=0.1;
x21_0=0.3;x22_0=0.1;x23_0=0.1;
%滤波器初值
phi111_0=0.2;phi112_0=0.2;phi121_0=0.2;phi122_0=0.2;
phi211_0=0.2;phi212_0=0.2;phi221_0=0.2;phi222_0=0.2;
%补偿信号初值
xi11_0=0;xi12_0=0;xi13_0=0;xi21_0=0;xi22_0=0;xi23_0=0;
%自适应滤初值
rho1_0=0.6;rho2_0=0.6;

%初始化数组
x11=zeros(N,1);x12=zeros(N,1);x13=zeros(N,1);x21=zeros(N,1);x22=zeros(N,1);x23=zeros(N,1);
phi111=zeros(N,1);phi112=zeros(N,1);phi121=zeros(N,1);phi122=zeros(N,1);
phi211=zeros(N,1);phi212=zeros(N,1);phi221=zeros(N,1);phi222=zeros(N,1);
xi11=zeros(N,1);xi12=zeros(N,1);xi13=zeros(N,1);xi21=zeros(N,1);xi22=zeros(N,1);xi23=zeros(N,1);
rho1=zeros(N,1);rho2=zeros(N,1);
%赋值
x11(1)=x11_0;x12(1)=x12_0;x13(1)=x13_0;x21(1)=x21_0;x22(1)=x22_0;x23(1)=x23_0;
phi111(1)=phi111_0;phi112(1)=phi112_0;phi121(1)=phi121_0;phi122(1)=phi122_0;
phi211(1)=phi211_0;phi212(1)=phi212_0;phi221(1)=phi221_0;phi222(1)=phi222_0;
xi11(1)=xi11_0;xi12(1)=xi12_0;xi13(1)=xi13_0;
xi21(1)=xi21_0;xi22(1)=xi22_0;xi23(1)=xi23_0;
rho1(1)=rho1_0;rho2(1)=rho2_0;

for i=1:N-1
    w1(i)=(2*(x11(i)-x11c(t(i))-xi11(i)))*(H111(x11(i))+H211(x11(i)))/((x11(i)-x11c(t(i))-xi11(i))*(x11(i)-x11c(t(i))-xi11(i))+psi1);
    w2(i)=(2*(x21(i)-x21c(t(i))-xi21(i)))*(H112(x21(i)))/((x21(i)-x21c(t(i))-xi21(i))*(x21(i)-x21c(t(i))-xi21(i))+psi2);
    a11(i)=-k11*(x11(i)-x11c(t(i)))+cos(t(i))-s11*(x11(i)-x11c(t(i))-xi11(i))^r-rho1(i)*w1(i);
    v111(i)=-r111*((abs(phi111(i)-a11(i)))^(1/2))*sign(phi111(i)-a11(i))+phi112(i);
    phi111(i+1)=phi111(i)+v111(i)*dt;
    phi112(i+1)=phi112(i)+(-r112*sign(phi112(i)-v111(i)))*dt;
    a12(i)=-k12*(x12(i)-phi111(i))+v111(i)-(-sin(x12(i))-x11(i))-(x11(i)-x11c(t(i)))-s12*(x12(i)-phi111(i)-xi12(i))^r;
    v121(i)=-r121*((abs(phi121(i)-a12(i)))^(1/2))*sign(phi121(i)-a12(i))+phi122(i);
    phi121(i+1)=phi121(i)+v121(i)*dt;
    phi122(i+1)=phi122(i)+(-r122*sign(phi122(i)-v121(i)))*dt;
    xi11(i+1)=xi11(i)+(-k11*xi11(i)+(phi111(i)-a11(i))+xi12(i)-l11*sign(xi11(i)))*dt;
    xi12(i+1)=xi12(i)+(-k12*xi12(i)+(phi121(i)-a12(i))-xi11(i)+xi13(i)-l12*sign(xi12(i)))*dt;
    xi13(i+1)=xi13(i)+(-k13*xi13(i)-xi12(i)-l13*sign(xi13(i)))*dt;
    rho1(i+1)=rho1(i)+((x11(i)-x11c(t(i))-xi11(i))*w1(i)-2*rho1(i))*dt;
    u1(i)=-k13*(x13(i)-phi121(i))+v121(i)-(exp(-x12(i))*x11(i)+x13(i))-(x12(i)-phi111(i))-s13*(x13(i)-phi121(i)-xi13(i))^r;
    
    a21(i)=-k21*(x21(i)-x21c(t(i)))+0.5*cos(t(i))+0.25*cos(0.5*t(i))-s21*(x21(i)-x21c(t(i))-xi21(i))^r-rho2(i)*w2(i);
    v211(i)=-r211*((abs(phi211(i)-a21(i)))^(1/2))*sign(phi211(i)-a21(i))+phi212(i);
    phi211(i+1)=phi211(i)+v211(i)*dt;
    phi212(i+1)=phi212(i)+(-r212*sign(phi212(i)-v211(i)))*dt;
    a22(i)=-k22*(x22(i)-phi211(i))+v211(i)-(x21(i)*sin(x22(i)))-(x21(i)-x21c(t(i)))-s22*(x22(i)-phi211(i)-xi22(i))^r;
    v221(i)=-r221*((abs(phi221(i)-a22(i)))^(1/2))*sign(phi221(i)-a22(i))+phi222(i);
    phi221(i+1)=phi221(i)+v221(i)*dt;
    phi222(i+1)=phi222(i)+(-r222*sign(phi222(i)-v221(i)))*dt;
    xi21(i+1)=xi21(i)+(-k21*xi21(i)+(phi211(i)-a21(i))+xi22(i)-l21*sign(xi21(i)))*dt;
    xi22(i+1)=xi22(i)+(-k22*xi22(i)+phi221(i)-a22(i)-xi21(i)+xi23(i)-l22*sign(xi22(i)))*dt;
    xi23(i+1)=xi23(i)+(-k23*xi23(i)-xi22(i)-l23*sign(xi23(i)))*dt;
    rho2(i+1)=rho2(i)+((x21(i)-x21c(t(i))-xi21(i))*w2(i)-2*rho2(i))*dt;
    u2(i)=-k23*(x23(i)-phi221(i))+v221(i)-(-x22(i)-x23(i))-(x22(i)-phi211(i))-s23*(x23(i)-phi221(i)-xi23(i))^r;
    
    x11(i+1)=x11(i)+x12(i)*dt;
    x12(i+1)=x12(i)+(x13(i)-sin(x12(i))-x11(i)+x21(i)*x11(i))*dt;
    x13(i+1)=x13(i)+(u1(i)+exp(-x12(i))*x11(i)+x13(i))*dt;
    x21(i+1)=x21(i)+x22(i)*dt;
    x22(i+1)=x22(i)+(x23(i)+x21(i)*sin(x22(i))+sin(x11(i)))*dt;
    x23(i+1)=x23(i)+(u2(i)-x22(i)-x23(i))*dt;    
end
figure;
plot(t,sin(t),'--b','Linewidth',2,'DisplayName','$sin(t)$');
hold on;
plot(t,x11,'-r','Linewidth',2,'DisplayName','$x_{11}$');
xlabel('Time (s)','Fontname','Times New Roman')
ylabel('Tracking performance','Fontname','Times New Roman')
legend;
figure;
plot(t,0.5*sin(t)+0.5*sin(0.5*t),'--b','Linewidth',2,'DisplayName','$0.5sin(t)+0.5sin(0.5t)$');
hold on;
plot(t,x21,'-r','Linewidth',2,'DisplayName','$x_{21}$');
xlabel('Time (s)','Fontname','Times New Roman')
ylabel('Tracking performance','Fontname','Times New Roman')
legend;
figure;
plot(t1,u1,'DisplayName','$u_1$');
xlabel('Time (s)','Fontname','Times New Roman')
legend;
figure;
plot(t1,u2,'DisplayName','$u_2$');
xlabel('Time (s)','Fontname','Times New Roman')
legend;




