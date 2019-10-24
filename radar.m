%===================================雷达系统第一题==========================%
clc;clear all
%============================参数部分==========================%
c=3e8; %m/s 光速
k=1.38e-23; % 玻尔兹曼常数
T0=290; %K 标准噪声温度
%-------------------------题目---------------------%
Rmax=80000; %m 不模糊探测距离
Lambda=0.03; %m 波长
D=0.25; %m 天线等效孔径（直径）
F=3; %dB 噪声系数
L=4; %dB 系统损耗
sita_3dB=6; %degree 天线波束宽度6°
RCS=1500; %m^2 目标RCS
Vs=20; %m/s 目标航速
Va=600; %m/s 导弹运动速度
alpha=30; %degree 目标航向与弹轴方向之间的夹角30°
beta=1*pi/180; %rad 目标偏离弹轴方向的角度
PRT_search=800e-6; %s 搜索状态的脉冲重复周期
Tp=160e-6; %s 脉冲宽度1us
B=1e6; %Hz 调频带宽1MHz

f0=c/Lambda;
Fs=2*B; % 采样频率
Ts=1/Fs; % 采样周期
mu=B/Tp; % 调频斜率
Ae=pi*(D/2)^2; %m^2 天线有效面积
G=4*pi*Ae/Lambda^2;%天线增益

%% ============================1.模糊函数==========================%
EPS=10^-10; % 防止出现奇异点
tau=-Tp:Tp/1600:Tp-Tp/1600; % tau的离散
fd=-B:B/1000:B-B/1000; % fd的离散
[tau_index,fd_index]=meshgrid(tau,fd);
temp1=1-abs(tau_index)./Tp;
temp2=pi*Tp*(mu*tau_index+fd_index).*temp1+eps;
ambiguity=abs(temp1.*sin(temp2)./temp2); % 模糊函数的模
tau0=find(tau==0);
fd0=find(fd==0);
figure(1)%模糊函数
mesh(tau*1e6,fd*1e-6,ambiguity),grid on;colormap('jet');
xlabel('\bf \tau / us');ylabel('\bf fd / MHz');
title('\bf 模糊函数图');
figure(2) % 速度模糊函数
plot(fd*1e-6,10*log10(ambiguity(:,tau0)),'b','LineWidth',1),grid on
xlabel('\bf fd / MHz');
title('\bf 速度模糊函数图');
axis([-1 1 -60 0]);
figure(3)% 距离模糊函数
plot(tau*1e6,10*log10(ambiguity(fd0,:)),'b','LineWidth',1),grid on
xlabel('\bf \tau / us');
title('\bf 距离模糊函数图');
axis([-160 160 -60 0]);
figure(4)%-4dB等高线图
contour(tau*1e6,fd*1e-6,ambiguity,[10^(-4/20),10^(-4/20)],'bl');
xlabel('\bf \tau / us'),ylabel('\bf fd / MHz');
title('\bf 模糊函数的-4dB切割等高线图');
[fd_index_3dB,tau_index_3dB]=find(abs(20*log10(ambiguity)-(-3)<0.1));
tau_3dB=abs(tau(tau_index_3dB(end))-tau(tau_index_3dB(1)))
B_3dB=abs(fd(fd_index_3dB(end))-fd(fd_index_3dB(1)))

%% ============================2.相干积累==========================%
V_scan=60*pi/180; %rad 天线扫描速度60°/s
N=sita_3dB/V_scan/PRT_search; % 积累脉冲数
Pt=30; %w 发射功率
R_snr=[1/800:1/800:1]*Rmax;
temp3=Pt*G^2*Lambda^2*RCS/( (4*pi)^3*10^(L/10) );
temp4=k*B*T0*10^(F/10); % ?????????????????单位
SNR_1=10*log10(temp3/temp4*1./(R_snr.^4));
SNR_N=SNR_1+10*log10(N);
figure(5);
plot(R_snr*1e-3,SNR_1,'b','LineWidth',2),hold on
plot(R_snr*1e-3,SNR_N,'r','LineWidth',2),grid on
xlabel('\bf R / km'),ylabel('\bf SNR / dB');
title('\bf 相干积累前后信噪比-距离关系曲线');
legend('\bf 相干积累前','\bf 相干积累后');

%-------------------------雷达威力图---------------------%
ae=8494; %km 地球有效半径
phi0=13; %degree 波束中心13°
phi=[0:0.1:90]; %degree 俯仰角离散
f=exp(-1.3863*((phi-phi0)/10).^2); %高斯型方向图
Gt_dB=10*log10(4*pi*Ae/Lambda^2);
Pt=1000; %w 发射的峰值功率
Gt=f*10^(Gt_dB/10); % 发射方向图
Gr_dB=Gt_dB; %dB 接收增益
Gr=f*10^(Gr_dB/10); % 发射方向图

%---等距线---%
R=[0:10:100]; %km 等斜距
xr_range=cos(phi'/180*pi)*R; %km 等距线的x轴
yr_range=sin(phi'/180*pi)*R; %km 等距线的y轴
%---等高线---%
Ht=[5:5:35]; %km 等高度
xr_high=zeros(length(Ht),length([min(Ht):1:100]));
yr_high=zeros(length(Ht),length([min(Ht):1:100]));
for num_high=1:length(Ht)
    r=[Ht(num_high):1:100]; %km 等高线的起始斜距应为高度##############################
    lenr_high(num_high)=length(r);
    phi_high=asin(Ht(num_high)./r-r/ae/2)/pi*180; %rad 仰角
    xr_high(num_high,1:lenr_high(num_high))=cos(phi_high/180*pi).*r; %km 等高线的x轴
    yr_high(num_high,1:1:lenr_high(num_high))=sin(phi_high/180*pi).*r; %km 等高线的y轴
end
%---等仰角线---%
r_ele=[0:1:100]; %km 等仰角线的距离离散
Lenr_ele=length(r_ele);
phi1=[0:1:15]';
phi2=cat(1,[16:1:18]',[20:2:30]',[40:10:90]');
Phi=cat(1,phi1,phi2); %degree 等仰角
xr_ele=cos(Phi/180*pi)*r_ele; %km 等高线的x轴
yr_ele=sin(Phi/180*pi)*r_ele; %km 等高线的y轴
%--------------灵敏度-------------%
k=1.38*10^(-23); % 波尔兹曼常数
T0=290; %K 绝对温度
B=1e6; %MHz 带宽
F=3; %dB 噪声系数
D0=5; %dB 检测因子
simin=10*log10(k*T0*B*10^(F/10)*10^(D0/10));
%
%---威力---%
R4=Pt*Gt.*Gr*Lambda^2*RCS/( (4*pi)^3*10^(simin/10)*10^(L/10) ); % 威力
Rcover=R4.^(1/4)/1000; %km 雷达威力图
xr_cover=cos(phi/180*pi).*Rcover; %km 威力图的x轴
yr_cover=sin(phi/180*pi).*Rcover; %km 威力图的y轴
%画图
figure(6)
for num_range=1:length(R) % 等距线
    plot(xr_range(:,num_range),yr_range(:,num_range),'g--'),hold on
end
for num_high=1:length(Ht) % 等高线
    plot(xr_high(num_high,1:lenr_high(num_high)),yr_high(num_high,1:lenr_high(num_high)),'k--'),hold on
end
num_text=zeros(1,length(Phi));
for num_ele=1:length(Phi) % 等仰角线
    for m=1:Lenr_ele
        if (xr_ele(num_ele,m)<90 && xr_ele(num_ele,m+1)>=90) || (yr_ele(num_ele,m)<24 && yr_ele(num_ele,m+1)>=24)
            num_text(num_ele)=m;
            break;
        end
    end
    plot(xr_ele(num_ele,:),yr_ele(num_ele,:),'b--'),hold on
end
plot(xr_cover,yr_cover,'r','LineWidth',2)
xlabel('\bf 水平距离 / km'),ylabel('\bf 高度 / km')
title('\bf 雷达威力图')
axis([0 85 0 24])

for k=1:length(phi1)
    text((85+0.5),(85+0.5)*tand(phi1(k)),[num2str(phi1(k)) '\circ'])
end
for k=1:length(phi2)
    text((24+0.5)/tand(phi2(k))-1.5,24+0.55,[num2str(phi2(k)) '\circ']);  %*onse(1,length(phi1))
end
%% ============================3.匹配滤波==========================%
FFT_POINTS=2^nextpow2(Tp*Fs);
t_3=[-1/2:1/FFT_POINTS:1/2-1/FFT_POINTS]*Tp; %s 时间离散，2倍带宽采样
f_indice=[-1/2:1/FFT_POINTS:1/2-1/FFT_POINTS]*Fs; % 频率轴
st_lfm=exp(1j*pi*mu*t_3.^2); % LFM信号
%-------------------------(1)ht和Hf---------------------%
%--------------(a)ht和Hf-----------------%
ht_lfm=conj(fliplr(st_lfm)); % 匹配滤波器
Hf_lfm=fftshift(fft(ht_lfm,FFT_POINTS)); % 匹配滤波器的频率响应
figure(7) %匹配滤波器的实部和虚部
plot(t_3*1e6,real(ht_lfm),'b'),hold on
plot(t_3*1e6,imag(ht_lfm),'r'),grid on
xlabel('\bf t / us'),ylabel('\bf h(t)');
title('\bf 匹配滤波器h(t)')
legend('\bf 实部','\bf 虚部')
figure(8) %匹配滤波器的频率响应
plot(f_indice*1e-6,abs(Hf_lfm),'b','LineWidth',2),grid on
xlabel('\bf f / MHz'),ylabel('\bf H(f)');
title('\bf 匹配滤波器的频率响应H(f)')
%-----------------(b)加窗-----------------%
st_lfm_f=fft(st_lfm,FFT_POINTS); % 基带信号的FFT
% 不加窗
Hf_lfm_NoWin=fft(ht_lfm,FFT_POINTS);
sr_lfm_NoWin=abs(ifftshift(ifft(st_lfm_f.*Hf_lfm_NoWin,FFT_POINTS)));
% 加汉明窗
Hf_lfm_Hamming=fft(ht_lfm.*hamming(length(ht_lfm)).',FFT_POINTS);
sr_lfm_Hamming=abs(ifftshift(ifft(st_lfm_f.*Hf_lfm_Hamming,FFT_POINTS)));
% 加泰勒窗
Hf_lfm_Taylor=fft(ht_lfm.*taylorwin(length(ht_lfm)).',FFT_POINTS);
sr_lfm_Taylor=abs(ifftshift(ifft(st_lfm_f.*Hf_lfm_Taylor,FFT_POINTS)));

figure(9)
plot((-Tp/2:Tp/(FFT_POINTS-1):Tp/2)*1e6,20*log10(sr_lfm_NoWin/max(sr_lfm_NoWin)),'k','linewidth',1),hold on
plot((-Tp/2:Tp/(FFT_POINTS-1):Tp/2)*1e6,20*log10(sr_lfm_Hamming/max(sr_lfm_Hamming)),'b','linewidth',1),hold on
plot((-Tp/2:Tp/(FFT_POINTS-1):Tp/2)*1e6,20*log10(sr_lfm_Taylor/max(sr_lfm_Taylor)),'r','linewidth',1),grid on
xlabel('\bf \tau / us'),ylabel('\bf 归一化幅度 / dB')
title('\bf 加窗与不加窗的脉压结果')
legend({'\bf 不加窗','\bf Hamming窗','\bf Taylor窗'})
ylim([-80 0])
%-----------(2)相位编码信号的多普勒敏感性-------------%
Tp_code=1e-6;
Ts_code=0.25e-6;
R0=[40,70]*1000;
vr=[0,0];
SNR=[14,10];
Rmin=30e3;
Rrec=60e3;
bos=2*pi/Lambda;
M=round(Tp_code/Ts_code);
n=2^7-1;
code=idinput(n,'prbs');
code2=kron(code',ones(1,M));
NR0=ceil(log2(2*Rrec/c/Ts_code)); 
NR1=2^NR0; % DFT的点数
M2=M*length(code);
t=(0:M2-1)*Ts_code; % 时间离散
sp=0.707*(randn(1,NR1)+1i*randn(1,NR1)); % 噪声
for k=1:length(R0)
    NR=fix(2*(R0(k)-Rmin)/c/Ts_code);
    Ri=2*(R0(k)-vr(k)*t);
    spt=(10^(SNR(k)/20))*exp(-1i*bos*Ri).*code2; % 信号
    sp(NR:NR+M2-1)=sp(NR:NR+M2-1)+spt;
end
spf=fft(sp,NR1);
Wf_t=fft(code2,NR1);
y=abs(ifft(spf.*conj(Wf_t),NR1))/NR0;
figure(10)% M序列
plot(t*1e6,real(code2)),grid on
xlabel('\bf 时间 / us'),ylabel('\bf 匹配滤波器信号实部')
% title(['码长为' num2str(n) '的M序列' ])
axis([0 128 -1.2 1.2])
figure(11)% M序列的非周期自相关
p=xcorr(code2,'coeff');
plot((1:length(p))-ceil(length(p)/2),20*log10(abs(p))),grid on
xlabel('\bf 时间 / us'),ylabel('\bf 非周期自相关函数 / dB')
% title('\bf M序列的非周期自相关函数')
axis([-120 120 -40 0.2])
figure(12)% 脉压输入信号实部
plot(real(sp)),grid on
xlabel('\bf 时域采样点'),ylabel('\bf 脉压输入信号实部')
% title('\bf 脉压输入信号实部')
xlim([0 2000]);
figure(13)% 速度为0的脉压结果
plot((0:NR1-1)*c*Ts_code/2/1000+Rmin/1000,20*log10(y)),grid on
xlabel('\bf 距离 / km'),ylabel('\bf 脉压输出 / dB')
% title('\bf 脉压结果(目标在[40 70]km处,速度为[0 0])')
axis([Rmin/1000 Rmin/1000+Rrec/1000 -20 50]);

vr=[100,40];
sp=0.707*(randn(1,NR1)+1i*randn(1,NR1));             %噪声
for k=1:length(R0)
    NR=fix(2*(R0(k)-Rmin)/c/Ts_code);              %fix为向下取整；
    Ri=2*(R0(k)-vr(k)*t);
    spt=(10^(SNR(k)/20))*exp(-1i*bos*Ri).*code2;     %信号
    sp(NR:NR+M2-1)=sp(NR:NR+M2-1)+spt;
end
spf=fft(sp,NR1);
Wf_t=fft(code2,NR1);
y=abs(ifft(spf.*conj(Wf_t),NR1))/NR0;
figure(14)% 速度不为0的脉压结果
plot((0:NR1-1)*c*Ts_code/2/1000+Rmin/1000,20*log10(y)),grid on
xlabel('\bf 距离/km'),ylabel('\bf 脉压输出/dB');
title('\bf 脉压结果(目标在[40 70]km处,速度为[100 40]m/s)')
axis([Rmin/1000 Rmin/1000+Rrec/1000 -20 50]);
%% =========================4.搜索状态的信号处理=======================%
SNR_4=-12;
Vr=600+Vs*cosd(30);
Rmin=c/2*Tp;
Rrec=100e3; % 接收窗
M=round(Tp/Ts);
t_4=0:Ts:Tp-Ts;
NR0=ceil(log2(2*Rrec/c/Ts));
NR1=2^NR0;
window=taylorwin(length(t_4));
st_lfm4=exp(1i*pi*mu*t_4.^2).*window';
sp=0.707*(randn(64,NR1)+1i*randn(64,NR1));
sp_noise=sp;
Pn_PPbef=10*log10(mean(real(sp_noise(1,:)).^2)); % 脉压前的噪声功率
R=0:75:(M-1)*c*Ts/2; %在0和Rmin之间采样，采样间隔为75m
NR=fix(2*(Rmax-Rmin)/c/Ts);
Ri=zeros(64,M);
for i=1:M
    for k=1:64
        Ri(k,i)=R(i)-Vr*Ts*(k-1); %Ri的列代表不同的采样位置，行代表同一初始采样位置下速度引起的位置不同
    end
end
taoi=2*Ri/c;
spt=zeros(64,M);
for i=1:64 % 时域信号
    spt(i,:)=10^(SNR_4/20)*exp(1i*2*pi*f0*taoi(i,:)+1i*pi*mu*taoi(i,:).^2);
end
for i=1:64 % 时域噪声加信号
    sp(i,NR:NR+M-1)=sp(i,NR:NR+M-1)+spt(i,:);
end
Ps_PPbef=10*log10(mean(mean(real(spt).^2))); % 脉压前的信号功率
SNR_PPbef=Ps_PPbef-Pn_PPbef; % 脉压前的信噪比
y=zeros(64,NR1);
for i=1:64
    y(i,:)=ifft(fft(sp(i,:),NR1).*conj(fft(st_lfm4,NR1)),NR1); % 脉压输出
end
%最高的峰位在744~747这个位置；
Pn_PPaft=10*log10(mean(mean(real([y(:,1:739),y(:,751:2048)]).^2))); %dB 脉压后的噪声功率
Ps_PPaft=10*log10(mean(mean(real(y(:,744:747)).^2)));
SNR_PPaft=Ps_PPaft-Pn_PPaft; %dB 脉压后的信噪比
%画图
figure(15)% 原始IQ信号
plot(real(sp')),grid on
xlabel('\bf 时域采样点'),ylabel('\bf 幅度')
title('\bf 原始I、Q信号')
xlim([0 2000]);
figure(16) % 脉压输出
plot((0:NR1-1)*c*Ts/2/1000+Rmin/1000,10*log10(abs(y').^2)),grid on
xlabel('\bf 距离 / km'),ylabel('\bf 功率 / dB')
title('\bf 脉压输出结果（Taylor窗）')
ylim([20 45]);

sct=zeros(256,NR1);
for i=1:2048
    sct(:,i)=abs(fftshift(fft(y(:,i),256)));
end
r=((0:NR1-1)*c*Ts/2+Rmin)./1e3;
dp=(-128:127)*(B/128)/1e3;
figure(17) % 距离-多普勒图
mesh(r,dp,sct);
xlabel('\bf 距离 / km'),ylabel('\bf 多普勒频率 / kHz'),zlabel('\bf 幅度')
title('\bf 距离-多普勒图')
figure(18) % 距离-多普勒等高线
contour(r,dp,sct);
xlabel('\bf 距离 / km'),ylabel('\bf 多普勒频率 / kHz')
title('\bf 等高线图')
xlim([79.5 80.2]);ylim([-70 -10]);
% 对最开始的噪声单独进行脉压并相干积累
sp_y=zeros(64,NR1);
for i=1:64
    sp_y(i,:)=ifft(fft(sp_noise(i,:),NR1).*conj(fft(st_lfm4,NR1)),NR1);
end
sp_sct=zeros(256,NR1);
for i=1:2048
    sp_sct(:,i)=abs(fftshift(fft(sp_y(:,i),256)));
end
Pn_CPI=10*log10(mean(mean(real(sp_sct).^2))); %dB 积累后的噪声功率
Ps_CPI=10*log10(mean(mean(real(sct(120:126,744:747)).^2))); %dB 积累后的信号功率
SNR_CPI=Ps_CPI-Pn_CPI; %dB 积累后的信噪比

%-------------------------CFAR---------------------%
P_fa=10.^(-6); % 虚警概率
M_ref=6;%参考单元数
L_slipper=M_ref+1;%滑窗长度
L_move=1;%滑窗间隔
L_num=floor((NR1-L_slipper)/L_move)+1;%滑窗次数
Z=zeros(256,L_num);
sct2=abs(sct).^2;
for k=1:256
    for i=1:L_num
        for j=1:L_slipper
            Z(k,i)=Z(k,i)+sct2(k,(i-1)*L_move+j);
        end
        Z(k,i)=Z(k,i)-sct2(k,(i-1)*L_move+M_ref/2+1);
    end
end
K=P_fa.^(-1/M_ref)-1;
S=mean(Z.*K);
figure(19)% CFAR检测门限
plot((0:NR1-1)*c*Ts/2/1000+Rmin/1000,20.*log10(abs(sct(123,:)')),'k'),hold on
plot((M_ref/2+1:NR1-M_ref/2)*c*Ts/2/1000+Rmin/1000,10.*log10(abs(S))),grid on
xlabel('\bf 距离 / km');ylabel('\bf 功率 / dB')
title('\bf CFAR')
xlim([60 140]);ylim([30 80]);
legend('\bf 目标所在多普勒通道','\bf CFAR门限');
%% =========================5.单脉冲测角=======================%
%-------------------------（1）单脉冲测角---------------------%
sita_3dB=6;%3dB波束宽度,度
delta_sita=sita_3dB/2;%半波束宽度，度
sita_vec=linspace(-sita_3dB,sita_3dB,2000)*pi/180;
F1=exp(-2.778/2*(sita_vec+delta_sita*pi/180).^2/(sita_3dB*pi/180)^2);%方向图1
F2=exp(-2.778/2*(sita_vec-delta_sita*pi/180).^2/(sita_3dB*pi/180)^2);%方向图2
SIGMA=(F1+F2); % 和通道
DELTA=-(F1-F2); % 差通道
normalized_error=real(DELTA.*conj(SIGMA))./(SIGMA.*SIGMA);
figure(20)
plot(sita_vec*180/pi,((F1)),'b-','linewidth',2),hold on
plot(sita_vec*180/pi,((F2)),'r--','linewidth',2),hold on
plot(sita_vec*180/pi,((SIGMA)),'k-','linewidth',2),hold on
plot(sita_vec*180/pi,(abs(DELTA)),'g-.','linewidth',2),hold on
plot(sita_vec*180/pi,((normalized_error)),'y.'),grid on
xlabel('\bf 角度 / \circ'),ylabel('\bf 相对幅度');
title('\bf 和、差波束及其归一化误差信号')
axis tight
legend({'波束1','波束2','和波束','差波束','误差信号'})
%对normalized_error进行曲线拟合
degree=4;
coe=polyfit(normalized_error,sita_vec.*180/pi,degree );
x=linspace(min(normalized_error),max(normalized_error),1000);
y=polyval(coe,x);
figure(21)
plot(sita_vec*180/pi,((normalized_error)),'b--'),hold on
plot(y,x,'r-'),grid on
xlabel('\bf 相对波束中心方位 / \circ'),ylabel('\bf 归一化方位误差信号')
title('\bf 误差曲线与拟合曲线比较')
axis tight
legend({'误差曲线',['拟合曲线（',num2str(degree ),'阶）']})
rang=fix(length(sita_vec)/2-100):fix(length(sita_vec)/2+100);
K_theta=polyfit(sita_vec(rang)*180/pi,normalized_error(rang),1);
fprintf('误差信号斜率为 %5.4f\n',1/K_theta(1));

%-------------------------(2)Monto Carlo测角精度---------------------%
theta=1*pi/180;%目标角度，rad
f1=exp(-2.778/2*(theta+delta_sita*pi/180).^2/(sita_3dB*pi/180)^2);%方向图1
f2=exp(-2.778/2*(theta-delta_sita*pi/180).^2/(sita_3dB*pi/180)^2);%方向图2
SNR_vec=5:1:30;
theta_est=zeros(length(SNR_vec),100);
for monte=1:100
    noise1=(randn(1,1)+1j*randn(1,1))/sqrt(2);
    noise2=(randn(1,1)+1j*randn(1,1))/sqrt(2);
    for SNR_cnt=1:length(SNR_vec)
        SNR=SNR_vec(SNR_cnt);
        x_sum=10^(SNR/20)*(f1+f2)+noise1;
        x_delta=10^(SNR/20)*(-f1+f2)+noise2;
        x_err=real(conj(x_delta).*x_sum)./(x_sum.*conj(x_sum));%real(x_delta./x_sum);
        theta_est(SNR_cnt,monte)=polyval(coe,x_err);
    end      
end
theta_mse=zeros(1,length(SNR_vec));
for SNR_cnt=1:length(SNR_vec)
    theta_mse(SNR_cnt)=sqrt(mean((theta_est(SNR_cnt,:)-theta*180/pi).^2));
end
theroy_val=sita_3dB/2/1.6./sqrt(2*10.^(SNR_vec/10));
figure(22)
plot(SNR_vec,theta_est.'-theta*180/pi,'r.'),grid on
xlabel('\bf SNR / dB'),ylabel('\bf 单次测量误差 / \circ）')
title('\bf SNR与单次测角误差');
axis tight
figure(23)
plot(SNR_vec,theta_mse,'o-','linewidth',2),hold on
plot(SNR_vec,theroy_val,'r--','linewidth',2),grid on
xlabel('\bf SNR / dB'),ylabel('\bf RMSE / \circ')
title('\bf 测角精度')
axis tight
legend('\bf Monto Carlo分析结果','\bf 理论计算值')

%-------------------------(3)和差通道时域脉压---------------------%
TimeWidth=10e-6;%脉冲宽度
PRT=250e-6;%脉冲重复周期，s
PRF=1/PRT;%雷达脉冲重复频率,Hz
BandWidth=10e6;%信号带宽,Hz
mu=BandWidth/TimeWidth;%调频率,Hz/s
fs=2*BandWidth;%采样率，Hz
samplenumber=fix(fs*PRT);%一个脉冲重复周期内的采样点数
number=fix(fs*TimeWidth);%回波在一个脉冲宽度内的采样点数
number=number+mod(number,2);
t_vec=linspace(-TimeWidth/2,TimeWidth/2,number);%快时间域矢量,s
%目标参数
target_Range=20e3;%目标距离，m
tau0=2*target_Range/c;%回波固定的延迟时间
DelayNumber=fix(fs*tau0);% 把目标距离换算成采样点（距离门）
%产生目标信号
Chirp=exp(1j*(pi*mu*t_vec.^2));%产生线性调频信号
coeff=conj(fliplr(Chirp));%产生脉压系数
Signal=zeros(1,samplenumber);% 一个脉冲
Signal(DelayNumber+1:DelayNumber+number)=Chirp;%目标的一个脉冲（未加多普勒速度）
Signal=Signal/max(abs(Signal));
%产生系统噪声信号                        %
SNR=20;%信噪比20dB
SystemNoise1=10^(-SNR/20)/sqrt(2)*(randn(1,samplenumber)+1j*randn(1,samplenumber));
SystemNoise2=10^(-SNR/20)/sqrt(2)*(randn(1,samplenumber)+1j*randn(1,samplenumber));
theta_sumdelta=1.25*pi/180;%方位角,rad
f1=exp(-2.778/2*(theta_sumdelta+delta_sita*pi/180).^2/(sita_3dB*pi/180)^2);%偏角为1.25度时波束2方向图函数
f2=exp(-2.778/2*(theta_sumdelta-delta_sita*pi/180).^2/(sita_3dB*pi/180)^2);%偏角为1.25度时波束2方向图函数
Echo1=Signal+SystemNoise1;%波束1回波
Echo2=Signal+SystemNoise2;%波束1回波
Echo_sum=f1*Echo1+f2*Echo2;%偏角为1.25度时和波束
Echo_dif=-f1*Echo1+f2*Echo2;%偏角为1.25度时差波束
Echo_sum_pc=conv(Echo_sum,coeff);Echo_sum_pc=Echo_sum_pc(number:number+samplenumber-1);%和通道时域脉压
Echo_dif_pc=conv(Echo_dif,coeff);Echo_dif_pc=Echo_dif_pc(number:number+samplenumber-1);%差通道时域脉压
sum=Echo_sum_pc(abs(Echo_sum_pc)==max(abs(Echo_sum_pc)));
delta=Echo_dif_pc(abs(Echo_dif_pc)==max(abs(Echo_dif_pc)));
error=real(sum.*conj(delta))./(sum.*conj(sum));
Echo_err=real(Echo_sum_pc.*conj(Echo_dif_pc))./(Echo_sum_pc.*conj(Echo_sum_pc));%归一化误差
figure(24)
plot((0:1/fs:PRT-1/fs)*c/2*1e-3,20*log10(abs(Echo_sum_pc))),hold on
plot((0:1/fs:PRT-1/fs)*c/2*1e-3,20*log10(abs(Echo_dif_pc)),'r--'),grid on
xlabel('\bf 距离 / km'),ylabel('\bf 幅度 / dB')
title('\bf 弹目距离20Km时和、差通道时域脉压结果')
legend('\bf 和通道','\bf 差通道')
fprintf('设定20km处方位角%5.4f°，误差信号为 %5.4f ，测量值为 %5.4f°\n',theta_sumdelta*180/pi,error,polyval(coe,error));
% figure;plot((0:1/Fs:PRT-1/Fs)*C/2*1e-3,20*log10(Echo_err));title('弹目距离20Km时归一化误差');
% xlabel('距离/Km');ylabel('归一化误差');grid on;
