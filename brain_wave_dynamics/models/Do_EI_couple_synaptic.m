%I-population Ne_1=1000 and E-population Ne_2=4000
%Each neuron is a nicely modified theta model

%probability of synaptic connections is 0.2
%synaptic coupling is averaged. So each neuron in population
%receives the same stimulus .

%XX_G are synaptic dynamics, v_1(phase),v_2(phase) are membrane potentials
%note that phase 0 means 0 rad. Phase "nx/2" means Pi rad
%(This is different from other programs)

%dt=0.01ms num is total number of time steps.(therefore, here is a 0.3s
%simulation in total)
%Raster plot and firing rate are plotted
%Red in figure means excitatory neurons

Ne_1=1000; %I
Ne_2=4000; %E
meanrate = [];
inout=[];
dt = 0.01; %0.005 0.0025
a_master=2;
a=a_master;
sig = 2;%sqrt(2*D);
num=30000;%60000

%for w=1 %=mu
theta=pi;
V_R=-62;
V_T=-55;
I_V_syn=-70;
E_V_syn=0;
c_1=2/(V_T-V_R);

I_g_L=0.1;
E_g_L=0.05;

II_c_2=(2*I_V_syn-V_T-V_R)/(V_T-V_R);
IE_c_2=(2*E_V_syn-V_T-V_R)/(V_T-V_R);
EI_c_2=(2*I_V_syn-V_T-V_R)/(V_T-V_R);
EE_c_2=(2*E_V_syn-V_T-V_R)/(V_T-V_R);


t_1 = 100; %delay step to transmit
II_t_r = 0.5; %refractory step 20
II_t_d= 5;%tau for inhibit
II_gbar = 0.138062206;
II_rate = 200;
IE_t_r = 0.5;
IE_t_d= 2;
IE_gbar = 0.010400214; % 0; 0.010400214;
IE_rate = 800;
EI_t_r = 0.5;
EI_t_d= 5;%tau for inhibit
EI_gbar = 0.172577757; %0.172577757; 0;
EI_rate = 200;
EE_t_r = 0.5;
EE_t_d= 2;
EE_gbar = 0.01291816; %0.01291816; 0;
EE_rate = 800;


II_c_3=-1/II_t_r/II_t_d;
II_c_4=-(II_t_r+II_t_d)/II_t_r/II_t_d;
II_c_5=II_gbar*II_rate/II_t_r/II_t_d;
IE_c_3=-1/IE_t_r/IE_t_d;
IE_c_4=-(IE_t_r+IE_t_d)/IE_t_r/IE_t_d;
IE_c_5=IE_gbar*IE_rate/IE_t_r/IE_t_d;
EI_c_3=-1/EI_t_r/EI_t_d;
EI_c_4=-(EI_t_r+EI_t_d)/EI_t_r/EI_t_d;
EI_c_5=EI_gbar*EI_rate/EI_t_r/EI_t_d;
EE_c_3=-1/EE_t_r/EE_t_d;
EE_c_4=-(EE_t_r+EE_t_d)/EE_t_r/EE_t_d;
EE_c_5=EE_gbar*EE_rate/EE_t_r/EE_t_d;

 EE_G= zeros(2,num);
 EI_G= zeros(2,num);
 IE_G= zeros(2,num);
 II_G= zeros(2,num);
%tmp_I2_1= 0;
%diff_I2_1=0;
%tmp_diff_I2_1=0;

%tmp_I2_2= 0;
%diff_I2_2=0;
%tmp_diff_I2_2=0;

K = -0.0; %-0.3

% R = 0.1 *power(10,1);
%S=[0.5*rand(Ne,Ne)];
for i=1:Ne_1
for j=1:Ne_1
    tmp = rand(1);
    if tmp > -5; %0.9;
        S(i,j)=1;
    else S(i,j)=0;
    end
end
end

%S=[ones(Ne,Ne)];

v_1=2*pi*rand(Ne_1,1)-pi;    % Initial values of v
v_2=2*pi*rand(Ne_2,1)-pi;    % Initial values of v

firings_1=[];             % spike timings
firings_2=[];             % spike timings
A_1 = zeros(num,1);
A_2 = zeros(num,1);
out1 = zeros(num,1);
out2 = zeros(num,1);
Nave=200;
in = zeros(num,1);
c=0.0;

I2_1=zeros(num,1);
I2_2=zeros(num,1);


for k=1:num
%  I_common= sigma * randn(Ne_1,1); % only for noisesyn
  I_common= randn(1,1); % only for noisesyn
  I_1= sig*randn(Ne_1,1); % input
  In(k) = 0.0;
  %In(t) = 0.5 * sin(2 * pi * f * t / 2000);
%  I_1= sigma * (sqrt(1-c) * I_1 + sqrt(c) * I_common);%1/sqrt(2) * (I_1 + I_common);
%   I_1= I_1 + In(i);
  if k>0
  fired_1=find(v_1>=theta);    % indices of spikes
  A_1(k)=size(fired_1,1)/dt / Ne_1;
  firings_1=[firings_1; k+0*fired_1,fired_1]; %i-1= actual firing time
 v_1(fired_1)=v_1(fired_1)-2*pi;
% ref_1 = find(firings_1(:,1)>(i-1-tauR));
% v_1(firings_1(ref_1,2))=0;
  end

  I_2= sig*randn(Ne_2,1); % input
  In(k) = 0.0;
  %In(t) = 0.5 * sin(2 * pi * f * t / 2000);
%  I_2=  sigma * (sqrt(1-c) * I_2 + sqrt(c) * I_common);%1/sqrt(2) * (I_2 + I_common);
%  I_2= I_2 + In(i);
if k>1
fired_2=find(v_2>=theta);    % indices of spikes
  A_2(k)=size(fired_2,1)/dt / Ne_2;
  firings_2=[firings_2; k+0*fired_2,fired_2]; %i-1= actual firing time
 v_2(fired_2)=v_2(fired_2)-2*pi;
% ref_2 = find(firings_2(:,1)>(i-1-tauR));
% v_2(firings_2(ref_2,2))=0;
end

 %  feedback = find(firings(:,1)==(i-5));
%  I2 =  sum(S(:,firings(feedback,2)),2);
  if k > t_1
    II_G(1,k+1)= II_G(1,k)+ dt*II_G(2,k);
    II_G(2,k+1)=II_c_3*dt*II_G(1,k)+(1+II_c_4*dt)*II_G(2,k)+II_c_5*A_1(k-t_1)*dt;%flux no need
    IE_G(1,k+1)= IE_G(1,k)+ dt*IE_G(2,k);
    IE_G(2,k+1)=II_c_3*dt*IE_G(1,k)+(1+IE_c_4*dt)*IE_G(2,k)+IE_c_5*A_2(k-t_1)*dt;
    EI_G(1,k+1)= EI_G(1,k)+ dt*EI_G(2,k);
    EI_G(2,k+1)=II_c_3*dt*EI_G(1,k)+(1+EI_c_4*dt)*EI_G(2,k)+EI_c_5*A_1(k-t_1)*dt;
    EE_G(1,k+1)= EE_G(1,k)+ dt*EE_G(2,k);
    EE_G(2,k+1)=II_c_3*dt*EE_G(1,k)+(1+EE_c_4*dt)*EE_G(2,k)+EE_c_5*A_2(k-t_1)*dt;
  end

  v_1=v_1+(-I_g_L * cos(v_1) + c_1 * (1+cos(v_1)).*(a - sig*sig/2*sin(v_1)) +II_G(1,k)*(II_c_2 * (1+cos(v_1))-sin(v_1)) +IE_G(1,k)*(IE_c_2 * (1+cos(v_1))-sin(v_1)))*dt + c_1 * (1+cos(v_1)).* I_1 * sqrt(dt);
  % g_L minus? sin-> minus


  v_2=v_2+(-E_g_L * cos(v_2) + c_1 * (1+cos(v_2)).*(a - sig*sig/2*sin(v_2)) +EI_G(1,k)*(EI_c_2 * (1+cos(v_2))-sin(v_2)) +EE_G(1,k)*(EE_c_2 * (1+cos(v_2))-sin(v_2)))*dt + c_1 * (1+cos(v_2)).* I_2 * sqrt(dt);

end
figure;
plot(firings_1(:,1),firings_1(:,2),'.');
% axis([27000 30000 1 500]);
hold on;
plot(firings_2(:,1),firings_2(:,2)+1000,'r.');
 axis([0 20000 1 5000]);



for i=1:num
    if i>Nave
    out1(i-(Nave)/2) = mean(A_1(i-Nave:i));
    out2(i-(Nave)/2) = mean(A_2(i-Nave:i));
    end
end


figure;
plot(A_1);
hold on;
plot(A_2,'r');
axis([0 num 0 1]);


%plot(In);
%meanrate(w+1)=size(firings,1)/Ne/(T*dt);
meanrate = [meanrate; mean(out1(1000:10000))];

%tmp = corrcoef(out1(1000:10000),In(1000:10000))
%inout = [inout; tmp(1,2)];

%end
