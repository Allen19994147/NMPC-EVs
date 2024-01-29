%%% Nonlinear Model Predictive Control (NMPC) of Electric Vehicle
%%% Allen Lee
clc
clear all

%%% Global Variables
global AccessaryLoad BrakeRegenFraction Crr DrivelineModelG DrivelineModelSpinLoss EnergyCapacity MassVehInertia MotorModelC
global MotorModelKc MotorModelKi MotorModelKw Rint Voc massVeh wheelRadius MaxTorqueOutput MaxPowerOutput gravity MaximumBreakForce
global N Ts
global n m vr v Ptr FBr SOC SOCr Vtrajectory
%__________________________________________________________________________
AccessaryLoad = 600;BrakeRegenFraction=0.55;Crr=0.01;DrivelineModelG=3.55;
DrivelineModelSpinLoss=6;EnergyCapacity=12.6;MassVehInertia=2392;MaxPowerOutput=1e5;
MaxTorqueOutput=500;MaximumBreakForce=10000;MotorModelC=628.2974;
MotorModelKc=0.045237;MotorModelKi=0.0167;MotorModelKw=5.06640055940796E-5;
Q=3000;Rint=0.1;Ts=0.01;Voc=340;gravity=9.81;massVeh=2300;wheelRadius=0.34;
%__________________________________________________________________________
N = 5; % Control/Prediction Horizen
n = 2; % #state
m = 2; % # input

%%% Equilibrium References

data = load("US06.mat");


data1 = data.VelocityReference.Reference_Speed;
data2 = data.VelocityReference.signal2;
% figure
% hold on
% plot(data1.Time,data1.Data)
% plot(data2.Time,data2.Data)
% hold off

Run_Sample = size(data1.Data,1);

TotalVr = data1.Data./2.237; %mph2m/s
Vr_preview = [TotalVr;zeros(N,1)];
Vr_real = [zeros(N,1);TotalVr];

Vahead = 0;
%%% References and IC
x0 = [0,100]';
u0 = [0,0]';
Xall = zeros(n,N+Run_Sample);
Uall = zeros(m,N+Run_Sample);
Xall(:,1) = x0;
x = x0;
Vtrajectory = zeros(N+1,1);

Xk = zeros(n*(N+1),1); % all future x from step k
Uk = zeros(m*N,1); % all future u from step k
Xk(1:n,1) = x0;
Uk(1:m,1) = u0;

zek = [Xk;Uk];

G_rw = DrivelineModelG/wheelRadius;

% Build QX,RU,H
xmin = [-1e5;5]; % 0 to fullfill inequality
xmax = [1e5;100]; % 0 to fulfill inequality
umin = [0;-1e5];

Fx = [0 0;0 1;0 0;0 -1]; % No constraints on velcity
Fu = [eye(m);-eye(m)];

Q_ = [Q 0;0 1];
QN = Q_;
R = eye(m);
Qx = Q_;
RU = R;
FX = Fx;
FU = Fu;
for i = 1:N-1
    Qx = blkdiag(Qx,Q_);
    RU = blkdiag(RU,R);
    FX = blkdiag(FX,Fx);
    FU = blkdiag(FU,Fu);
end
FX = blkdiag(FX,Fx);
Qx = blkdiag(Qx,QN);
H = blkdiag(Qx,RU.*0);

% Simulation with NMPC
ObjFunc = @(ze) ze'*H*ze;
F = blkdiag(FX,FU);
Feq = []; Geq = []; lb = []; ub = []; % Linear constraints

nonlcon = @NonlinearConstraints;
options = optimoptions('fmincon','Display','none','Algorithm','sqp');

SOCr=100; vr=0;v=0;Ptr=zeros(N,1);FBr=zeros(N,1);SOC=100;Pt=0;FB=0;
% for i = 1:Run_Sample+N
for i = 1:3000

    Vahead = Vr_preview(i,1); 
    vr = Vr_real(i,1);
    SOCr = SOC;
    Vtrajectory(1:end-1,1) = Vtrajectory(2:end,1);
    Vtrajectory(end,1) = Vahead;
    x = Xall(:,i); 
    MotorSpeed = max(G_rw*abs(v),1e-5); % MotorSpeed must be positive
    umax = [min(MaxTorqueOutput,MaxPowerOutput/MotorSpeed);0];

    Gxe = [];Gue=[];G=[];Ptr=zeros(size(Ptr));FBr=zeros(size(FBr));
    for j = 1:N
        xr = [Vtrajectory(j,1);SOCr]; %future references
        [Pt,FB] = desiredControl(xr(1),Vtrajectory(j+1),Pt,FB); % Find equilibrium control using fd
        Ptr(j,1)=Pt;FBr(j,1)=FB;

        gxe = [xmax;-xmin] - Fx*xr;
        gue = [umax;-umin] - Fu*[Pt;FB];
        Gxe = [Gxe;gxe];
        Gue = [Gue;gue];
    end
    xr = [Vtrajectory(end,1);SOCr]; %future references
    [Ptr(end,1),FBr(end,1)] = desiredControl(xr(1),xr(1),Pt,FB);
    gxe = [xmax;-xmin] - Fx*xr;
    Gxe = [Gxe;gxe];
    G = [Gxe;Gue];
    zek = fmincon(ObjFunc,zek,F,G,Feq,Geq,lb,ub,nonlcon,options);
    Uall(:,i) = zek(n*(N+1)+1:n*(N+1)+m,1) + [Ptr(1,1);FBr(1,1)];
    if(vr>v)
        Pt=Uall(1,i);FB=0;Uall(2,i)=0;
    else
         Pt=0;FB=Uall(2,i);Uall(1,i)=0;
    end
    [v_k1,soc_k1] = fd2(v,SOC,Pt,FB);
    Xall(1,i+1) = v_k1;
    Xall(2,i+1) = soc_k1;
    SOC = soc_k1;
    v = v_k1;
end

%%
figure
hold on
plot(Vr_real,'r-')
plot(Xall(1,:),'b-')
% plot(Xall(2,:),'r.-')
legend('vel','true')
hold off

figure
hold on
for i = 1:m
    plot(Uall(i,:),'.-')
end
hold off
%% Save Data
filename = 'The Best one';
save(filename,"Vr_real","Xall","Uall")

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ptr, FBr] = desiredControl(vr,vr1,Pt,FB)
global massVeh gravity Crr DrivelineModelSpinLoss DrivelineModelG wheelRadius MassVehInertia
global BrakeRegenFraction AccessaryLoad  MotorModelKc MotorModelKw MotorModelKi MotorModelC Ts

G_rw = DrivelineModelG/wheelRadius;
C1 = 0.5*vr^2 + massVeh*gravity*Crr + DrivelineModelSpinLoss*G_rw;
C2 = AccessaryLoad + MotorModelKw*(G_rw*vr)^3 + MotorModelKi*G_rw*vr + MotorModelC;
root = roots([MotorModelKc,G_rw*vr,C2]);
if(isreal(root))
    TM_r = max(root);

    u2r = C1 - TM_r*G_rw + (vr1-vr)*MassVehInertia/Ts;
    u1r = TM_r + BrakeRegenFraction/G_rw*u2r;
    Ptr = u1r;
    FBr = u2r;
else
    Ptr = Pt;
    FBr = FB;
end

end

function [c,ceq] = NonlinearConstraints(z)
global v vr SOC SOCr Ptr FBr Vtrajectory
global N n m
c=[];
shift = n*(N+1);
c1 = zeros(n*N,1);
% v_r=0;Pt_r=0;FB_r=0;
for i=1:N
    vz = z(n*(i-1)+1,1);
    socz = z(i*n,1);
    Ptz = z(shift+m*(i-1)+1,1);
    FBz = z(shift+i*m,1);
    Pt_r = Ptr(i);FB_r=FBr(i); v_r = Vtrajectory(i,1);
    [v1,soc1] = fd2(vz+v_r,socz+SOCr,Ptz+Pt_r,FBz+FB_r);
    c1(n*(i-1)+1) = z(i*n+1,1)+v_r - v1;
    c1(n*i) = z((i+1)*n,1)+SOCr - soc1;
end
ceq = [z(1)+Vtrajectory(1,1)-v;z(2)+SOCr-SOC;c1];
end


function [v_prime,soc_prime] = fd2(v,SOC,PositiveTorque,FB)
global Ts
[dv1,dsoc1] = fc(v,PositiveTorque,FB);
v_prime = v + Ts*dv1;
soc_prime = SOC + Ts*dsoc1;
end

function [dv,dsoc] = fc(v,PositiveTorque,FB)
global MassVehInertia massVeh gravity Crr DrivelineModelSpinLoss DrivelineModelG
global wheelRadius BrakeRegenFraction Voc EnergyCapacity MaxTorqueOutput MaxPowerOutput

G_wr = DrivelineModelG/wheelRadius;
MotorSpeed = max(G_wr*abs(v),1e-5);
if(abs(v)<5/2.237)
    RegenTorque = 0;
else
    RegenTorque = BrakeRegenFraction*FB/G_wr; % Negative
    RegenTorque = max(RegenTorque,-0.5*min(MaxTorqueOutput,MaxPowerOutput/MotorSpeed));
end
BrakeForce = max((1-BrakeRegenFraction)*FB,-1e4);
MotorTorque = PositiveTorque + RegenTorque;
if(abs(v)~=0)
    TractiveForce = (MotorTorque - DrivelineModelSpinLoss)*G_wr  + BrakeForce;
else
    TractiveForce = MotorTorque*G_wr  + BrakeForce;
end
dv = (TractiveForce - 0.5*v^2 - massVeh*gravity*Crr)/MassVehInertia;
P = MotorPower(MotorTorque,MotorSpeed);
I = Current(P);
dsoc = -Voc*I/(EnergyCapacity*1000*3600)*100;% in %
end

function P = MotorPower(MotorTorque,MotorSpeed)
global AccessaryLoad  MotorModelKc MotorModelKw MotorModelKi MotorModelC v
if(MotorSpeed>1e-3)
    Powerloss = MotorModelKc*MotorTorque^2 +... 
         MotorModelKw*(MotorSpeed)^3 + MotorModelKi*MotorSpeed + MotorModelC;
else
    Powerloss =  0;
end
P = Powerloss + MotorSpeed*MotorTorque + AccessaryLoad;
end

function I = Current(P)
global Voc Rint
dummy = Voc^2 - 4*Rint*P;
if(dummy<0)
    error("Power is too large")
end
I = (Voc - sqrt(dummy))/(2*Rint);
end


