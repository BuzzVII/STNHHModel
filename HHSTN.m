
function [T I V]=HHSTN(tmax,Imaxin,Istartin,Istopin)
% Runs a HH model of an STN neuron for tmax ms, where a rectangle current can be
% applied with Imaxin as the current value and Istartin and Istopin as the
% rising and falling edge times.


%% VARIABLES
% Followin "Activity Patterns in a Model for the Subthalamopallidal Network
% of the Basal Ganglia" - Rubin et al. The J. of Neurosci. 22(7) 2002

global C gNabar gKbar gLbar gCabar gTbar gAHP ENa EK EL ECa tauh1 taun1 taur1 tauh0 taun0 taur0 phih phin phir k1 kCa epsilon thetam thetah thetan thetar thetaa thetab thetas thetaht thetant thetart sigm sigh sign sigr siga sigb sigs sight signt sigrt Istart Istop Imax
%make global so ODE15s can reach them

C = 1.0; 
gNabar = 37.5;
gKbar = 45;
gLbar = 2.25;
gCabar = 0.5;
gTbar =  0.5;
gAHP = 9;
ENa = 55;
EK = -80;
EL = -60;
ECa = 140;
tauh1 = 500;
taun1 = 100;
taur1 = 17.5;
tauh0 = 1;
taun0 = 1;
taur0 = 40;
phih  = 0.75;
phin =0.75;
phir = 0.2;
k1= 15;
kCa = 22.5;
epsilon = 3.75e-5;
thetam = -30;
thetah = -39; 
thetan = -32;
thetar =-67;
thetaa = -63;
thetab= 0.4;
thetas=-39;
thetaht=-57;
thetant=-80;
thetart=68;
sigm=15;
sigh=-3.1;
sign=8;
sigr=-2;
siga=7.8;
sigb=-0.1;
sigs=8;
sight=-3;
signt=-26;
sigrt=-2.2;
vrest = -60;

%% Initial Values 
% set v, h, n, r and cCa to rest state

v = vrest;
h = xinf(v,thetah,sigh)/10;
n = xinf(v,thetan,sign)*2;
r = xinf(v,thetar,sigr);
cCa = 0.04825;


%% Solve
% using ode45

% Initial Conditions 
s0 = [h n r cCa v];

% Time settings
tspan = [0,tmax];
Imax=Imaxin;
Istart=Istartin;
Istop=Istopin;

tic;
options=odeset('InitialStep',10^(-3),'MaxStep',10^(-1));
[T,S] = ode15s(@fxn,tspan,s0,options);
toc

figure(1)
clf
plot(T,S(:,5)), xlabel('time'), ylabel('Voltage')
hold on
plot([Istart,Istop],[Imax,Imax],'r','linewidth',4)

figure(2)
clf
plot(T,S(:,1:4)), legend('h','n','r','[Ca]')
hold on
plot([Istart,Istop],[Imax/100,Imax/100],'r','linewidth',4)

figure(3)
clf
subplot(2,2,1)
plot(S(:,5),S(:,1)),xlabel('v'), ylabel('h')
subplot(2,2,2)
plot(S(:,5),S(:,2)),xlabel('v'), ylabel('n')
subplot(2,2,3)
plot(S(:,5),S(:,3)),xlabel('v'), ylabel('r')
subplot(2,2,4)
plot(S(:,5),S(:,4)),xlabel('v'), ylabel('[Ca]')

figure(4)
I2=-diff(S(:,5))./diff(T)*C;
I=(gLbar*(S(:,5)-EL) +gK*(S(:,5)-EK) +gNa*(S(:,5)-ENa) +   gT*(S(:,5)-ECa)+ gCa*(S(:,5)-ECa)+gAHP*(S(:,5)-EK).*(S(:,4)./(S(:,4)+k1)));
V=S(:,5);
tI=T(1:length(T)-1);
plot(tI,I2),xlabel('Time'),ylabel('Current');
hold on
plot(T,I,'red');



    function xinfout = xinf(vin, thetain, sigmain)
        xinfout = 1/(1+exp(-(vin-thetain)/sigmain));
    end
    
    function binfout = binf(rin)
        binfout = 1/(1+exp((rin-thetab)/sigb))-1/(1+exp(-thetab/sigb));
    end
    
    function tauout = tau(vin, tau0, tau1, thetat, sigt)
        tauout= tau0+tau1/(1+exp(-(vin-thetat)/sigt));
    end
        
function ds = fxn(t,s)
    ds = zeros(5,1);
    ds(1) = phih*((xinf(s(5),thetah,sigh)-s(1))/tau(s(5),tauh0,tauh1,thetaht,sight));
    ds(2) = phin*((xinf(s(5),thetan,sign)-s(2))/tau(s(5),taun0,taun1,thetant,signt));
    ds(3) = phir*((xinf(s(5),thetar,sigr)-s(3))/tau(s(5),taur0,taur1,thetart,sigrt));
    gNa = gNabar*(xinf(s(5),thetam,sigm)^3)*s(1);
    gK = gKbar*(s(2)^4);
    gT = gTbar*xinf(s(5),thetaa,siga)^3*binf(s(3))^2;
    gCa= gCabar*xinf(s(5),thetas,sigs)^2;
    ds(4)=  epsilon*(-gCa*(s(5)-ECa)-gT*(s(5)-ECa)-kCa*s(4));
    Iin = Imax*(heaviside(t-Istart)-heaviside(t-Istop));
    ds(5) = -(1/C)*(gLbar*(s(5)-EL) +gK*(s(5)-EK) +gNa*(s(5)-ENa) +   gT*(s(5)-ECa)+ gCa*(s(5)-ECa)+gAHP*(s(5)-EK)*(s(4)/(s(4)+k1))) + Iin/C;
    %keyboard;
end
    
    
end