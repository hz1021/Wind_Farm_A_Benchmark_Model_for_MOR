clc;
clear all;

parSimulacio.ti=0;
parSimulacio.tf=10;
parSimulacio.tol=1e-9;
parSimulacio.pas=1/(60*3e3);
parSimulacio.captura.Tmostreig=0.5e-4;
parSimulacio.captura.Delmat=1;

parFont.f=50;
parFont.V_D=66e3;
parFont.V_I=200;
parFont.phi_D=0;
parFont.phi_I=22*pi/180;
parFont.tt_D=[0,1e5];
parFont.tA_D=[1,1];
parFont.tt_I=[0,1e5];
parFont.tA_I=[0,0];

parMecanica.R=40; %radi
parMecanica.taupitch=0.1;
parMecanica.rho=1.225; %densistat de l'aire
parMecanica.nu=90; %nu
parMecanica.It=4e6;%inercia
parMecanica.parCp.c1=1; 
parMecanica.parCp.c2=39.52;
parMecanica.parCp.c6=2.04;
parMecanica.parCp.c7=14.47;
parMecanica.parCp.cpopt=parMecanica.parCp.c1*parMecanica.parCp.c2*exp(-(parMecanica.parCp.c2+parMecanica.parCp.c6*parMecanica.parCp.c7)/parMecanica.parCp.c2)/parMecanica.parCp.c7;
parMecanica.parCp.lambdaopt=parMecanica.parCp.c2*parMecanica.parCp.c7/(parMecanica.parCp.c2+parMecanica.parCp.c6*parMecanica.parCp.c7);
parMecanica.vwn=9; %velocitat vent nominal
parMecanica.vw0=7;
parMecanica.tt_vw=[0,1e5];
parMecanica.tA_vw=[7,7];
parMecanica.omegatn=parMecanica.parCp.lambdaopt*parMecanica.vwn/parMecanica.R;
parMecanica.omegat0=parMecanica.parCp.lambdaopt*parMecanica.vw0/parMecanica.R;
parMecanica.gammatopt=1/2*parMecanica.rho*pi*parMecanica.R^2*parMecanica.vwn^3*parMecanica.parCp.cpopt/parMecanica.omegatn; %possible kcp!!
parMecanica.KCP=parMecanica.gammatopt/(parMecanica.omegatn^2);
parMecanica.Prated=1/2*parMecanica.rho*pi*parMecanica.R^2*parMecanica.vwn^3*parMecanica.parCp.cpopt;

parGenerador.rs=0.015;
parGenerador.Lq=0.04/(2*pi*50);
parGenerador.Ld=0.0401/(2*pi*50);
parGenerador.lambda_m=2.35;
parGenerador.P=2;
parGenerador.isqd0=[0;0];

parControl.Electric.Generador.rs=parGenerador.rs;
parControl.Electric.Generador.Lq=parGenerador.Lq;
parControl.Electric.Generador.Ld=parGenerador.Ld;
parControl.Electric.Generador.lambda_m=parGenerador.lambda_m;
parControl.Electric.Generador.P=parGenerador.P;
parControl.Electric.Generador.qalpha=5e2;
parControl.Electric.Generador.dalpha=parControl.Electric.Generador.qalpha;
parControl.Electric.Generador.T=1/3e3;
parControl.Electric.Generador.qKp=parControl.Electric.Generador.qalpha*parControl.Electric.Generador.rs*parControl.Electric.Generador.T/2*((1+exp(-parControl.Electric.Generador.rs/parControl.Electric.Generador.Lq*parControl.Electric.Generador.T))/(1-exp(-parControl.Electric.Generador.rs/parControl.Electric.Generador.Lq*parControl.Electric.Generador.T)));
parControl.Electric.Generador.qKi=parControl.Electric.Generador.qalpha*parControl.Electric.Generador.rs;
parControl.Electric.Generador.dKp=parControl.Electric.Generador.dalpha*parControl.Electric.Generador.rs*parControl.Electric.Generador.T/2*((1+exp(-parControl.Electric.Generador.rs/parControl.Electric.Generador.Ld*parControl.Electric.Generador.T))/(1-exp(-parControl.Electric.Generador.rs/parControl.Electric.Generador.Ld*parControl.Electric.Generador.T)));
parControl.Electric.Generador.dKi=parControl.Electric.Generador.dalpha*parControl.Electric.Generador.rs;
parControl.Electric.Generador.qGcnum=[(parControl.Electric.Generador.qKi*parControl.Electric.Generador.T+2*parControl.Electric.Generador.qKp)/2, (parControl.Electric.Generador.qKi*parControl.Electric.Generador.T-2*parControl.Electric.Generador.qKp)/2];
parControl.Electric.Generador.qGcden=[1 -1];
parControl.Electric.Generador.dGcnum=[(parControl.Electric.Generador.dKi*parControl.Electric.Generador.T+2*parControl.Electric.Generador.dKp)/2, (parControl.Electric.Generador.dKi*parControl.Electric.Generador.T-2*parControl.Electric.Generador.dKp)/2];
parControl.Electric.Generador.dGcden=[1 -1];
parControl.Electric.Generador.SVPWM.T=parControl.Electric.Generador.T;
parControl.Electric.Generador.SVPWM.resolucio=64;

parConvertidor.C=1e-2;
parConvertidor.E0=2600;
parConvertidor.Ll=1e-3;
parConvertidor.rl=0.02;
parConvertidor.Urated=970;
parConvertidor.Pbase=parMecanica.Prated;
parConvertidor.Zbase=(parConvertidor.Urated/sqrt(3))^2/(parConvertidor.Pbase/3);
parConvertidor.frated=50;
parConvertidor.Xlpu=0.2;
parConvertidor.Xl=parConvertidor.Xlpu*parConvertidor.Zbase;
parConvertidor.Ll=parConvertidor.Xl/(2*pi*parConvertidor.frated);
parConvertidor.Xcfilt=(parConvertidor.Urated/sqrt(3))^2/(0.1*parConvertidor.Pbase/3);
parConvertidor.Cfilt=1/(2*pi*parConvertidor.frated*parConvertidor.Xcfilt);

% parameters of the collection grid
parGrid.WTtransf.Pbase=parMecanica.Prated;
parGrid.WTtransf.U2rated=970;
parGrid.WTtransf.frated=50;
parGrid.WTtransf.Zbase2=(parGrid.WTtransf.U2rated/sqrt(3))^2/(parGrid.WTtransf.Pbase/3);
parGrid.WTtransf.rt=0.01*parGrid.WTtransf.Zbase2;
parGrid.WTtransf.Xtpu=0.05;
parGrid.WTtransf.Xt=parGrid.WTtransf.Xtpu*parGrid.WTtransf.Zbase2;
parGrid.WTtransf.Lt=parGrid.WTtransf.Xt/(2*pi*parGrid.WTtransf.frated);

parGrid.WTtransf.rt2=0;
parGrid.WTtransf.Lt2=0;


parGrid.WTstring.Pbase=parMecanica.Prated*10;
parGrid.WTstring.Urated=66e3;
parGrid.WTstring.frated=50;
parGrid.WTstring.Zbase=(parGrid.WTstring.Urated/sqrt(3))^2/(parGrid.WTstring.Pbase/3);
parGrid.WTstring.rcable=0.01*parGrid.WTstring.Zbase;
parGrid.WTstring.Xcablepu=0.02/100;
parGrid.WTstring.Xcable=parGrid.WTstring.Xcablepu*parGrid.WTstring.Zbase;
parGrid.WTstring.Lcable=parGrid.WTstring. Xcable/(2*pi*parGrid.WTstring.frated);

parGrid.WFtransf.f=50;
parGrid.WFtransf.V_D=66e3;
parGrid.WFtransf.V_I=200;
parGrid.WFtransf.phi_D=0;
parGrid.WFtransf.phi_I=22*pi/180;
parGrid.WFtransf.tt_D=[0,1e5];
parGrid.WFtransf.tA_D=[1,1];
parGrid.WFtransf.tt_I=[0,1e5];
parGrid.WFtransf.tA_I=[0,0];
parGrid.WFtransf.Pbase=parGrid.WTstring.Pbase*10;
parGrid.WFtransf.Xpu=0.15;
parGrid.WFtransf.rpu=0.01;
parGrid.WFtransf.Zbase=(parGrid.WFtransf.V_D/sqrt(3))^2/(parGrid.WFtransf.Pbase/3);
parGrid.WFtransf.X=parGrid.WFtransf.Xpu*parGrid.WFtransf.Zbase;
parGrid.WFtransf.L=parGrid.WFtransf.X/(2*pi*parConvertidor.frated);
parGrid.WFtransf.r=parGrid.WFtransf.rpu*parGrid.WFtransf.Zbase;
parGrid.WFtransf.Xcfilt=(parGrid.WFtransf.V_D/sqrt(3))^2/(0.05*parGrid.WFtransf.Pbase/3);
parGrid.WFtransf.rloss=(parGrid.WFtransf.V_D/sqrt(3))^2/(0.005*parGrid.WFtransf.Pbase/3);
parGrid.WFtransf.Cfilt=1/(2*pi*parConvertidor.frated*parGrid.WFtransf.Xcfilt);

parGrid.WTstring.Xc_cable=(parGrid.WTstring.Urated/sqrt(3))^2/(0.1/5*parGrid.WTstring.Pbase/3);
parGrid.WTstring.Rparallel_cable=(parGrid.WTstring.Urated/sqrt(3))^2/(0.02*parGrid.WTstring.Pbase/3);
parGrid.WTstring.C_cable=1/(2*pi*parGrid.WTstring.frated*parGrid.WTstring.Xc_cable); 
parGrid.N1N2=parGrid.WTstring.Urated/parGrid.WTtransf.U2rated; 

parControl.Electric.BusContinua.C=parConvertidor.C;
parControl.Electric.BusContinua.zeta=0.8;
parControl.Electric.BusContinua.omegan=2*pi*6;
parControl.Electric.BusContinua.T=1/3e3;
parControl.Electric.BusContinua.Kp=2*parControl.Electric.BusContinua.zeta*parControl.Electric.BusContinua.omegan*parControl.Electric.BusContinua.C;
parControl.Electric.BusContinua.Ki=parControl.Electric.BusContinua.omegan^2*parControl.Electric.BusContinua.C;
parControl.Electric.BusContinua.Gcnum=[(parControl.Electric.BusContinua.Ki*parControl.Electric.BusContinua.T+2*parControl.Electric.BusContinua.Kp)/2, (parControl.Electric.BusContinua.Ki*parControl.Electric.BusContinua.T-2*parControl.Electric.BusContinua.Kp)/2];
parControl.Electric.BusContinua.Gcden=[1 -1];
parControl.Electric.BusContinua.Easterisc=2600;

parControl.Electric.BandaXarxa.Ll=parConvertidor.Ll;
parControl.Electric.BandaXarxa.rl=parConvertidor.rl;
parControl.Electric.BandaXarxa.alpha=5e2;
parControl.Electric.BandaXarxa.T=1/3e3;
parControl.Electric.BandaXarxa.Kp=parControl.Electric.BandaXarxa.alpha*parControl.Electric.BandaXarxa.rl*parControl.Electric.BandaXarxa.T/2*((1+exp(-parControl.Electric.BandaXarxa.rl/parControl.Electric.BandaXarxa.Ll*parControl.Electric.BandaXarxa.T))/(1-exp(-parControl.Electric.BandaXarxa.rl/parControl.Electric.BandaXarxa.Ll*parControl.Electric.BandaXarxa.T)));
parControl.Electric.BandaXarxa.Ki=parControl.Electric.BandaXarxa.alpha*parControl.Electric.BandaXarxa.rl;
parControl.Electric.BandaXarxa.Gcnum=[(parControl.Electric.BandaXarxa.Ki*parControl.Electric.BandaXarxa.T+2*parControl.Electric.BandaXarxa.Kp)/2, (parControl.Electric.BandaXarxa.Ki*parControl.Electric.BandaXarxa.T-2*parControl.Electric.BandaXarxa.Kp)/2];
parControl.Electric.BandaXarxa.Gcden=[1 -1];
parControl.Electric.BandaXarxa.SVPWM.T=parControl.Electric.Generador.T;
parControl.Electric.BandaXarxa.SVPWM.resolucio=64;


parControl.Mecanic.T=1e-2;
parControl.Mecanic.KCP=parMecanica.gammatopt/(parMecanica.nu^3*parMecanica.omegatn^2);
parControl.Mecanic.omegamn=parMecanica.omegatn*parMecanica.nu;
parControl.Mecanic.Kp=0.1;
parControl.Mecanic.Ki=0.02;

windspeed.string1=7;
windspeed.string2=7; 
windspeed.string3=7;

Length=1;
Rbase=43.56;
Lbase=0.1385;
Cbase=0.73e-6;
Vbase=66000;

%% Signal Generator
% Define signal generator here
eigenS = [0j, 7.2722e-5j, 0.0026j, 0.0150j, 0.6912j, 1.2566j];
eigenS_multiplicity = ones(1, size(eigenS, 2));

% Compute S and L
[S, L] = cal_QR(eigenS, eigenS_multiplicity);
L = L';
nu = size(S, 1);

% Initial condition SG
load w0.mat;


% Initial condition ROM
ksai0_value = 0.0;
ksai0 = ksai0_value * rand(nu, 1);

gain = 220 .* ones(nu, 1);











