%% test coupling 
clear all;
Config();
format longg
clear all; close all;% clc
Sys.HBharms = 1:3;
Sys.freq = 1;
% Sys.OverSampling = 4;

Mtrl = GetMtrlParamsFFT(Sys);

int = 1;
idx = 1:2:max(Sys.HBharms);
E = int*ones(floor(Mtrl.nHarms/2)+1,1);
a = 10;
b = 5;
c = 1;

%%% generic implementation
Mtrl.func = @(x)(a+b*x+c*x.^2);
% Mtrl.func = @(x)(int./(int+x.^2));
Field.E = int*ones(2*Mtrl.nHarms,1);
Field.E(2:2:end) = 0;
Field.Sin = sin(2*pi*Sys.freq*Sys.HBharms.'*Mtrl.t);
Field.Cos = cos(2*pi*Sys.freq*Sys.HBharms.'*Mtrl.t);
Mtrl.Sin = sin(2*pi*Mtrl.fnl.'*Mtrl.t);
Mtrl.Cos = cos(2*pi*Mtrl.fnl.'*Mtrl.t);

% Dref = GetCouplKerrAnalyt(length(Sys.HBharms(idx)), E, a, b);
% Dsin = GetCouplKerrFFT(length(Sys.HBharms(idx)), E, a, b, Mtrl.Cos, Field.Sin(idx,:), ...
%     Mtrl.posR(idx), Mtrl.posL(idx));
[Mtrl, Field] = GetCouplFFT(Mtrl, Field);
% disp(round(Mtrl.Dref*1e2)*1e-2)
% disp(round(Mtrl.Dsin*1e2)*1e-2)
% disp(round(Mtrl.D*1e2)*1e-2)
% disp(round(Mtrl.Dref))
% disp(round(Mtrl.Dsin))
% disp(round(Mtrl.D*1e6)*1e-6)

disp(Mtrl.D*Field.E)

% disp(D-Dref)


