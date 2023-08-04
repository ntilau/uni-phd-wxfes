%% test coupling 
Config();
format shorte
clear all; close all;% clc
Sys.HBharms = 1:2:3;
Sys.freq = 1;

Mtrl = GetMtrlParamsFFT(Sys);

E = 1e3*ones(Mtrl.nHarms,1);
a = 10;
b = 1;

%%% generic implementation
Mtrl.func = @(x)(a+b*x.^2);
Field.E = 1e3*ones(Mtrl.nHarms,1);
Field.Sin = sin(2*pi*Sys.freq*Sys.HBharms.'*Mtrl.t);
Field.Cos = cos(2*pi*Sys.freq*Sys.HBharms.'*Mtrl.t);
Mtrl.Sin = sin(2*pi*Mtrl.fnl.'*Mtrl.t);
Mtrl.Cos = cos(2*pi*Mtrl.fnl.'*Mtrl.t);


Dref = GetCouplKerrAnalyt(Mtrl.nHarms, E, a, b);
Dsin = GetCouplKerrFFT(Mtrl.nHarms, E, a, b, Mtrl.Cos, Field.Sin, ...
    Mtrl.posR, Mtrl.posL);
D = GetCouplFFT(Mtrl, Field)

% disp(D-Dref)


