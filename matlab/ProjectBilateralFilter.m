clear all; %clc;
close all
Config();
Sys.pOrd = 2;
Sys.hOrd = 1;
prjName = 'BilatFilter';
Mesh = IOrPoly(prjName, 'q34a1A', Sys.hOrd, 1e-6);
% Sys.Einc = 15e3 * sqrt(1.651e-03/2);
% Sys.Einc = 15e3;
Sys.Height = 1.651e-3/2;
Sys.WPnModes = 15;
Sys.WPportPlot = 1;
Sys.WPmodePlot = 1;
Sys.WPpow = 1;
Mesh.epsr = [1 2.1];% [1 2.1*(1-1i*0.02)];
Mesh.BC.Dir = 1;
Mesh.BC.WP = [11 12];
% PlotMesh(Mesh,1)

[Sys,Mesh] = AssembLin(Sys, Mesh);
%%
Sys.nFreqs = 1;
Sys.freqs = 146e9;
Sys.fPlot = Sys.freqs(1);
Sys.nFreqs = 81;
Sys.freqs = linspace(138e9,158e9,Sys.nFreqs);
Sys.Sparams = zeros(length(Mesh.BC.WP)*Sys.WPnModes, Sys.nFreqs);
for kf = 1:Sys.nFreqs
    freq = Sys.freqs(kf);
    fprintf('freq = %g GHz\n',freq/1e9);
    Sys =  AssembWP(Sys, freq);
    
    X = Sys.A\Sys.B;
    
    sp = X(1:length(Sys.WP)*Sys.WPnModes,1:length(Sys.WP)*Sys.WPnModes) - ...
        eye(length(Sys.WP)*Sys.WPnModes);
    Sys.Sparams(:,kf) = sp(:,(Sys.WPportPlot-1)*Sys.WPnModes+Sys.WPmodePlot);
    fprintf('  losses = %2.2g%%\n', (1 - norm(Sys.Sparams(:,kf)))*100);

    Sys.u = zeros(Sys.NDOFs,(Sys.WPportPlot-1)*Sys.WPnModes+Sys.WPmodePlot);
    Sys.u(Sys.nnWP) = X(length(Sys.WP)*Sys.WPnModes+1:end, ...
        (Sys.WPportPlot-1)*Sys.WPnModes+Sys.WPmodePlot);
    for ip=1:length(Sys.WP)
        Sys.u(Sys.WP{ip}) = Sys.WPgvec{ip}*X((1:Sys.WPnModes)+(ip-1)*Sys.WPnModes, ...
            (Sys.WPportPlot-1)*Sys.WPnModes+Sys.WPmodePlot);
    end
    if freq == Sys.fPlot
        Sys.uPlot = Sys.u;
    end
end
%%% Scattering
figure; plot(Sys.freqs*1e-9,Sys.db(Sys.Sparams(Sys.WPmodePlot:Sys.WPnModes:end,:).'));
axis([138 158 -40 0])
%%% Field
Sys.u = Sys.uPlot;
if exist('pdeplot','file')
    figure; pdeplot(Mesh.refNode.',[],Mesh.refEle.',...
        'xydata',(abs(Sys.u)),...
        'mesh','off',...
        'colormap','jet','xygrid','off');
    axis equal; axis tight; camlight left; lighting phong;
else
    IOwVTK(Sys, Mesh, prjName);
end

max(abs(Sys.u))