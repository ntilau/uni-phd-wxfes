clear all; clc;
close all
Config();
Sys.pOrd = 4;
Sys.hOrd = 1;
prjName = 'WaveGuide';
% WriteWaveGuide(19.05,1.4*19.05);
WriteWaveGuide(22.86,22.86);
Mesh = IOrPoly(prjName, 'q34aA', Sys.hOrd, 1e-3);
Mesh.a = 22.86e-3;
Mesh.b = Mesh.a/2;
Sys.Height = Mesh.b;
Sys.WPnModes = 1;
Sys.WPportPlot = 1;
Sys.WPmodePlot = 1;
% Sys.Einc = 1e3;
Sys.WPpow = 1;
Mesh.BC.Dir = 1;
Mesh.BC.WP = [11 12];
% PlotMesh(Mesh,1)
[Sys,Mesh] = AssembLin(Sys, Mesh);
%%
Sys.nFreqs = 1;
Sys.freqs = 1e10;
% Sys.nFreqs = 21;
% Sys.freqs = linspace(5e9,15e9,Sys.nFreqs);
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
    if freq == 10e9
        Sys.uPlot = Sys.u;
    end
end

% %% magnetic field
% [Shape1, Shape1Deriv] = CalcShapeFunctions(1, Sys.pOrd);
% [Shape2, Shape2DerivX, Shape2DerivY] = CalcShapeFunctions(2, Sys.pOrd);
% Hconst = 1/(1i*2*pi*freq*Sys.mu0);
% Sys.u = zeros(Mesh.NELE,2);
% for ie=1:Mesh.NELE
%     [gIs] = CalcGlobIndex(2, Sys.pOrd, Mesh, ie);
%     [detJ, invJt] = CalcJacobian(Mesh.node(Mesh.ele(ie,:),:));
%     sol = Sys.uPlot(gIs);
%     dE = invJt*[Shape2DerivX(.5,.5); Shape2DerivY(.5,.5)]*sol;
%     Sys.u(ie,:) = (Hconst*[dE(2), -dE(1)]);
% end         
% IOwVTKH(Sys, Mesh, [prjName,'H']);

%%% Scattering
% figure; plot(Sys.freqs,Sys.db(Sys.Sparams(Sys.WPmodePlot:Sys.WPnModes:end,:).'));
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