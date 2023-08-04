clear all; clc;
% close all
Config();
Sys.pOrd = 2;
Sys.hOrd = 1;
a = 19.05;
b = a/2;
r = 0.05*a;
c = 0.4*a;
L = 1.4*a;
EpsR = 112.5;
nt = 60;
WriteTwoPosts(a,c,r,L,nt);
prjName = 'TwoPosts';
Mesh = IOrPoly(prjName, 'q34a1A', Sys.hOrd, 1e-3);
prjName = [ prjName, 'NL'];
Sys.WPnModes = 1;
Sys.WPportPlot = 1;
Sys.WPmodePlot = 1;
Sys.WPpow = 1;
Sys.Einc = 10e3 * sqrt(1.651e-03/2);
Mesh.epsr = [1 EpsR];
Mesh.kerr = [0 1.625e-10];
Mesh.BC.Dir = 1;
Mesh.BC.WP = [11 12];
% PlotMesh(Mesh,1)

% [Sys,Mesh] = AssembHB(Sys, Mesh);
Sys = CalcDoFsNumber(Sys, Mesh);
Sys.u = zeros(Sys.NDOFs,1);

%%
fc = 7.868577546294576e+09;
Sys.nFreqs = 1;
Sys.freqs = 1.43*fc;
Sys.fPlot = Sys.freqs(1);
% Sys.nFreqs = 81;
% Sys.freqs = linspace(1.38*fc, 1.48*fc, Sys.nFreqs);
Sys.Sparams = zeros(length(Mesh.BC.WP)*Sys.WPnModes, Sys.nFreqs);
for kf = 1:Sys.nFreqs
    freq = Sys.freqs(kf);
    fprintf('freq = %g GHz\n',freq/1e9);
    error = 1;
    Sys.u0 = Sys.u;
    while error > 1e-4
        [Sys,Mesh] = AssembNL(Sys, Mesh);
        Sys =  AssembWP(Sys, freq);

        X = Sys.A\Sys.B;

        sp = X(1:length(Sys.WP)*Sys.WPnModes,1:length(Sys.WP)*Sys.WPnModes) - ...
            eye(length(Sys.WP)*Sys.WPnModes);
        Sys.Sparams(:,kf) = sp(:,(Sys.WPportPlot-1)*Sys.WPnModes+Sys.WPmodePlot);
%         fprintf('  losses = %2.2g%%\n', (1 - norm(Sys.Sparams(:,kf)))*100);
        
        Sys.u = zeros(Sys.NDOFs,(Sys.WPportPlot-1)*Sys.WPnModes+Sys.WPmodePlot);
        Sys.u(Sys.nnWP) = X(length(Sys.WP)*Sys.WPnModes+1:end, ...
            (Sys.WPportPlot-1)*Sys.WPnModes+Sys.WPmodePlot);
        for ip=1:length(Sys.WP)
            Sys.u(Sys.WP{ip}) = Sys.WPgvec{ip}*X((1:Sys.WPnModes)+(ip-1)*Sys.WPnModes, ...
                (Sys.WPportPlot-1)*Sys.WPnModes+Sys.WPmodePlot);
        end
        error = norm(Sys.u - Sys.u0) / norm(Sys.u);
        Sys.u0 = Sys.u;
        fprintf(' %g\n',error);
    end
    if freq == Sys.fPlot
        Sys.uPlot = Sys.u;
    end
end

%%% Scattering
figure; plot(Sys.freqs,Sys.db(Sys.Sparams(Sys.WPmodePlot:Sys.WPnModes:end,:).'));
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

