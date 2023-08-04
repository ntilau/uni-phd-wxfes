clear all; %clc;
close all
Config();
Sys.pOrd = 3;
Sys.hOrd = 1;
prjName = 'BilatFilter';
Mesh = IOrPoly(prjName, 'q34aA', Sys.hOrd, 1e-6);
prjName = [ prjName, 'HB'];
Sys.WPnModes = 15;
Sys.WPportPlot = 1;
Sys.WPmodePlot = 1;
% Sys.WPpow = 1;
Sys.Height = 1.651e-03/2;
Sys.Einc = 10e3;
Mesh.epsr = [1 2.1];%*(1-1i*0.02)];
Mesh.kerr = [0 1.625e-10];
Mesh.NLlab = 2;
Mesh.BC.Dir = 1;
Mesh.BC.WP = [11 12];
Sys.HBharmPlot = 1;
Sys.HBharms = [1 3 5 7];
Sys.SinOnly = true;
Sys.nHarms = length(Sys.HBharms);
%%% comment Sys.FFT to go analytic odd frequencies + kerr ()
% Sys.FFT = true;
Sys.OverSampling = 2;
if ~isfield(Sys,'SinOnly')
    Sys.nHarms = Sys.nHarms*2;
else
    if ~Sys.SinOnly
        Sys.HBharmPlot = 2*Sys.HBharmPlot-1;
    end
end
% PlotMesh(Mesh,1)
Sys.WPnum = 10;
Sys.Pfund = 1;
Sys.Pitrf = 0;
Sys = CalcDoFsNumber(Sys, Mesh);
Sys.u = zeros(Sys.NDOFs*Sys.nHarms,1);

%%
Sys.nFreqs = 1;
Sys.freqs = 144e9;
% Sys.nFreqs = 21;
% Sys.freqs = linspace(138e9,158e9,Sys.nFreqs);
Sys.Sparams = zeros(length(Mesh.BC.WP)*Sys.WPnModes*Sys.nHarms, Sys.nFreqs);
for kf = 1:Sys.nFreqs
    Sys.freq = Sys.freqs(kf);
    fprintf('freq = %g GHz\n',Sys.freq/1e9);
    error = 1;
    Sys.u0 = Sys.u;
    while error > 1e-9
        if isfield(Sys, 'FFT')
            [Sys,Mesh] = AssembHBKerr(Sys, Mesh);
        else
            [Sys,Mesh] = AssembHBKerrAnalytic(Sys, Mesh);
        end
        Sys =  AssembWPHB(Sys);
%         Sys.B(:,1) = Sys.B(:,1) + Sys.B(:,Sys.WPnModes+1);
        X = Sys.A\Sys.B;

        sp = X(1:length(Sys.WP)*Sys.WPnModes*Sys.nHarms,...
            1:length(Sys.WP)*Sys.WPnModes*Sys.nHarms) - ...
            eye(length(Sys.WP)*Sys.WPnModes*Sys.nHarms);
        Sys.Sparams(:,kf) = ...
            sp(:,(Sys.WPportPlot-1)*Sys.WPnModes*Sys.nHarms+...
            Sys.WPmodePlot);
        fprintf('  losses = %2.2g%%\n', abs(1 - norm(Sys.Sparams(:,kf)))*100);
        
        Sys.u = zeros(Sys.NDOFs*Sys.nHarms,1);
        Sys.u(Sys.nnWP) = X(length(Sys.WP)*Sys.WPnModes*Sys.nHarms+1:end, ...
            (Sys.WPportPlot-1)*Sys.WPnModes*(Sys.HBharmPlot*Sys.nHarms)+Sys.WPmodePlot);
        for jh=1:length(Sys.Harms)
            for ip=1:length(Sys.WP)
                Sys.u((jh-1)*Sys.NDOFs+Sys.WP{ip}) = ...
                    Sys.WPgvec{ip,jh}*X((1:Sys.WPnModes)+((jh-1)+(ip-1)*Sys.nHarms)*Sys.WPnModes, ...
                    (Sys.WPportPlot-1)*Sys.WPnModes+Sys.WPmodePlot);
            end
        end
        
        error = norm(Sys.u - Sys.u0) / norm(Sys.u);
        Sys.u0 = Sys.u;
        fprintf(' %g\n',error);
    end
    if Sys.freq == Sys.freqs
        Sys.uPlot = Sys.u;
    end
end

%%% Scattering
if Sys.nFreqs == 1
    Sys.db(Sys.Sparams(Sys.WPmodePlot:Sys.WPnModes:end,:))
else
    figure; plot(Sys.freqs,Sys.db(Sys.Sparams(Sys.WPmodePlot:Sys.WPnModes:end,:).'));
end

%%% Field
for HBharmPlot = 1:3%length(Sys.HBharms)
    Sys.u = Sys.uPlot((HBharmPlot-1)*Sys.NDOFs+(1:Sys.NDOFs));
%     if exist('pdeplot','file')
%         figure; pdeplot(Mesh.refNode.',[],Mesh.refEle.',...
%             'xydata',(abs(Sys.u)),...
%             'mesh','off',...
%             'colormap','jet','xygrid','off');
%         axis equal; axis tight; camlight left; lighting phong;
%     else
        IOwVTK(Sys, Mesh, [prjName, num2str(HBharmPlot)]);
%     end
end

% max(abs(Sys.u))