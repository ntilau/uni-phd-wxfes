clear; close all; % clc;
Config()
% profile on;
% db = @(x)20*log10(x);
functPlot = @(x)abs(x);
theta = 0;
Sys.pOrd = 3;
Sys.hOrd = 1;
Sys.k = 2*pi;%*1.5e9/299792458;
Sys.Z0 = 120*pi;
Sys.kEinc = [cosd(theta) sind(theta) 0];
Sys.toll = 0.001;


Mesh = IOrPoly('ModelScatteringDD', 'q34a0.001A', Sys.hOrd, 1);
% Mesh = IOrPoly('ModelScatteringSquare', 'q34aA', Sys.hOrd);
% PlotMesh(Mesh,1)
Mesh.BC.Dir = 1;
Mesh.BC.ABC = 133;
Mesh.epsr = [1 1];
%Mesh.mur = [1 1];
% Mesh.mur = {eye(2), [mur 1i*kr; -1i*kr mur]};

% PlotMesh(Mesh, 1)
% return

[Sys, Mesh] = AssembLin(Sys, Mesh);

Mesh.BC.DD = 13;
[Sys1, Mesh1] = AssembDD(Sys, Mesh, 1);
[Sys2, Mesh2] = AssembDD(Sys, Mesh, 2);

Sys.A = (Sys.S-Sys.k^2*Sys.T+1i*Sys.k*Sys.ABC);
Sys.B = (Sys.S-Sys.k^2*Sys.T);

if isfield(Sys,'Dir')
    Sys.A(Sys.Dir{1},:) = 0; Sys.A(:,Sys.Dir{1}) = 0;
    Sys.A(Sys.Dir{1},Sys.Dir{1}) = eye(length(Sys.Dir{1}));
    Sys.b = Sys.B(:,Sys.Dir{1}) * Sys.fsEinc(Sys.Dir{1}) ;
    Sys.b(Sys.Dir{1}) = - Sys.fsEinc(Sys.Dir{1});
    Sys.B(Sys.Dir{1},:) = 0; Sys.B(:,Sys.Dir{1}) = 0;
    Sys.B(Sys.Dir{1},Sys.Dir{1}) = eye(length(Sys.Dir{1}));
else
    Sys.b = - Sys.B * Sys.fsEinc +  Sys.f;
end

Sys1.Ddd = Sys1.DD(:,Sys1.DirDD);
Sys1.Tdd = 1i/Sys.k*Sys1.DD(Sys1.DirDD,Sys1.DirDD);
Sys1.A = [(Sys1.S-Sys1.k^2*Sys1.T+1i*Sys1.k*Sys1.ABC), Sys1.Ddd;...
    Sys1.Ddd.', Sys1.Tdd];
Sys1.B = blkdiag((Sys1.S-Sys1.k^2*Sys1.T), eye(length(Sys1.DirDD)));
Sys1.T12 = [zeros(Sys.NDOFs), zeros(Sys.NDOFs,length(Sys1.DirDD));...
    Sys1.Ddd.', -Sys1.Tdd];
% Sys1.T12 = blkdiag(zeros(Sys.NDOFs), Sys1.Tdd);
if isfield(Sys1,'Dir')
    Sys1.A(Sys1.Dir,:) = 0; Sys1.A(:,Sys1.Dir) = 0;
    Sys1.A(Sys1.Dir,Sys1.Dir) = eye(length(Sys1.Dir));
    Sys1.b = Sys1.B(:,Sys1.Dir) * Sys1.fEinc(Sys1.Dir) ;
    Sys1.b(Sys1.Dir) = - Sys1.fEinc(Sys1.Dir);
    Sys1.B(Sys1.Dir,:) = 0; Sys1.B(:,Sys1.Dir) = 0;
    Sys1.B(Sys1.Dir,Sys1.Dir) = eye(length(Sys1.Dir));
%     Sys1.T12(Sys1.Dir,:) = 0; Sys1.T12(:,Sys1.Dir) = 0;
%     Sys1.T12(Sys1.Dir,Sys1.Dir) = eye(length(Sys1.Dir));
else
    Sys1.b = - Sys1.B * Sys1.fEinc +  Sys1.f;
end
idx1 = [Sys1.DirReg; Sys.NDOFs+(1:length(Sys1.DirDD)).'];
Sys1.A = Sys1.A(idx1,idx1);
Sys1.b = Sys1.b(idx1,1);

% figure; spy(Sys1.T12)
Sys2.Ddd = Sys2.DD(:,Sys2.DirDD);
Sys2.Tdd = 1i/Sys.k*Sys2.DD(Sys2.DirDD,Sys2.DirDD);
Sys2.A = [(Sys2.S-Sys2.k^2*Sys2.T+1i*Sys2.k*Sys2.ABC), Sys2.Ddd;...
    Sys2.Ddd.', Sys2.Tdd];
Sys2.B = blkdiag((Sys2.S-Sys2.k^2*Sys2.T), eye(length(Sys2.DirDD)));
Sys2.T21 = [zeros(Sys.NDOFs), zeros(Sys.NDOFs,length(Sys2.DirDD));...
    Sys2.Ddd.', -Sys2.Tdd];
% Sys2.T21 = blkdiag(zeros(Sys.NDOFs), Sys2.Tdd);
if isfield(Sys2,'Dir')
    Sys2.A(Sys2.Dir,:) = 0; Sys2.A(:,Sys2.Dir) = 0;
    Sys2.A(Sys2.Dir,Sys2.Dir) = eye(length(Sys2.Dir));
    Sys2.b = Sys2.B(:,Sys2.Dir) * Sys2.fEinc(Sys2.Dir) ;
    Sys2.b(Sys2.Dir) = - Sys2.fEinc(Sys2.Dir);
    Sys2.B(Sys2.Dir,:) = 0; Sys2.B(:,Sys2.Dir) = 0;
    Sys2.B(Sys2.Dir,Sys2.Dir) = eye(length(Sys2.Dir));
%     Sys2.T21(Sys2.Dir,:) = 0; Sys2.T21(:,Sys2.Dir) = 0;
%     Sys2.T21(Sys2.Dir,Sys2.Dir) = eye(length(Sys2.Dir));
else
    Sys2.b = - Sys2.B * Sys2.fEinc +  Sys2.f;
end
idx2 = [Sys2.DirReg; Sys.NDOFs+(1:length(Sys2.DirDD)).'];
Sys2.A = Sys2.A(idx2,idx2);
Sys2.b = Sys2.b(idx2,1);

% figure;plot(Mesh.refNode(Sys1.DirReg,1),Mesh.refNode(Sys1.DirReg,2),'.'); axis equal
% figure;plot(Mesh.refNode(Sys2.DirReg,1),Mesh.refNode(Sys2.DirReg,2),'.'); axis equal
% return
% figure; spy(Sys.A)
% figure; spy(Sys1.A)
% figure; spy(Sys2.A)
% return

tic
Sys.u = Sys.A\Sys.b;
fprintf('Direct solver: %g s\n',toc);
Sys.u = Sys.fsEinc + Sys.u(1:Sys.NDOFs);

% Sys1.u = Sys1.A\Sys1.b;
% Sys2.u = Sys2.A\Sys2.b;

Sys1.uf = zeros(Sys.NDOFs,1);
Sys2.uf = zeros(Sys.NDOFs,1);
Sys1.u = zeros(length(Sys1.DirReg)+length(Sys1.DirDD),1);
Sys2.u = zeros(length(Sys2.DirReg)+length(Sys1.DirDD),1);
Sys1.uft = zeros(Sys.NDOFs,1);
Sys2.uft = zeros(Sys.NDOFs,1);

% Sys1.uf(Sys1.DirReg,1) = Sys1.u(1:length(Sys1.DirReg));
% Sys2.uf(Sys2.DirReg,1) = Sys2.u(1:length(Sys2.DirReg));

Sys1.T12 = Sys1.T12(idx1,idx2);
Sys2.T21 = Sys2.T21(idx2,idx1);

% figure(1); pdeplot(Mesh.refNode.',[],Mesh.refEle.',...
%     'xydata',(functPlot(Sys.u)),...
%     'mesh','off',...
%     'colormap','jet','xygrid','on');
% axis equal; camlight left; lighting phong; 

% Sys1.T21 = Sys1.T21(idx2,:);
% size(Sys1.T12,2)
% size(Sys2.T21,2)
error=1;
i=1;
while error > Sys.toll
    Sys1.u = Sys1.A\(Sys1.b+Sys1.T12*Sys2.u);
    Sys2.u = Sys2.A\(Sys2.b+Sys2.T21*Sys1.u);
    Sys1.ufp = Sys1.uf;
    Sys2.ufp = Sys2.uf;
    Sys1.uf(Sys1.DirReg,1) = Sys1.u(1:length(Sys1.DirReg));
    Sys2.uf(Sys2.DirReg,1) = Sys2.u(1:length(Sys2.DirReg));
    error1 = norm( Sys1.uf - Sys1.ufp ) / norm( Sys1.uf);
    error2 = norm( Sys2.uf - Sys2.ufp ) / norm( Sys2.uf);
    error(i) = max(error1, error2);
    fprintf('%g\n', error(i));
    i=i+1;
%     break
%     Sys1.uft(Sys1.DirReg,1) = Sys.fEinc(Sys1.DirReg) + Sys1.u(1:length(Sys1.DirReg));
%     Sys2.uft(Sys2.DirReg,1) = Sys.fEinc(Sys2.DirReg) + Sys2.u(1:length(Sys2.DirReg));
%     Sys1.uft(Sys1.DirDD) = 0;
% %     Sys2.uf(Sys2.DirDD) = 0;
%     uDD = Sys1.uft + Sys2.uft;
%     figure(2); pdeplot(Mesh.refNode.',[],Mesh.refEle.',...
%         'xydata',(functPlot(uDD)),...
%         'mesh','off',...
%         'colormap','jet','xygrid','on');
%     axis equal; camlight left; lighting phong; 
end
fprintf('No of iterations = %d\n', i-1);

% return
%%
Sys1.uf(Sys1.DirReg,1) = Sys.fsEinc(Sys1.DirReg) + Sys1.u(1:length(Sys1.DirReg));
Sys2.uf(Sys2.DirReg,1) = Sys.fsEinc(Sys2.DirReg) + Sys2.u(1:length(Sys2.DirReg));
% Sys1.uf = Sys1.uf + Sys.fEinc;
% Sys2.uf = Sys2.uf + Sys.fEinc;

Sys1.uf(Sys1.DirDD) = 0;
% Sys1.uf(Sys.Dir{1}DD) = 0;
uDD = Sys1.uf + Sys2.uf;
Sys.u = uDD;

% IOwVTK(Sys, Mesh, 'Scattering')
%%% Field
% Sys.u = Sys.uPlot;
if exist('pdeplot','file')
    figure; pdeplot(Mesh.refNode.',[],Mesh.refEle.',...
        'xydata',(abs(Sys.u)),...
        'mesh','off',...
        'colormap','jet','xygrid','off');
    axis equal; axis tight; camlight left; lighting phong;
else
    IOwVTK(Sys, Mesh, prjName);
end

return
figure; pdeplot(Mesh.refNode.',[],Mesh.refEle.',...
    'xydata',(functPlot(Sys.u)),...
    'mesh','off',...
    'colormap','jet','xygrid','on');
axis equal; camlight left; lighting phong; 

figure; pdeplot(Mesh.refNode.',[],Mesh.refEle.',...
    'xydata',(functPlot(uDD)),...
    'mesh','on',...
    'colormap','jet','xygrid','off');
axis equal; camlight left; lighting phong;

fprintf('Total error = %2.4g\n', norm(Sys.u-uDD)/norm(Sys.u));
figure; semilogy(error(2:end)); axis tight

return

figure; pdeplot(Mesh.refNode.',[],Mesh.refEle.',...
    'xydata',(functPlot(Sys1.uf)),...
    'zdata',(functPlot(Sys1.uf)),'mesh','off',...
    'colormap','jet','xygrid','on');
axis equal; camlight left; lighting phong; 

figure; pdeplot(Mesh.refNode.',[],Mesh.refEle.',...
    'xydata',(functPlot(Sys2.uf)),...
    'zdata',(functPlot(Sys2.uf)),'mesh','off',...
    'colormap','jet','xygrid','on');
axis equal; camlight left; lighting phong; 

% IOwVTK(Sys, Mesh, 'Scattering')
% return
% profile viewer;
% profile off; 
return
%%
% figure; pdeplot(Mesh.refNode.',[],Mesh.refEle.',...
%     'xydata',(functPlot(Sys.u)),...
%     'mesh','off',...
%     'colormap','jet','xygrid','on');
% axis equal; camlight left; lighting phong; 
% return
figure; pdeplot(Mesh.refNode.',[],Mesh.refEle.',...
    'xydata',(functPlot(Sys.u)),...
    'zdata',(functPlot(Sys.u)),'mesh','off',...
    'colormap','jet','xygrid','on');
axis equal; camlight left; lighting phong; 
%%
return
tic;
[M1,M2] = ilu(Sys.A,struct('type','ilutp','droptol',1e-3));
fprintf('ILU preconditioner: %g s\n',toc);
%%
% tic
% [Sys.u2,flag,relres,iter,resvec]  = qmr(Sys.A,Sys.b,1e-6,1000);
% fprintf('qmr: %g s\n',toc);

% tic
% [Sys.u2,flag,relres,iter,resvec]  = gmres(Sys.A,Sys.b,[],1e-6,1000);
% fprintf('gmres: %g s\n',toc);
% figure; semilogy(resvec)

tic
[Sys.u2,flag,relres,iter,resvec]  = gmres(Sys.A,Sys.b,[],1e-11,1000,M1,M2);
fprintf('GMRES Solver: %g s\n',toc);

Sys.u2 = Sys.fEinc + Sys.u2;
fprintf('Relative error: %2.4g\n', norm(Sys.u-Sys.u2)/norm(Sys.u));
% Sys.u2 = Sys.u2 - Sys.u;

% figure; pdeplot(Mesh.refNode.',[],Mesh.refEle.',...
%     'xydata',(functPlot(Sys.u2)),...
%     'zdata',(functPlot(Sys.u2)),'mesh','off',...
%     'colormap','jet','xygrid','on');
% axis equal; camlight left; lighting phong

figure; semilogy(resvec,'.-')
% profile viewer;
% profile off;
