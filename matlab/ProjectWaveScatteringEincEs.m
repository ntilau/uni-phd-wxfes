clear;% close all; % clc;
Config()
% profile on;
% db = @(x)20*log10(x);
functPlot = @(x)abs(x);
theta = 0;
Sys.pOrd = 2;
Sys.hOrd = 1;
Sys.k = 2*pi;%*1.5e9/299792458;
Sys.Z0 = 120*pi;
Sys.kEinc = [cosd(theta) sind(theta) 0];
FEBI = false;
% npts = 10;
% Mesh = BuildRegularSquare(npts);
% Mesh.BC.ABC = 1;
% Mesh.node = Mesh.node*2-1;
% Sys.NDOF = CalcDoFsNumber(Sys.pOrd,Mesh.NNODE,Mesh.NELE,Mesh.NSPIG);
% Mesh = IOrPoly('ModelScatteringVoid', 'q34a0.005AQ');
% Mesh.BC.ABC = 133;

%Mesh = IOrPoly('ModelScatteringDD', 'q34aA', Sys.hOrd);

Mesh = IOrPoly('ModelScatteringSquare', 'q34a0.01A', Sys.hOrd, 1);
Mesh.BC.Dir = 1;
Mesh.BC.DD = 13;

Mesh.BC.ABC = 133;
Mesh.epsr = [1 1];
Mesh.mur = [1 2];


% PlotMesh(Mesh, 1)
% return

[Sys Mesh] = AssembLin(Sys, Mesh);

Sys.A = (Sys.S-Sys.k^2*Sys.T+1i*Sys.k*Sys.ABC);
Sys.B = (Sys.S-Sys.k^2*Sys.T);

if isfield(Sys,'Dir')
    Sys.A(Sys.Dir{1},:) = 0; Sys.A(:,Sys.Dir{1}) = 0;
    Sys.A(Sys.Dir{1},Sys.Dir{1}) = eye(length(Sys.Dir{1}));
    Sys.b = Sys.B(:,Sys.Dir{1}) * Sys.fEinc(Sys.Dir{1}) ;
    Sys.b(Sys.Dir{1}) = - Sys.fEinc(Sys.Dir{1});
%     Sys.b = - Sys.B * Sys.fEinc + Sys.f;
%     Sys.b(Sys.Dir{1}) = - Sys.fEinc(Sys.Dir{1});
    Sys.B(Sys.Dir{1},:) = 0; Sys.B(:,Sys.Dir{1}) = 0;
    Sys.B(Sys.Dir{1},Sys.Dir{1}) = eye(length(Sys.Dir{1}));
%     Sys.b = Sys.b + Sys.f;
else
    Sys.b = - Sys.B * Sys.fEinc +  Sys.f;
end

if FEBI
    ns = length(Sys.DirABC);
    Sys.Z = zeros(ns,ns);
    xy = Mesh.node(Sys.DirABC,:).';
    ds = zeros(ns,1);
    for m = 1:ns
        idSpigs = find(Mesh.spig2(:,1) == Sys.DirABC(m) | Mesh.spig2(:,2) == Sys.DirABC(m));
        tmp_i=1;
        for i=1:length(idSpigs)
            tmp = find(Sys.idsABC == idSpigs(i));
            if tmp
                id2Spigs(tmp_i) = tmp;
                tmp_i = tmp_i + 1;
            end
        end
        ds(m) = norm(diff(Mesh.node(Mesh.spig2(id2Spigs(1),:),:)))/2 +  ...
            norm(diff(Mesh.node(Mesh.spig2(id2Spigs(2),:),:)))/2;
      for n = m+1:ns
        tmp = Sys.k*Sys.Z0*ds(m)/4* ...
          besselh(0,2,Sys.k*sqrt((xy(1,m)-xy(1,n))^2+(xy(2,m)-xy(2,n))^2));
        Sys.Z(m,n) = tmp;
        Sys.Z(n,m) = tmp;
      end
    end
    for n = 1:ns
      Sys.Z(n,n) = Sys.k*Sys.Z0*ds(n)/4*(1-1i*2/pi*log(Sys.k*1.781*ds(n)/4/exp(1)));
    end
    Sys.b = greenFFopt(Sys.k, xy(1,:), xy(2,:), 0).';
%     Sys.C = sparse(Sys.NDOF, ns);
    Sys.C = sparse(Sys.DirABC, 1:ns, Sys.f(Sys.DirABC), Sys.NDOF, ns );

    Sys.Mat = [Sys.A, Sys.C; Sys.C.', Sys.Z];
    Sys.Rhs = [zeros(Sys.NDOF,1); Sys.b]; 
    sol = Sys.Mat\Sys.Rhs;
    Sys.u = sol(1:Sys.NDOF);
    Sys.J = sol(Sys.NDOF+1:Sys.NDOF+ns);


    phi2 = linspace(-pi,pi,1001);
    Es = greenFFopt(Sys.k, xy(1,:), xy(2,:), phi2)*(Sys.J.*ds);
    figure;plot(phi2*180/pi,10*log10(Sys.k*Sys.Z0^2/4*abs(Es).^2)); axis tight
else
    tic
    Sys.u = Sys.A\Sys.b;
    fprintf('Direct solver: %g s\n',toc);
%     plot(abs(Sys.u(Sys.NDOF+1:end)))
    Sys.u = Sys.fEinc + Sys.u(1:Sys.NDOFs);
end

% IOwVTK(Sys, Mesh, 'Scattering')
% return
% profile viewer;
% profile off; 
% return
%%
figure; pdeplot(Mesh.refNode.',[],Mesh.refEle.',...
    'xydata',(functPlot(Sys.u)),...
    'zdata',(functPlot(Sys.u)),...
    'mesh','off',...
    'colormap','jet','xygrid','on');
axis equal; camlight left; lighting phong; 
return
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
