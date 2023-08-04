clear all;% clc;
close all
% warning('off','all');
log = sprintf('#Begin: %d/%d/%d, %d:%d:%.2g\n',clock);
a = 0.9*25.4; b = 0.4*25.4;
% a = 10.1;
% b = 5.1;
% a = 1.5; b = 1;1
epsr = 1;
npts = 6;
Sys.pOrd = 2;
Sys.hOrd = 1;
TM11 = false;
% prjName = 'ModelWR90';
prjName = 'ModelWR90Strip';
% Mesh = BuildRegularSquare(5,3);
% Mesh = BuildRegularSquare(202,102);
% Mesh.node = Mesh.node*[a 0; 0 b]*1e-3;
Mesh = IOrPoly('ModelWR90Strip', 'q34a0.00001A', Sys.hOrd, 1);
% Mesh = IOrPoly( prjName, 'q34a0.00001A', Sys.hOrd);

% PlotPoly(prjName, figure, 1)
% return

% Mesh.node = Mesh.node*1000;
nmodes = 1;

Mesh.BC.Dir = 1;
% PlotMesh(Mesh,1)
% save ACESMesh Mesh;
% load ACESMesh;
% return
Sys.Hcurl = 1;
[Sys,Mesh] = AssembLin(Sys, Mesh);

% if TM11
    idx = 1:length(Sys.T);
    idx(Sys.Dir{1}) = 0;
    idx = find(idx);
%     Tte = Sys.T(idx,idx);
%     Ste = Sys.S(idx,idx);
% else
    Tte = Sys.T;
    Ste = Sys.S;
% end
pOrd = Sys.pOrd;
% return
% TM modes
% switch pOrd 
%   case 1
%     idx = find(nlab ~= 1);
%   case 2
%     idx = [find(nlab ~= 1) NNODE+find(slab ~= 1)]; % non Dirichlet nodes
%   case 3
%     idx = [find(nlab ~= 1) NNODE+find(slab ~= 1) ...
%       NNODE+NSPIG+find(slab ~= 1) NNODE+2*NSPIG+1:NNODE+2*NSPIG+NELE]; % non Dirichlet nodes
% end
% Stm = Ste(idx,idx);
% Ttm = Tte(idx,idx);
%% solution of eigenvalue problem
% N=2; toll= 1e-7;
% c0  =   299792458;
% mu0 =   pi*4e-7;
% Z0  =   c0*mu0;
% geom = 1; % 1 rect, 2 single ridge, 3 double ridge
% fieldLines = 0; 
% mesh = 0;
% coeff = 0; % S/b
% itp = 0; % number of modes to plot
% Vtm = zeros(NDOF,N);
% [Vte,Dte] = eig(full(Ste),full(Tte));
scal = 1;%(max(max(Tte))/min(min(Tte(Tte>0))));


% opts.tol = 1e-9;
[e] = eigs(Ste,epsr*Tte,2*nmodes,(2*pi*5e9/3e8)^2);
e = e(abs(e)>1e-5);
val0 = sqrt(e(1:nmodes))*3e8/2/pi
if true
    if Sys.pOrd == 1
        dir = find(~Mesh.slab);
    elseif Sys.pOrd == 2
        dir = [find(~Mesh.slab); Mesh.NSPIG+find(~Mesh.slab); (2*Mesh.NSPIG+1:2*Mesh.NSPIG+2*Mesh.NELE).'];
    elseif Sys.pOrd == 3
        dir = [find(~Mesh.slab); Mesh.NSPIG+find(~Mesh.slab); ...
            (2*Mesh.NSPIG+1:2*Mesh.NSPIG+2*Mesh.NELE).';...
            2*Mesh.NSPIG+2*Mesh.NELE+find(~Mesh.slab); ...
            (3*Mesh.NSPIG+2*Mesh.NELE+1:3*Mesh.NSPIG+6*Mesh.NELE).'];
    end
    dirn = 1:Sys.NDOFs;
    dirn(Sys.Dir{1})=0; %find(~Mesh.nlab);
    %idx(Sys.Dir) = 0;
    dirn = find(dirn);
    St = (Sys.St(dir,dir));
    Tt = (Sys.Tt(dir,dir));
    G = (Sys.G(dir,dirn)).';
    Sz = (Sys.S(dirn,dirn));
    Tz = (Sys.T(dirn,dirn));
    e = eigs(St,epsr*Tt,2*nmodes,(2*pi*5e9/3e8)^2);
    e = e(abs(e)>1e-5);
    val1 = sqrt(e(1:nmodes))*3e8/2/pi
end
m=1; n=0;
fc = sqrt((m*pi/a).^2+(n*pi/b).^2)*3e11/2/pi;
disp(abs(val0(end)-fc)/fc)
disp(abs(val1(end)-fc)/fc)

%%
nmodes = 4;
ak = linspace(3e9,20e9,51)*2*pi/3e8;
g = [];
B = Tt;
C = G;
opts.disp = 0;
tic
for i=1:length(ak)
    k=ak(i)*sqrt(epsr);
    A = St-k^2*Tt;
    D = Sz-k^2*Tz;
    if 0
        Sf = [B C'; C D];
        Tf = [B+A/(k^2) C'; C D];
%         Sf = Sf+1i*eps*Sf;
%         Tf = Tf+1i*eps*Tf;

%         e = eigs(sparse(Sf),sparse(Tf),nmodes,'lm');
        e = eigs(Tf\Sf,nmodes,'LM',opts);
        e = real(sqrt(k^2*(1-1./e)));
    else
        Sf = sparse(blkdiag(A,zeros(Sys.NDOFs-length(Sys.Dir{1}))));
        Tf = sparse([B C'; C D]);
        %Sf = Sf+1i*eps*Sf*100;
        %Tf = Tf+1i*eps*Tf*100;
        %e1 = eig(full(Sf),full(Tf));
        [vec,e] = eigs(Tf\Sf,nmodes,'SR',opts);
        %e = eigs(sparse(Sf),sparse(Tf),nmodes,'LM');
        %e1 = real(diag(e1(abs(e1))>1e-11));
        %e = real(sqrt(-e))
        %e1 = sqrt(-e1);
        %e = real(e) - imag(e);
        %e1 = real(e1) - imag(e1);
        %e1 = sort(e1,'descend');
        
        %[vec,e] = eigs(Sf,Tf,nmodes,'LM',opts);
        % e = real(sqrt((1/e + k^2)));
        e = sort(diag(real(sqrt(-e))),'descend');
    end
    g(:,i)= e;
%     g1(:,i)= e1(1);
end
time = toc/length(ak)


figure
    plot(ak*3e8/2*pi,g(1:nmodes,:).');
    %,'-', ak*3e8/2*pi,g1(1:nmodes,:).',':'); 
    axis tight
% disp(g)



return
%%
% k=sort(sqrt(diag(Dte))); k=k(2);
ks=sort(sqrt(diag(Ds)));
ks=ks(ks>0.001);% ks=ks(1);
if TM11
    m=1;n=1;
else
    m=1;n=0;
end
kref = sqrt((m*pi/a).^2+(n*pi/b).^2);
fprintf('ks: %2.16g\n',ks);
fprintf('kref-ks: %2.4g\n',kref-ks(1));
fprintf('rel err: %2.4g\n',(ks(1)-kref)/kref);
disp(ks*3e8/2/pi)
% k




return
[Vtm(idx,:),Dtm] = eigs(Stm,Ttm,N,'sr');
%%% remove trivial solutions
Vte = Vte( :,diag(Dte)>toll);
Dte = diag(Dte(Dte>toll));
Vtm = Vtm(:,diag(Dtm)>toll);
Dtm = diag(Dtm(Dtm>toll));
type = [0*ones(1,length(Dte)), 1*ones(1,length(Dtm))];
V = [Vte, Vtm];
D = blkdiag(Dte, Dtm);
f = sqrt(diag(D))*c0/2/pi;
kFem = sort(sqrt(diag(D)));
[f,ord] = sort(f);
V = V(:,ord);
type = type(1,ord);
%%
% close all;
if geom == 2 && fieldLines
  itp = 1;
  for i=1:length(D)
    if itp
      line(pts(:,1),pts(:,2),'color','k');
      hold on
      pdeplot(xy,[],ele(2:4,:),'xydata',V(:,i),'xystyle','off','contour',...
        'on','levels',50,'colorbar','off','mesh','off');
      title('Field lines of fundamental mode');
      axis equal;axis tight; 
      axis off
      itp = itp - 1;
    end
  end
  break
end

for i=1:length(D)
  if itp
    figure(1);
    line(pts(:,1),pts(:,2),'color','k');
    hold on;
    if type(i) == 1
      [Fxy(:,1),Fxy(:,2)] = pdegrad(xy,ele(2:4,:),V(:,i));
      Ei = pdeprtni(xy,ele(2:4,:),sqrt(sum(Fxy.^2,2)).');
      pdeplot(xy,[],ele(2:4,:),'xydata',Ei,...
        'flowdata',Fxy,'flowstyle','arrow','colorbar','off');
      title(['TM ',sprintf('@ %02.5g GHz\n',f(i)*1e-9)])
    else
      [Fxy(:,2),Fxy(:,1)] = pdegrad(xy,ele(2:4,:),V(:,i));
      Ei = pdeprtni(xy,ele(2:4,:),sqrt(sum(Fxy.^2,2)).');
      Fxy(:,2) = -Fxy(:,2);
      pdeplot( xy,[],ele(2:4,:),'xydata',Ei, ...
       'flowdata',Fxy,'flowstyle','arrow','colorbar','off');     
      title(['TE ',sprintf('@ %02.5g GHz\n',f(i)*1e-9)])
    end
    axis equal;
    axis(1e-3*[-a*0.6 a*0.6 -b*0.6 b*0.6]);
%     axis off;
    itp = itp - 1;
    if ~itp
      close all;
      break
    end
    fprintf('freq: %g GHz\n',f(i)*1e-9)
    pause;
    clf(1);
  end
end
%% eigenfrequencies
m = [1 2 0 1 1 2 2 3 3 3 4 0 1 1 4 4 2 2 5 3 3 5 5 4 4 6 6 0 1 1];
n = [0 0 1 1 1 1 1 0 1 1 0 2 2 2 1 1 2 2 0 2 2 1 1 2 2 0 1 3 3 3]; 
k = 1e3*sqrt((m*pi/a).^2+(n*pi/b).^2);
k = sort(k).';
fc = c0/2/pi*k;
fc = fc*1e-9;
fcn = f.'*1e-9;

fprintf('Relative error: %2.4g\n',(kFem(1:4)-k(1:4))./k(1:4));
return

if geom == 1
%   fprintf('Relative error fcTE10: %2.4g\n',(fcn(1)-fc(1))/fc(1));
  figure(1); plot(0:20,fcn(1:21),'x',0:20,fc(1:21),'o'); 
  legend('FEM','Analytical','location','best');
  figure(2); semilogy(0:20,abs((fcn(1:21)-fc(1:21))./fc(1:21)),'x');
  legend('Err','location','best');
else
  figure;
  plot(0:20,fcn(1:21),'o');
  legend('FEM','location','best');
end
ylabel('f_c[GHz]')
xlabel('mode [-]')
fprintf('Cut-off 1st higher mode: %2.4g fc1\n',fcn(2)/fcn(1));
fprintf('Monomodal bandwidth: %2.4g fc1\n',(fcn(2)-fcn(1))/fcn(1));
%   plot(1:30,(fc(1:30)-fcn(1:30))./fc(1:30),'.-')