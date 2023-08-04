function [Sys, Mesh, errorD, spD] = AssembSolveDD(Sys,Mesh)

if ~isfield(Sys,'first')
    Sys.RegToCompute = 1:Mesh.NLlab-1;
    [Sys,Mesh] = AssembHBDDFerrite(Sys, Mesh);
    Sys.Slin = Sys.S;
    Sys.Tlin = Sys.T;
    Sys.first = true;
end
Sys.RegToCompute = Mesh.NLlab;
[Sys,Mesh] = AssembHBDDFerrite(Sys, Mesh);
Sys.S = Sys.S + Sys.Slin;
Sys.T = Sys.T + Sys.Tlin;
Sys =  AssembWPHBDDschur(Sys);
[XD, Sys] =  SolvHBDDschur(Sys);
Sys.Compute =  Sys.Compute*0;
Sys.Compute(Mesh.NLlab) = 1;
Sys.u = zeros(Sys.NDOFs*Sys.nHarms,1);
Sys.u(Sys.nnWP(Sys.BndDoFRed(length(Sys.WP)*Sys.WPnModes*Sys.nHarms+1:end)-length(Sys.WP)*Sys.WPnModes*Sys.nHarms)) = XD(length(Sys.WP)*Sys.WPnModes*Sys.nHarms+1:end,1);
for jh=1:length(Sys.Harms)
    for ip=1:length(Sys.WP)
        Sys.u((jh-1)*Sys.NDOFs+Sys.WP{ip}) = Sys.WPgvec{ip,jh}*...
            XD((1:Sys.WPnModes)+((jh-1)+(ip-1)*Sys.nHarms)*Sys.WPnModes,1);
    end
end
X5 = Sys.AII{Mesh.NLlab}\(-Sys.AIF{Mesh.NLlab}*XD(:,1));
Sys.u(Sys.nnWP(Sys.RegDoFRed{Mesh.NLlab}-length(Sys.WP)*Sys.WPnModes*Sys.nHarms)) = X5(:,1);
errorD = norm(Sys.u - Sys.u0) / norm(Sys.u);
Sys.u0 = Sys.u;
spD = XD(1:length(Sys.WP)*Sys.WPnModes*Sys.nHarms,1);
spD(1,1) = (spD(1,1) - 1) * sqrt(Sys.Pfund);
spD(Sys.WPnum,1) = (spD(Sys.WPnum,1) - 1) * sqrt(Sys.Pitrf);