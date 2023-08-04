function [SysF, MeshF, errorF, spF] = AssembSolveFull(SysF,MeshF)

[SysF,MeshF] = AssembHBFerrite(SysF, MeshF);
SysF =  AssembWPHB(SysF);      
XF =  SolvDir(SysF);
SysF.u = zeros(SysF.NDOFs*SysF.nHarms,1);
SysF.u(SysF.nnWP) = XF(length(SysF.WP)*SysF.WPnModes*SysF.nHarms+1:end, 1);
for jh=1:length(SysF.Harms)
    for ip=1:length(SysF.WP)
        SysF.u((jh-1)*SysF.NDOFs+SysF.WP{ip}) = SysF.WPgvec{ip,jh}*...
            (XF((1:SysF.WPnModes)+((jh-1)+(ip-1)*SysF.nHarms)*SysF.WPnModes,1).*...
            SysF.Power((1:SysF.WPnModes)+((jh-1)+(ip-1)*SysF.nHarms)*SysF.WPnModes).');
    end
end
errorF = norm(SysF.u - SysF.u0) / norm(SysF.u);
SysF.u0 = SysF.u;
spF = XF(1:length(SysF.WP)*SysF.WPnModes*SysF.nHarms,1);
spF(1,1) = (spF(1,1) - 1) * sqrt(SysF.Pfund);
spF(SysF.WPnum,1) = (spF(SysF.WPnum,1) - 1) * sqrt(SysF.Pitrf);