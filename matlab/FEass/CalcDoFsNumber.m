function Sys = CalcDoFsNumber(Sys, Mesh)

Sys.NDOFs = Mesh.NNODE + Mesh.NSPIG*(Sys.pOrd-1) + ...
    Mesh.NELE*(0<(Sys.pOrd-2))*(Sys.pOrd-2)*(Sys.pOrd-1)*.5;

Sys.NDOFv = Sys.pOrd*Mesh.NSPIG*(1 + (Sys.pOrd>2)) + ...
    2*Mesh.NELE*(Sys.pOrd>1)*(Sys.pOrd-1);

end