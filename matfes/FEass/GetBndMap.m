function [Mesh, Sys] = GetBndMap(Sys, Mesh)
idsDD = find(Mesh.slab == Mesh.BC.DDschur);
RegEle = cell(length(unique(Mesh.elab)),1);
RegDoF = cell(length(unique(Mesh.elab)),1);
RegDoFmap = cell(length(unique(Mesh.elab)),1);
BndDoF = [];
for ie=1:Mesh.NELE
    gIs = CalcGlobIndex(2, Sys.pOrd, Mesh, ie);
    tmp = Mesh.elab(ie);
    RegEle{tmp} = unique([RegEle{tmp}; ie]);
    RegDoF{tmp} = unique([RegDoF{tmp}; gIs.']);
    onDD = [sum(abs(Mesh.spig(ie,1)) == idsDD) ...
        sum(abs(Mesh.spig(ie,2)) == idsDD) ...
        sum(abs(Mesh.spig(ie,3)) == idsDD)];
    if sum(onDD) > 0
        idOnDD = find(onDD==1);
        for i=1:length(idOnDD)
            spigId = Mesh.spig(ie,idOnDD(i));
            nodeId = Mesh.spig2(abs(spigId),:);
            gIs = CalcGlobIndex(1, Sys.pOrd, Mesh, ie, idOnDD(i));
            BndDoF = [BndDoF gIs];
        end
    end
end
BndDoF = unique(BndDoF);
Sys = CalcDoFsNumber(Sys, Mesh);
DDmap = zeros(Sys.NDOFs,1);
roof = length(BndDoF); 
% roof = 0;
% regOrder = [Mesh.NLlab setdiff(1:length(RegDoF), Mesh.NLlab)];
regOrder = 1:length(RegDoF);
for ir=regOrder
    RegDoF{ir} = setdiff(RegDoF{ir},BndDoF);
    RegDoFmap{ir} = roof+(1:length(RegDoF{ir})).';
    DDmap(RegDoF{ir}) = RegDoFmap{ir};
    roof = roof+length(RegDoF{ir});
end

DDmap(BndDoF) = 1:length(BndDoF); 
% DDmap(BndDoF) = roof+(1:length(BndDoF));
Sys.DDmap = DDmap;
Sys.BndDoF = BndDoF.';
Sys.RegEle = RegEle;
Sys.RegDoF = RegDoF;
Sys.RegDoFmap = RegDoFmap;