function [uF, Sys] = SolvHBDDschur(Sys)
%     tic
    SF = Sys.AFF*0;
    for ir=1:length(Sys.RegDoF)
        if Sys.Compute(ir)
%             fprintf('reg %d computed\n',ir);
            Sys.SFprev{ir} = - Sys.AFI{ir}*(Sys.AII{ir}\Sys.AIF{ir});
        end
            
%         if norm(full(SFtmp-Sys.SFprev{ir}))==0
%             fprintf('reg %d did not change\n',ir);
%         end
        SF = SF + Sys.SFprev{ir};
%         figure;spy(SF)
    end
    SF = SF+Sys.AFF;
%     figure;spy(SF)
    uF = SF\Sys.gF;
%     toc
end