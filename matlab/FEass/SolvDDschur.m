function uF = SolvDDschur(Sys)
    SF = Sys.AFF*0;
    for ir=1:length(Sys.RegDoF)
        fprintf('%g\n',length(Sys.AII{ir}));
        tic
        SF = SF - Sys.AFI{ir}*(Sys.AII{ir}\Sys.AIF{ir});
        toc
    end
    tic
    SF = SF+Sys.AFF;
    uF = SF\Sys.gF;
    toc
end