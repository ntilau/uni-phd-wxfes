clear all; clc; close all;
Config();
prjName = 'CircKoshiba26_8_plot';
% HarmsRetained = 8;
load(prjName);%, '_', num2str(HarmsRetained)]);
%% Field
for ih=1:4%Sys.HBharmPlot
    Sys.u = Sys.uPlot((ih-1)*Sys.NDOFs+(1:Sys.NDOFs));
    if exist('pdeplot','file')
        figure; pdeplot(Mesh.refNode.',[],Mesh.refEle.',...
            'xydata',(abs(Sys.u)),...
            'mesh','off','contour','on','levels',10,...
            'colormap','jet','xygrid','off');
        axis equal; axis tight; camlight left; lighting phong;
        title([num2str(Sys.HBharms(ih)*Sys.freq*1e-9), 'GHz'])
    else
        IOwVTK(Sys, Mesh, prjName);
    end
end