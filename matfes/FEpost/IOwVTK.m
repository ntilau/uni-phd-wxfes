function IOwVTK(Sys, Mesh, filename)
% write vtk file for mesh plot
tic;
fid = fopen(['FEpost/',filename, '.vtk'], 'w');
fprintf(fid, '# vtk DataFile Version 2.0\nUnstructured Grid\n');
fprintf(fid, 'ASCII\nDATASET UNSTRUCTURED_GRID\n');
fprintf(fid, 'POINTS %d double\n', length(Mesh.refNode));
for i=1:length(Mesh.refNode)
  fprintf(fid, '%g %g 0\n', Mesh.refNode(i,1), Mesh.refNode(i,2));
end
fprintf(fid, 'CELLS %d %d\n', length(Mesh.refEle), 4*length(Mesh.refEle));
for i=1:length(Mesh.refEle)
  fprintf(fid, '3 %d %d %d\n', Mesh.refEle(i,1)-1, ...
    Mesh.refEle(i,2)-1, Mesh.refEle(i,3)-1);
end
fprintf(fid, 'CELL_TYPES %d\n', length(Mesh.refEle));
for i=1:length(Mesh.refEle)
  fprintf(fid, '9\n');
end
% fprintf(fid, 'CELL_DATA %d\n', 1);
fprintf(fid, 'POINT_DATA %d\n', length(Mesh.refNode));
% fprintf(fid, 'CELL_DATA %d\n', 1);
% fprintf(fid, 'SCALARS real float 1\n');
% fprintf(fid, 'LOOKUP_TABLE my_table\n');
% for i=1:length(Mesh.refNode)
%   fprintf(fid, '%f\n', real(data((i-1)*1+1)) );
% end
fprintf(fid, 'SCALARS Eabs float 1\n');
fprintf(fid, 'LOOKUP_TABLE jet\n');
for i=1:length(Mesh.refNode)
    fprintf(fid, '%f\n', abs(Sys.u(i)) );
end
fprintf(fid, 'SCALARS Ereal float 1\n');
fprintf(fid, 'LOOKUP_TABLE jet\n');
for i=1:length(Mesh.refNode)
    fprintf(fid, '%f\n', real(Sys.u(i)) );
end
fprintf(fid, 'SCALARS Eimag float 1\n');
fprintf(fid, 'LOOKUP_TABLE jet\n');
for i=1:length(Mesh.refNode)
    fprintf(fid, '%f\n', imag(Sys.u(i)) );
end
% for i=1:nHarm
%   fprintf(fid, ['SCALARS field', num2str(i), ' float 1\n']);
%   fprintf(fid, 'LOOKUP_TABLE my_table\n');
%   for j=1:length(Mesh.refNode)
%    fprintf(fid, '%f\n', abs(data(nHarm*(j-1)+i)) );
%   end
% end

%fprintf(fid, 'VECTORS vectors float\n');

fclose(fid);
fprintf('Writing VTK: %2.4g s\n',toc);