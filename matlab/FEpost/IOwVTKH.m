function IOwVTKH(Sys, Mesh, filename)
% write vtk file for mesh plot
tic;
fid = fopen(['FEpost/',filename, '.vtk'], 'w');
fprintf(fid, '# vtk DataFile Version 2.0\nUnstructured Grid\n');
fprintf(fid, 'ASCII\nDATASET UNSTRUCTURED_GRID\n');
fprintf(fid, 'POINTS %d double\n', Mesh.NNODE);
for i=1:Mesh.NNODE
  fprintf(fid, '%g %g 0\n', Mesh.node(i,1), Mesh.node(i,2));
end
fprintf(fid, 'CELLS %d %d\n', Mesh.NELE, 4*Mesh.NELE);
for i=1:Mesh.NELE
  fprintf(fid, '3 %d %d %d\n', Mesh.ele(i,1)-1, ...
    Mesh.ele(i,2)-1, Mesh.ele(i,3)-1);
end
fprintf(fid, 'CELL_TYPES %d\n', Mesh.NELE);
for i=1:Mesh.NELE
  fprintf(fid, '9\n');
end
% fprintf(fid, 'CELL_DATA %d\n', 1);
% fprintf(fid, 'POINT_DATA %d\n', Mesh.NNODE);
fprintf(fid, 'CELL_DATA %d\n', Mesh.NELE);
% fprintf(fid, 'SCALARS real float 1\n');
% fprintf(fid, 'LOOKUP_TABLE my_table\n');
% for i=1:Mesh.NNODE
%   fprintf(fid, '%f\n', real(data((i-1)*1+1)) );
% end
fprintf(fid, 'SCALARS Hx float 1\n');
fprintf(fid, 'LOOKUP_TABLE jet\n');
for i=1:Mesh.NELE
    fprintf(fid, '%f\n', abs(Sys.u(i,1)) );
end
fprintf(fid, 'SCALARS Hy float 1\n');
fprintf(fid, 'LOOKUP_TABLE jet\n');
for i=1:Mesh.NELE
    fprintf(fid, '%f\n', abs(Sys.u(i,2)) );
end
% fprintf(fid, 'SCALARS Ereal float 1\n');
% fprintf(fid, 'LOOKUP_TABLE jet\n');
% for i=1:Mesh.NNODE
%     fprintf(fid, '%f\n', real(Sys.u(i)) );
% end
% fprintf(fid, 'SCALARS Eimag float 1\n');
% fprintf(fid, 'LOOKUP_TABLE jet\n');
% for i=1:Mesh.NNODE
%     fprintf(fid, '%f\n', imag(Sys.u(i)) );
% end
% for i=1:nHarm
%   fprintf(fid, ['SCALARS field', num2str(i), ' float 1\n']);
%   fprintf(fid, 'LOOKUP_TABLE my_table\n');
%   for j=1:Mesh.NNODE
%    fprintf(fid, '%f\n', abs(data(nHarm*(j-1)+i)) );
%   end
% end

%fprintf(fid, 'VECTORS vectors float\n');

fclose(fid);
fprintf('Writing VTK: %2.4g s\n',toc);