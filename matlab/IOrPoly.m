function Mesh = IOrPoly(filename, args, hOrd, scal)
% Mesh = IOrPoly(filename, args, hOrd, scal)
% reads Triangle mesh generated files
% filename = name of the *.poly file without extension
% args =  selected arguments
% 
% The following args are enforced:
% -p  Triangulates a Planar Straight Line Graph (.poly file).
% -D  Conforming Delaunay:  all triangles are truly Delaunay.
% -e  Generates an edge list.
% 
% remains to trim with:
% -q  Quality mesh generation.  A minimum angle may be specified.
% -a  Applies a maximum triangle area constraint.
% -u  Applies a user-defined triangle constraint.
% -A  Applies attributes to identify triangles in certain regions.
% -o2 Generates second-order subparametric elements.
tic
system(['cd FEpre/ && ./IOrMesh.exe ./', filename, ' ', args]);
load(['FEpre/',filename,'.h', num2str(hOrd),'.mat']);
Mesh.node = node*scal;
Mesh.ele = ele;
Mesh.spig = spig;
Mesh.spig2 = spig2;
Mesh.nlab = nlab;
Mesh.elab = elab;
Mesh.slab = slab;
Mesh.NNODE = size(node,1);
Mesh.NELE = size(ele,1);
Mesh.NSPIG = size(spig2,1);
fprintf('Meshing geometry: %2.4g s\n',toc);
return

system(['triangle -p', args , 'De ',  filename, '.poly']);

% nodes
if exist([filename,'.1.node'],'file')
  tmp = importdata([filename,'.1.node'], ' ', 1);
  tmpInfo = getTmpInfo(tmp);
  node = tmp.data(:,2:2+tmpInfo(2)-1);
  % idx = 2+tmpInfo(2);
  % natt = tmp.data(:,idx:idx+tmpInfo(3)-1);
  idx = 2+tmpInfo(2)+tmpInfo(3);
  nlab = tmp.data(:,idx);
else
  node = 0;
  nlab = 0;
end
% elements
if exist([filename,'.1.ele'],'file')
  tmp = importdata([filename,'.1.ele'], ' ', 1);
  tmpInfo = getTmpInfo(tmp);
  ele = sort(tmp.data(:,2:2+tmpInfo(2)-1),2);
  idx = 2+tmpInfo(2);
  elab = tmp.data(:,idx:idx+tmpInfo(3)-1);
%   nNodesEle = tmpInfo(2);
else
  ele = 0;
  elab = 0;
end
% edges
if exist([filename,'.1.edge'],'file')
  tmp = importdata([filename,'.1.edge'], ' ', 1);
  tmpInfo = getTmpInfo(tmp);
  spig2 = sort(tmp.data(:,2:3),2);
  slab = tmp.data(:,4:4+tmpInfo(2)-1);
else
  Mesh.spig2 = 0;
  Mesh.slab = 0;
end

% Mesh.ele(1+(1:nNodesEle),:) = Mesh.ele;
% Mesh.ele(1,:) = nNodesEle;
% Mesh.spig = ele2spig(Mesh.ele, Mesh.spig2);

NNODE = size(node,1);
NELE = size(ele,1);
% nlab = zeros(NNODE,1);
% nlab(node(:,1) == 0 | node(:,1) == 1 | node(:,2) == 0 | node(:,2) == 1) = 1;
% elab = ones(NELE,1);
% tspig = [];
% % slab = [];
% for ie=1:NELE
%   ele(ie,:) = sort(ele(ie,:));
%   nodes = ele(ie,:);
%   for i=1:3
%     snodes = sort([nodes(i) nodes(mod(i,3)+1)]);
%     if isempty(tspig)
%       tspig = [tspig snodes.'];
% %       if (nlab(snodes(1)) == 1) && (nlab(snodes(2)) == 1) && sum(diff(node(snodes,:))==0)
% %           slab = [slab; 1];
% %       else
% %           slab = [slab; 0];
% %       end
%     else
%       if find(tspig(1,:) == snodes(1) & tspig(2,:) == snodes(2))
%         % disp('exists1')
%       elseif find(tspig(1,:) == snodes(2) & tspig(2,:) == snodes(1))
%         % disp('exists2')
%       else
%         tspig = [tspig snodes.'];
% %         if (nlab(snodes(1)) == 1) && (nlab(snodes(2)) == 1) && sum(diff(node(snodes,:))==0)
% %             slab = [slab; 1];
% %         else
% %             slab = [slab; 0];
% %         end
%       end
%     end
%   end  
% end
% spig2 = tspig.';
tspig = spig2.';
NSPIG = length(spig2);
% slab = zeros(1,NSPIG);
spig = zeros(NELE,3);
for ie=1:NELE
    ele(ie,:) = sort(ele(ie,:));
   nodes = ele(ie,:);
    for i=1:3
      snodes = [nodes(mod(i,3)+1) nodes(mod(i+1,3)+1)];
      %snodes = [nodes(i) nodes(mod(i,3)+1)];
      if isempty(find(tspig(1,:) == snodes(1) & tspig(2,:) == snodes(2), 1))
        spig(ie,i) = - find(tspig(1,:) == snodes(2) & tspig(2,:) == snodes(1));
      else
        spig(ie,i) = find(tspig(1,:) == snodes(1) & tspig(2,:) == snodes(2));
      end
    end
end
% trimesh(ele,node(:,1),node(:,2))
% hold on;
% plot(node(nlab==1,1),node(nlab==1,2),'o')
Mesh.node = node;
Mesh.ele = ele;
Mesh.spig = spig;
Mesh.spig2 = spig2;
Mesh.nlab = nlab;
Mesh.elab = elab;
Mesh.slab = slab;
Mesh.NNODE = NNODE;
Mesh.NELE = NELE;
Mesh.NSPIG = NSPIG;

system('del *.node *.ele *.edge *.1.poly');

end


function tmpInfo = getTmpInfo(tmp)
  nbr = length(tmp.textdata);
  if nbr == 1
    tmpInfo = str2num(tmp.textdata{1}); %#ok<*ST2NM>
  else
    tmpInfo = zeros(1,nbr);
    for i=1:nbr
      tmpInfo(i) =  str2num(tmp.textdata{i});
    end
  end
end

% function e2s = ele2spig(ele, tspig)
% 
% e2s = zeros( 4, length(ele));
% for ie=1:length(ele)
%   nodes = ele(2:4,ie);
%   for i=1:3
%     snodes = [nodes(i) nodes(mod(i,3)+1)];
%     if isempty(find(tspig(1,:) == snodes(1) & tspig(2,:) == snodes(2), 1))
%       e2s(i+1,ie) = - find(tspig(1,:) == snodes(2) & tspig(2,:) == snodes(1));
%     else
%       e2s(i+1,ie) = find(tspig(1,:) == snodes(1) & tspig(2,:) == snodes(2));
%     end
%   end
% end
% e2s(1,:) = 3;
% end