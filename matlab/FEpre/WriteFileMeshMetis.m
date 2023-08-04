function WriteFileMeshMetis(ele,k)
%% Scrive il file .mesh che utilizza Metis per partizionare
% ele = matrice elementi --> numeroElementi x (3 o 4)
% k = 1 se l'elemento è un triangolo, 2 se tetraedro
meshID = fopen('FileMeshMetis.mesh','w');

%% Prima riga: numero elementi; numero di sottodomini
fprintf(meshID, '%d %d \n', length(ele), k);

%% Le restanti righe sono i nodi degli elementi
for i=1:length(ele)
fprintf(meshID, '%d %d %d \n', ele(i,1), ele(i,2), ele(i,3));
end

fclose(meshID);

%WriteFileMeshMetis(Mesh.ele,1)
%system(['partdmesh.exe ', 'FileMeshMetis.mesh ', num2str(regioni)])
%Mesh.elab = importdata(['FileMeshMetis.mesh.epart.','',num2str(regioni)])+1;
%Mesh.nlab = importdata(['FileMeshMetis.mesh.npart.','',num2str(regioni)])+1;

%% Label spigoli d'interfaccia tra i sottodomini
%auxvect=(1:Mesh.NELE)';
%eleR1 = auxvect(Mesh.elab(:)==1); % indice degli elementi regione 1
%eleR2 = auxvect(Mesh.elab(:)==2); % indice degli elementi regione 2
%spigR1 = unique(sort(abs(Mesh.spig(eleR1,:))));
%spigR2 = unique(sort(abs(Mesh.spig(eleR2,:))));
%spigR1R2=intersect(spigR1,spigR2);
%Mesh.slab(spigR1R2)=12;