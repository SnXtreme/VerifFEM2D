function coor_Gauss=CallGaussPoints(omega)

% Calcul des points de Gauss dans le repère global
coor_Gauss=zeros(size(omega.elems,1),omega.dim);
[~,Xg]=gaussPoints(omega);
N=[Xg(:,1)  Xg(:,2)  1-Xg(:,1)-Xg(:,2)];
for id_elem=1:size(omega.elems,1)
%     Xe = omega.nodes(omega.elems(id_elem,:));
%     [~,Xg]=gaussPoints(omega);
%     N=shapesFunctions(omega,Xg,Xe,0);
    coor_Gauss(id_elem,:)=N*omega.nodes(omega.elems(id_elem,:),:);
end

end