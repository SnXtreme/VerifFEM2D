% Charge le maillage
%omega = Mesh(geo2msh('./meshes/plate_crack.geo',1,1/10));

% Extraction de maillage secondaires
% H = max(omega.nodes(:,2));
domega = omega.border;
% dtop = domega.restrict(@(x) x(:,2) == H);
% dleft = domega.restrict(@(x) x(:,1) == 0);
% dbottom = domega.restrict(@(x) x(:,2) == 0 & x(:,1) >= 1);

% Détermination des interfaces
[conn,eta,normal] = interfaceElems(omega);

figure('Name','Problème conn')
hold on
loading_interf_ids = ismember(sort(conn,2),sort(dtop.elems,2),'rows');
border_interf_ids = ismember(sort(conn,2),sort(domega.elems,2),'rows');
for id_interf=(find(border_interf_ids))'
	elems_ids = find(eta(id_interf,:));
    nodes_elems = unique(omega.elems(elems_ids,:));
    plot(omega.restrict(@(x) ismember(x,omega.nodes(omega.elems(elems_ids,:),:),'rows')))
    scatter(omega.nodes(nodes_elems,1),omega.nodes(nodes_elems,2),'filled');
end