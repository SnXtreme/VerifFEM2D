function sigmal = stress_ZZ2(omega,sigma,domega)
% Calcul le champ lissé via la méthode ZZ2
% Utilisation d'un nouveau repère (Noeud du patch,tilde_X,tilde_Y) pour
% améliorer le conditionnement de A

debug=false;
sigmal = zeros(omega.nbNodes,3);
% Calcul des points de Gauss de tous les éléments dans le repère global
coor_Gauss=CallGaussPoints(omega);

tilde_X=@(x,x_min,x_max) -1+((2*(x(:,1)-x_min(1)))./(x_max(1)-x_min(1)));
tilde_Y=@(x,x_min,x_max) -1+((2*(x(:,2)-x_min(2)))./(x_max(2)-x_min(2)));

P = @(x,x_min,x_max) [ones(size(x,1),1),...
    tilde_X(x,x_min,x_max),...
    tilde_Y(x,x_min,x_max),...
    tilde_X(x,x_min,x_max)*tilde_Y(x,x_min,x_max)];

border_ids = ismember(1:omega.nbNodes,sort(reshape(unique(domega.elems),[],1),2));

%% Cas noeud interne :
for id_node=find(~border_ids)
    [patch_elems_ids,~] = find(omega.elems == id_node);
    patch_id_nodes=unique(reshape(omega.elems(patch_elems_ids,:),[],1));
    
    % Coordonnées des noeuds du patch dans la base globale
    patch_coor_nodes=omega.nodes(patch_id_nodes,:);
    x_min=min(patch_coor_nodes(:,[1 2]));
    x_max=max(patch_coor_nodes(:,[1 2]));
    nb_comp_as=length(P([0,0],x_min,x_max));
    A = zeros(nb_comp_as);
    b = zeros(nb_comp_as,3);
    
    % Coordonnées de Gauss dans la base (noeud_patch,X,Y)
    patch_coor_Gauss=coor_Gauss(patch_elems_ids,:)-repmat(omega.nodes(id_node,:),length(patch_elems_ids),1);
    
    for id_patch_elem=1:length(patch_elems_ids)
        id_elem=patch_elems_ids(id_patch_elem);
        P_elem=P(patch_coor_Gauss(id_patch_elem,:),x_min,x_max);
        A = A + P_elem'*P_elem;
        for id_comp=1:3
            b(:,id_comp) = b(:,id_comp) + P_elem'*sigma(3*(id_elem-1)+id_comp);
        end
    end
    
    % Calcul des coefficients du polynôme et des valeurs sigma lissé du
    % noeud
    as=zeros(size(b,1),id_comp);
    for id_comp=1:3
        as(:,id_comp) = A\b(:,id_comp);
        sigmal(id_node,id_comp)=P([0,0],x_min,x_max)*as(:,id_comp);
    end
    
%% Spécial débug : Test du tracé du polynôme et comparaison avec les points connus
    if debug
        id_comp=1;
        % Grille exprimée dans le repère canonique
        [X,Y] = meshgrid(linspace(min(x_min),max(x_max),50));
        X=reshape(X,[],1);
        Y=reshape(Y,[],1);
        idbound=boundary(patch_coor_nodes(:,1),patch_coor_nodes(:,2));
        idinc=find(inpolygon(X,Y,patch_coor_nodes(idbound,1),patch_coor_nodes(idbound,2)));
        X=X(idinc);
        Y=Y(idinc);
        Z=arrayfun(@(id) P([X(id),Y(id)],x_min,x_max)*as(:,id_comp),(1:length(X))');
        hold on
        scatter3(X,Y,Z);
    for id_patch_elem=1:length(patch_elems_ids)
        id_elem=patch_elems_ids(id_patch_elem);
        % Coordonnée du point de Gauss dans le repère canonique
        cG=coor_Gauss(id_elem,:);
        scatter3(cG(1),cG(2),sigma(3*(id_elem-1)+id_comp),'filled');
    end
    hold off;
    end
end

%% Définition de la valeur des noeuds au bord comme la moyenne des éléments environnant
for id_node=find(border_ids)
    [patch_elems_ids,~] = find((omega.elems == id_node));
        for id_comp=1:3
        sigmal(id_node,id_comp)=(sum(sigma(3*([patch_elems_ids]-1)+id_comp)))/length(patch_elems_ids);
        end
end
sigmal=reshape(sigmal',[],1);
end