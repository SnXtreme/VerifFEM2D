function [Error_RExt,tk] = ResiduExt(omega,sigma,domega,dtop,Fd)
% Calcul de l'erreur via la méthode des résidus explicite

isvid=false;

[conn,eta,normal] = interfaceElems(omega);
nInterf = size(conn,1);
tk = zeros(omega.nbElems,2);
loading_interf_ids = ismember(sort(conn,2),sort(dtop.elems,2),'rows');
border_interf_ids = ismember(sort(conn,2),sort(domega.elems,2),'rows');
elemSize = omega.elementSize;
SigmaMat = @(id_sigma,id) [sigma(id_sigma(1,id)),sigma(id_sigma(3,id));...
                      sigma(id_sigma(3,id)),sigma(id_sigma(2,id))];

% Pour la réalisation d'une vidéo
if isvid;
v=VideoWriter('ResiduExt.avi');
open(v);
f=figure;
u = uicontrol('Style','slider','Position',[10 50 20 340],...
    'Min',1,'Max',nInterf,'Value',1);
end;

% Calcul du saut sigma à chaque interface. La normale est extérieur au
% premier élément
for id_interf=1:nInterf
    elems_ids = find(eta(id_interf,:));
        id_sigma = bsxfun(@(id,j) 3*(id-1)+j,elems_ids,(1:3)');
    if ~border_interf_ids(id_interf)
        % Interface interne
        saut_T = (SigmaMat(id_sigma,1)*eta(id_interf,elems_ids(1))...
            +SigmaMat(id_sigma,2)*eta(id_interf,elems_ids(2)))*normal(id_interf,:)';
        tk(elems_ids,:) = tk(elems_ids,:) + (saut_T*eta(id_interf,elems_ids))';
    elseif loading_interf_ids(id_interf)
        % Bord du problème chargé
        saut_T = ((SigmaMat(id_sigma,1)*normal(id_interf,:)')'-Fd)';
        tk(elems_ids,:) = tk(elems_ids,:) + (saut_T*eta(id_interf,elems_ids))';
    end
if isvid;
    coor_Gauss=CallGaussPoints(omega);
    oelems=omega.restrict(@(x) ismember(x,omega.nodes(omega.elems(elems_ids,:),:),'rows'));
    hold on
    plotElemField(oelems,sigma(id_sigma));
    scatter(coor_Gauss(elems_ids,1),coor_Gauss(elems_ids,2),'filled');
    quiver(sum(omega.nodes(conn(id_interf,:),1))/2,sum(omega.nodes(conn(id_interf,:),2))/2,...
        saut_T(1)/100,saut_T(2)/100);
    hold off
    u.Value = id_interf;
    M(id_interf) = getframe(gcf);
    writeVideo(v,getframe(gcf));
end
end
    % Integration numerique de Gauss
    [Wg,~] = gaussPoints(omega);
    tk_norm = Wg*arrayfun(@(id) tk(id,:)*tk(id,:)',(1:size(tk,1))');
    tk=[tk,zeros(omega.nbElems,1)];
    tk=reshape(tk',[],1);
    Error_RExt = sqrt(1/24)*sqrt(elemSize.*tk_norm);

if isvid
close(v);
end
end