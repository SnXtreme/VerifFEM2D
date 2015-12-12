function sigmal = stressl(omega,sigma)
% Calcule le champs de contrainte lissé a partir d'un champs de deplacement 2D
%
% Parametres:
%   - omega : le maillage support du champs de deplacement
%   - sigma : Terseur de contrainte exprimé aux points de Gauss
%
% Retourne le champs de contrainte (forme vectorielle)
% [\sigma_11,\sigma_22, \sigma12] en chaque noeud

sdeco=reshape(sigma,3,[])';

sigmal = zeros(3*length(omega.nodes),1);
for nnode=1:length(omega.nodes)
    Masq=any(omega.elems==nnode,2);
    sigmal(3*(nnode-1)+1)=(transpose(sdeco(:,1)) * Masq)/(transpose(Masq)*double(Masq));
    sigmal(3*(nnode-1)+2)=(transpose(sdeco(:,2)) * Masq)/(transpose(Masq)*double(Masq));
    sigmal(3*(nnode-1)+3)=(transpose(sdeco(:,3)) * Masq)/(transpose(Masq)*double(Masq));
end

end