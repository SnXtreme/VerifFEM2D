function ChElem = pt_Gauss(omega,ChNode)
% Définition de la valeur du champ par point aux points de Gauss
%
% Parametres:
%   - omega : le maillage support du champs de deplacement
%   - ChNode : le champ par point à modifier
%   - nbcomp : le nombre de composantes du ChNode
% Retourne le champs de contrainte (forme vectorielle)
% [\sigma_11,\sigma_22, \sigma12] en chaque points de Gauss de chaques
% elements
%
% ! : Le code ne fonctionne pour le moment qu'avec des TRI3 (Ordre 1)
    
    B = eye(3);
    nbcomp = 3;
    ncomp = length(B);
    % Integration numerique de Gauss
    [Wg,Xg] = gaussPoints(omega);
    
    Ng = numel(Wg);
    ChElem = zeros(ncomp*omega.nbElems*Ng,1);
    for i=1:omega.nbElems
        ids = omega.elems(i,:); % ids des noeuds
        mapu = bsxfun(@(id,j) (id-1)*nbcomp+j,ids(:)',(1:nbcomp)'); % index de l'inconnue u
        maps = (i-1)*Ng*ncomp + (1:Ng*ncomp); % index de l'inconnue sigma
        Xe = omega.nodes(ids,:); % coordonnees de l'element

        M1 = shapesFunctions(omega,Xg,Xe,0,3);

        D = kron(eye(Ng),B);

        ChElem(maps(:),1) = ChElem(maps(:),1) + M1*ChNode(mapu(:));
    end

end