function [ZZ_Value]=ZZ(omega,B,sigma,sigmal)

% Calcul de l'estimateur ZZ1 "en moyenne"
% Parametres:
%   - omega : le maillage support du champs de deplacement
%   - B : le tenseur de Hook dans sa notation de Voigt (matrice)
%   - sigma : Tenseur de contrainte exprimé aux points de Gauss
%   - sigmal : Tenseur de contrainte lissé aux noeuds
% Sortie de la valeur de ZZ1 et d'un calcul d'écart local

    B = eye(3);
    nbcomp = 3;
    ncomp = length(B);
    % Integration numerique de Gauss
    [Wg,Xg] = gaussPoints(omega);
    
    Ng = numel(Wg);
    ZZ_Value = zeros(omega.nbElems,1);
    for i=1:omega.nbElems
        dS = zeros(ncomp*Ng,1);
        ids = omega.elems(i,:); % ids des noeuds
        mapu = bsxfun(@(id,j) (id-1)*nbcomp+j,ids(:)',(1:nbcomp)'); % index de sigma lissé
        maps = (i-1)*Ng*ncomp + (1:Ng*ncomp); % index de l'inconnue sigma
        Xe = omega.nodes(ids,:); % coordonnees de l'element
        M1 = shapesFunctions(omega,Xg,Xe,0,3);
        dS(:,1) = M1*sigmal(mapu(:)) - sigma(maps(:),1);
        ZZ_Value(i,1) = dS'*inv(B)*dS;
    end
    