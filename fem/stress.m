function A = stress(omega,B,u)
% Calcule le champs de contrainte � partir d'un champs de d�placement 2D
%
% Param�tres:
%   - omega : le maillage support du champs de d�placement
%   - B : le tenseur de Hook dans sa notation de Voigt (matrice)
%   - u : le champs de d�placement 2D
%
% Retourne le champs de contrainte (forme vectorielle)
% [\sigma_11,\sigma_22, \sigma12] en chaque points de Gauss de chaques
% �l�ments

    % V�rification des entr�es
    assert(isa(omega,'Mesh'),'Objet maillage inconnu');
    assert(isnumeric(B) && all(size(B) == [3 3]),'Mauvaise repr�sentation du tenseur de Hook');
    assert(isnumeric(u) && numel(u) == 2*omega.nbNodes,'Mauvais champs de d�placement 2D');

    % Int�gration num�rique de Gauss
    [Wg,Xg] = gaussPoints(omega);
    
    Ng = numel(Wg);
    A = zeros(3*omega.nbElems*Ng,1);
    for i=1:omega.nbElems
        ids = omega.elems(i,:); % ids des noeuds
        mapu = bsxfun(@(id,j) (id-1)*2+j,ids(:)',(1:2)'); % index de l'inconnue u
        maps = (i-1)*Ng*3 + (1:Ng*3); % index de l'inconnue sigma
        Xe = omega.nodes(ids,:); % coordonn�es de l'�l�ment
        
        M1 = shapesFunctions(omega,Xg,Xe,1);
        
        D = kron(eye(Ng),B);

        A(maps(:),1) = A(maps(:),1) + D*M1*u(mapu(:));
    end
    
end