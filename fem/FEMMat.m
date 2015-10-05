function [A,M] = FEMMat(omega,B)
% Construit la matrice A �l�ments finis du syst�me Au=b, tel que :
%   A = \int_{omega} [\epsilon(N)]^T*B*[\epsilon(N)] dx
%   et optionnellement
%   M = \int_{omega} N^T*eye(2)*N dx
%
%   B est le tenseur de Hook, �crit sous forme matriciel avec les notations
%   de Voigt.

    % V�rification des donn�es
    assert(isa(omega,'Mesh'),'Param�tre #1 invalide : objet maillage type invalide');
    assert(isnumeric(B) && (all(size(B) == [3 3]) || all(size(B) == [2 2])),'Mauvaise repr�sentation du tenseur de Hook');
    
    A = matrixAssembly(omega,1,B);
    if nargout == 2
        M = matrixAssembly(omega,0,eye(2));
    end
end

function A = matrixAssembly(omega,order,B)
    % D�terminant de la matrice jacobienne
    detJ = @(J) J(1:2:end,1).*J(2:2:end,2) - J(2:2:end,1).*J(1:2:end,2);
    
    % Int�gration num�rique de Gauss
    [Wg,Xg] = gaussPoints(omega);
    
    % Cr�ation de la matrice
    A = sparse(2*omega.nbNodes,2*omega.nbNodes);
    for i=1:omega.nbElems
        ids = omega.elems(i,:); % ids des noeuds
        map = bsxfun(@(id,j) (id-1)*2+j,ids(:)',(1:2)'); % index de l'inconnue
        Xe = omega.nodes(ids,:); % coordonn�es de l'�l�ment
        
        [M1,J] = shapesFunctions(omega,Xg,Xe,order); % Evaluation des fonctions de formes
        
        D = kron(diag(Wg.*detJ(J)),B); % Matrice pour la quadrature de Gauss
        
        A(map(:),map(:)) = A(map(:),map(:))+(M1'*D*M1);
    end
end