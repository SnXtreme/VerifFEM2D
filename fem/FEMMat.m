function A = FEMMat(omega,od1,od2,B)
% Construit la matrice A �l�ments finis du syst�me Au=b, tel que :
%   A = \int_{omega} [d^{od1}N/dx]^T*B*[d^{od2}N/dx] dx
%
%   od1 et od2 sont les ordres de d�rivations �gaux � 0 ou 1. B est le
%   tenseur de Hook, �crit sous forme matriciel avec les notations de
%   Voigt.

    % V�rification des donn�es
    assert(isa(omega,'Mesh'),'Param�tre #1 invalide : objet maillage type invalide');
    assert((od1==0 || od1==1) && (od2==0 || od2==1),'L''ordre de d�rivation doit �tre 0 ou 1');
    assert(isnumeric(B) && all(size(B) == [3 3]),'Mauvaise repr�sentation du tenseur de Hook');
    
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
        
        [M1,J] = shapesFunctions(omega,Xg,Xe,od1); % Evaluation des fonctions de formes
        M2 = shapesFunctions(omega,Xg,Xe,od2); % Evaluation des fonctions de formes
        
        D = kron(diag(Wg.*detJ(J)),B); % Matrice pour la quadrature de Gauss
        
        A(map(:),map(:)) = A(map(:),map(:))+(M1'*D*M2);
    end
end