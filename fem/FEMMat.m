function [A,M] = FEMMat(omega,B)
% Construit la matrice A elements finis du systeme Au=b, tel que :
%   A = \int_{omega} [\epsilon(N)]^T*B*[\epsilon(N)] dx
%   et optionnellement
%   M = \int_{omega} N^T*eye(2)*N dx
%
%   B est le tenseur de Hook, ecrit sous forme matriciel avec les notations
%   de Voigt.

    % Verification des donnees
    assert(isa(omega,'Mesh'),'Parametre #1 invalide : objet maillage type invalide');
    assert(isnumeric(B) && (all(size(B) == [3 3]) || all(size(B) == [2 2])),'Mauvaise representation du tenseur de Hook');

    A = matrixAssembly(omega,1,B);
    if nargout == 2
        M = matrixAssembly(omega,0,eye(2));
    end
end
