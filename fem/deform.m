function new_mesh = deform(mesh,u,coeff)
% Cr�er un maillage d�form�e � partir d'un champs de d�placement 2D
%
% Param�tres:
%   - le maillage initial
%   - un champs de d�placement 2D
%   - un coefficient d'amplification (optionel)
    
    % V�rification des entr�es
    assert(isa(mesh,'Mesh'),'Objet maillage inconnu');
    assert(isnumeric(u) && numel(u) == 2*mesh.nbNodes,'Mauvais champs de d�placement 2D');

    if nargin == 2
        coeff = 1;
    else
        assert(isnumeric(coeff),'Mauvais facteur d''amplification');
        coeff = coeff(1);
    end
    
    if any(size(u) == 1)
        u = reshape(u(:),[],mesh.nbNodes)';
    end
    
    % Cr�ation du nouveau maillage
    new_mesh = Mesh(mesh.nodes+coeff*u,mesh.elems);
end