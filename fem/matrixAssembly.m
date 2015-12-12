function A = matrixAssembly(omega,order,B,dim)
    if nargin == 3;
        dim = 2;
    end;

    % Determinant de la matrice jacobienne
    detJ = @(J) J(1:2:end,1).*J(2:2:end,2) - J(2:2:end,1).*J(1:2:end,2);

    % Integration numerique de Gauss
    [Wg,Xg] = gaussPoints(omega,2*(omega.order-order));

    % Creation de la matrice
    A = sparse(dim*omega.nbNodes,dim*omega.nbNodes);
    for i=1:omega.nbElems
        ids = omega.elems(i,:); % ids des noeuds
        map = bsxfun(@(id,j) (id-1)*dim+j,ids(:)',(1:dim)'); % index de l'inconnue
        Xe = omega.nodes(ids,:); % coordonnees de l'element

        [M1,J] = shapesFunctions(omega,Xg,Xe,order,dim); % Evaluation des fonctions de formes

        D = kron(diag(Wg.*detJ(J)),B); % Matrice pour la quadrature de Gauss

        A(map(:),map(:)) = A(map(:),map(:))+(M1'*D*M1);
    end
end