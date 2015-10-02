function pos_filename = adaptation(omega,element_size_map,pos_filename)
% Construit le fichier pos_filename n�cessaire pour faire de l'adaptation

% Param�tres obligatoires:
%   - le maillage omega sur lequel les calculs ont �t� fait; 
%   - un vecteur element_size_map contenant la taille optimale des �l�ments
%   pour chaques �l�ments du maillage omega, de taille le nombre d'elements
%
% Param�tre optionnel :
%   - pos_filename : permet de specifier le nom du fichier de sortie. Par
%   d�faut, le chemin est 'adaptation.pos'.
%
% Retourne le chemin du fichier d'adaptation.

    % Verification des param�tres
    assert(isa(omega,'Mesh'),'Param�tre #1 invalide : objet de type invalide');
    assert(isnumeric(element_size_map) && numel(element_size_map) == omega.nbElems,'Param�tre #2 invalide : mauvais type ou mauvaise taille');
    
    % Chemin du fichier de sortie
    if nargin < 3
        pos_filename = fullfile(pwd,'adaptation.pos');
    else
        assert(ischar(pos_filename),'Chemin invalide.');
    end

    % Transformation d'une donn�e �l�ment par �l�ment � une donn�e aux noeuds
    node_map = zeros(omega.nbNodes,1);
    for i=1:numel(node_map)
        elems_ids = any(omega.elems == i,2); % trouve le patch
        node_map(i) = mean(element_size_map(elems_ids)); % la moyenne des tailles sur ce patch
    end
    
    % Cr�er le fichier d'adaptation n�cessaire � GMSH
    fID = fopen(pos_filename,'w');
    fprintf(fID,'View "background mesh" {\n'); % Write the header
    for i=1:omega.nbElems
        map = omega.elems(i,:);
        coords = sprintf('%f,%f,0,',omega.nodes(map(:),:)');
        lengths = sprintf('%f,',node_map(map(:)));
        fprintf(fID,['ST(' coords(1:end-1) '){' lengths(1:end-1) '};\n']);
    end
    fprintf(fID,'};');
    fclose(fID);

end