function [omega,domega,dtop,dleft,dbottom,dright]=mesh_load(filename,posname)

% Charge le maillage
if nargin==1;
    omega = Mesh(geo2msh(filename,1,1/10));
elseif nargin>1 && isnumeric(posname)
    omega = Mesh(geo2msh(filename,1,posname));
else
    omega = Mesh(geo2msh(filename,1,posname));
end

% Extraction de maillage secondaires
domega = omega.border;
dtop = domega.restrict(@(x) x(:,2) == max(omega.nodes(:,2)));
dleft = domega.restrict(@(x) x(:,1) == 0);
dright = domega.restrict(@(x) x(:,1) == max(omega.nodes(:,2)));
dbottom = domega.restrict(@(x) x(:,2) == 0);