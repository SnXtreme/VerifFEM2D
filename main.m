%% Programme d'étude des estimateurs
clearvars;
close all force;
clc;

%% Script d'exemple de resolution du probleme elements finis
addpath('./fem/');
addpath('./eet/');
addpath('./control/');

% Parametres du probleme
E = 200;
nu = 0.3;
B = E/(1+nu)/(1-2*nu)*[1-nu nu 0; nu 1-nu 0; 0 0 (1-2*nu)/2]; % Matrix de Hook (Avec les notations de Voigt)
Fd = 1000*[0 1];
a = 1;
N_adapt = 1;
filename = './meshes/plate_hole.geo';
posname = 1/10;
criterion = {'ZZ1_Node','ZZ1_Shape','ZZ2','ResiduExt'};
% criterion = {'ZZ2'};
id_adapt = 3;
isvid=false;
debug=false;
Nbelem_adapt=zeros(1,N_adapt);


% Bouclage sur l'adaptation
for id=1:N_adapt 

% Chargement du maillage
[omega,domega,dtop,dleft,dbottom]=mesh_load(filename,posname);
Nbelem_adapt(id)=omega.nbNodes;
% Résolution du problème
[u,sigma] = EFsolve(omega,B,dtop,Fd,dleft,0,[1 0],dbottom,0,[0 1]);

% Calcul de l'erreur 
[error,sigmal,time]=ValueError(criterion,omega,sigma,B,domega,dtop,Fd,dleft,dbottom,1);

% Adaptation :
NbC = length(criterion);
h_adapt=cell(NbC,1);
if N_adapt>1
for id_crit=1:NbC
    h_adapt{id_crit} = adapt(omega,error{id_crit,2},1E-3*(u'*u));
end
    posname = adaptation(omega,h_adapt{id_adapt});
end

% Visualisation
Visu_Error(omega,u,sigma,sigmal,error,time);

end