%% Script d'étude des estimateurs a priori
clear all;
close all force;

if exist('./save/1.mat', 'file')==0 % Si la boucle de calcul a déjà été effectuée,

addpath('./fem/');

% Parametres du probleme
E = 200;
nu = 0.3;
B = E/(1+nu)/(1-2*nu)*[1-nu nu 0; nu 1-nu 0; 0 0 (1-2*nu)/2]; % Matrix de Hook (Avec les notations de Voigt)
Fd = 1000*[0 1];
a = 1;

% Parametres de la bouche :
nh = 25. % Nombre de réduction du facteur d'échelle
nex = 60 % Facteur d'échelle définissant la solution exacte

NormeE = zeros(nh,1);
Fech = zeros(nh,1);
list = (1.:nh);
list = [list,nex];
%Boucle sur le facteur d'échelle :
for fech = list(1,:)
fech
% Charge le maillage
omega = Mesh(geo2msh('./meshes/plate_crack.geo',1,1/fech));


% Extraction de maillage secondaires
H = max(omega.nodes(:,2));
domega = omega.border;
dtop = domega.restrict(@(x) x(:,2) == H);
dleft = domega.restrict(@(x) x(:,1) == 0);
dbottom = domega.restrict(@(x) x(:,2) == 0 & x(:,1) >= a);

% Ecriture du systeme
K = FEMMat(omega,B);
F = FEMVec(dtop,Fd);

% Initialisation de l'inconnue
u = zeros(size(K,2),1);

% Imposition des conditions aux limites
[cl_index,u0] = CL(dleft,0,[1 0],dbottom,0,[0 1]);
u(cl_index) = u0;
   
% Resolution
u(~cl_index) = K(~cl_index,~cl_index)\(F(~cl_index) - K(~cl_index,cl_index)*u(cl_index));
NormeE(fech,1) = sqrt(u'*K*u); % Définition de la norme énergétique
Fech(fech,1) = max(omega.elementSize); % Expression de h

% Visualisation
figure('Name','Solution');
  subplot(1,2,1)
    plotElemField(deform(omega,u,1./max(abs(u))),stress(omega,B,u));
    xlabel('x');
    ylabel('y');
    colorbar;
    title('Stresses');
  subplot(1,2,2);
    plot(omega.border);
    hold on
    quiver(omega.nodes(:,1),omega.nodes(:,2),u(1:2:end),u(2:2:end));
    xlabel('x');
    ylabel('y');
    title('Displacement');
 
save('./save/1','u','omega','NormeE','nex', 'Fech');

end;
end;

load('./save/1');

Erreur=sqrt(NormeE(nex,1)^2-NormeE.^2);
figure('Name','Solution');
plot(log10(Fech),log10(Erreur));
npoly=floor(nex/10);
S=polyfit(log10(Fech(1:npoly)),log10(Erreur(1:npoly)),1);
hold on;
plot((-1:0.1:1),(S(1)*(-1:0.1:1)+S(2)));
S