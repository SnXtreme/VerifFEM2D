function Visu_Error(omega,u,sigma,sigmal,error,time,varargin)
% Visu_Error permet la visualisation des champs lissés permettant le calcul
% d'erreur, ainsi que l'erreur et le champ h d'adaptation pour comparer les
% différents critères d'erreur
% ARGIN :
% omega : Class Mesh comportant le maillage complet
% u : Champ de déplacement
% sigma : Champ de contrainte
% sigmal : Cell contenant les champs lissés sous la forme
% Par ligne, {'Nom du critère' 'Champ'}
% error : Cell contenant les champs d'erreur
% varargin : Pour eet - on attend sigma_adm,omega_adm,u_adm

nbE=omega.nbElems;
nbN=omega.nbNodes;
x=deform(omega,u,1./max(abs(u)));
nb_sigmal=size(sigmal,1);
nb_error=size(error,1);

%% Figure contenant les champs lissés
figure('Name','Solution');
  subplot(floor((nb_sigmal+1)/3)+1,mod((nb_sigmal+1),3)+1,1)
    plotElemField(x,sigma);
    xlabel('x');
    ylabel('y');
    axis equal tight;
    colorbar;
    title('Stresses');
  for id_sigmal=1:nb_sigmal
    subplot(floor((1+nb_sigmal)/3)+1,mod((1+nb_sigmal),3)+1,id_sigmal+1)
    if mod(length(sigmal{id_sigmal,2}),nbE)==0
        plotElemField(x,sigmal{id_sigmal,2});
    elseif mod(length(sigmal{id_sigmal,2}),nbN)==0
        plotNodeField(x,sigmal{id_sigmal,2});
    else
        warning('Mauvaise définition de sigmal? : Impossible de tracer le champ lissé');
    end
    xlabel('x');
    ylabel('y');
    axis equal tight;
    colorbar;
    title(['Linear Stresses ' sigmal{id_sigmal,1}]);
  end
  
%% Figure contenant les champs d'erreur  
figure('Name','Error');
  for id_error=1:nb_error
    subplot(floor(nb_error/3)+1,mod(nb_error,3)+1,id_error)
    if mod(length(error{id_error,2}),nbE)==0
        plotElemField(x,error{id_error,2});
    elseif mod(length(error{id_error,2}),nbN)==0
        plotNodeField(x,error{id_error,2});
    else
        warning('Mauvaise définition de VErreur? : Impossible de tracer le champ d''erreur');
    end
    xlabel('x');
    ylabel('y');
    axis equal tight;
    colorbar;
    title(['Erreur ' error{id_error,1}]);
  end

%% Histogramme en temps
figure('Name',['Temps de calcul pour ',num2str(omega.nbElems),' Elements'])
bar(time);
  
%% Pour EET

if nargin>7
    %Cas EET
    %vargin={sigma_adm,omega_adm,u_adm}

sigma_adm=varargout{1};
omega_adm=varargout{2};
u_adm=varargout{3};
n = max(u);
figure('Name','Interpolation');
  subplot(1,2,1);
    plot(omega.border);
    hold on
    quiver(omega.nodes(:,1),omega.nodes(:,2),u(1:2:end)./n,u(2:2:end)./n,2);
    xlabel('x');
    ylabel('y');
    title('Displacement');
  subplot(1,2,2);
    plot(omega_adm.border);
    hold on
    quiver(omega_adm.nodes(:,1),omega_adm.nodes(:,2),u_adm(1:2:end)./n,u_adm(2:2:end)./n,0);
    xlabel('x');
    ylabel('y');
    title('Interpolated Displacement');

figure('Name','CRE');
  subplot(1,2,1);
    plotElemField(deform(omega,u,1./max(abs(u))),sigma);
    xlabel('x');
    ylabel('y');
    colorbar;
    title('Stresses');
  subplot(1,2,2);
    plotAdmField(deform(omega_adm,u_adm,1./max(abs(u_adm))),sigma_adm);
    xlabel('x');
    ylabel('y');
    colorbar;
    title('Admissible Stresses');
end
end