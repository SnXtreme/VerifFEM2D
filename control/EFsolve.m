function [u,sigma,varargout]=EFsolve(omega,B,dtop,Fd,varargin)

% Ecriture du systeme
[K,M] = FEMMat(omega,B);
F = FEMVec(dtop,Fd);

% Initialisation de l'inconnue
u = zeros(size(K,2),1);

% Imposition des conditions aux limites
[cl_index,u0] = CL(varargin{:});
u(cl_index) = u0;
   
% Resolution
u(~cl_index) = K(~cl_index,~cl_index)\(F(~cl_index) - K(~cl_index,cl_index)*u(cl_index));

% Calcul des champs de contraintes
sigma=stress(omega,B,u); %Non lissé

if nargout>2
    varargout{1}=K;
    varargout{2}=M;
end

end