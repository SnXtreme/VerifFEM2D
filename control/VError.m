function [Value_Error,sigmal,varargout]=VError(type,omega,sigma,varargin)

switch type
    case 1
        % 'ZZ1_Node'
        %varargin{:}=B
        sigmal=stressl(omega,sigma); % Sigma lissé aux noeuds du maillage
        [Value_Error]=ZZ(omega,varargin{1},sigma,sigmal);
        
    case 2
        % 'ZZ1_Shape'
        %varargin{:}=B
        K=matrixAssembly(omega,0,eye(3),3);
        F=FEMVec2(omega,sigma);
        sigmal=K\F;
        sigmal=full(sigmal);
        [Value_Error]=ZZ(omega,varargin{1},sigma,sigmal);
        
    case 3
        % 'ZZ2'
        %varargin{:}=B,domega
        sigmal=stress_ZZ2(omega,sigma,varargin{2});
        [Value_Error]=ZZ(omega,varargin{1},sigma,sigmal);

    case 4
        % 'ResiduExt'
        % varargin{:} = dtop,Fd
        % sigmal=tk
        [Value_Error,sigmal] = ResiduExt(omega,sigma,varargin{:});

    case 5
        % 'EET'
        % varargin{:} = B,dtop,Fd,dleft,dbottom,1
        [omega_adm,sigma_adm] = eet(omega,sigma,varargin{:});
        u_adm = interpField(omega,u,omega_adm);
        if nargout>1
            sigmal=sigma_adm;
            varargout={omega_adm,u_adm};
        end
        
    otherwise
        error(['Le type d''erreur n''est pas dans la base de données \n' ...
            'Les erreurs codées sont : ZZ1_Node (1)\n,',...
            'ZZ1_Shape (2)\n,',...
            'ZZ2 (3)\n,',...
            'ResiduExt (4)\n,',...
            'EET (5)\n,']);
end
end