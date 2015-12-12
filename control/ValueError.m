function [Error,sigmal,time]=ValueError(Criterion,omega,sigma,B,domega,dtop,Fd,dleft,dbottom,od)
NbC=length(Criterion);
Coded_crit={'ZZ1_Node','ZZ1_Shape','ZZ2','ResiduExt','EET'};

Error=cell(NbC,2);
sigmal=cell(NbC,2);
time=zeros(NbC,1);
for id_crit=1:length(Criterion)
    crit=Criterion(id_crit);
    nb_crit=find(strcmp(crit,Coded_crit));
    if ~isempty(nb_crit)
        nb_crit=double(nb_crit);
        switch nb_crit
            case {1,2}
                % Cas ZZ1
                vararg={B};
            case 3
                % Cas ZZ2
                vararg={B,domega};
            case 4
                % Cas ResiduExt
                vararg={domega,dtop,Fd};
            case 5
                % Cas EET
                vararg={B,dtop,Fd,dleft,dbottom,od};
            otherwise
                % Critère non codé?
                error('Critère proposé : varargin non inscrit')
        end
    tic    
    [Error{id_crit,2},sigmal{id_crit,2}]=VError(nb_crit,omega,sigma,vararg{:});
    time(id_crit)=toc;
    Error(id_crit,1)=crit;
    sigmal(id_crit,1)=crit;
    else
        warning(['Critère ',crit,' non codé']);
    end
end
end