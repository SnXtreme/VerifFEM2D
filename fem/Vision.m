function Vision
figure('Name','Solution');
x=deform(omega,u,1./max(abs(u)));
  subplot(2,3,1)
    plotElemField(x,sigma);
    xlabel('x');
    ylabel('y');
    colorbar;
    title('Stresses');
  subplot(2,3,2)
    plotNodeField(x,sigmal);
    xlabel('x');
    ylabel('y');
    colorbar;
    title('Linear Stresses ZZ1 Local');
  subplot(2,3,3)
    plotNodeField(x,sigmal_2);
    xlabel('x');
    ylabel('y');
    colorbar;
    title('Linear Stresses ZZ1 Global');
   subplot(2,3,4)
    plotNodeField(x,sigmal_3);
    xlabel('x');
    ylabel('y');
    colorbar;
    title('Linear Stresses ZZ2');
figure('Name','Valeur de ZZ1');
   subplot(1,3,1)
    plotElemField(x,ZZ1_Value);
    xlabel('x');
    ylabel('y');
    colorbar;
    title('Compare ZZ1 Local');
   subplot(1,3,2)
    plotElemField(x,ZZ1_Value2);
    xlabel('x');
    ylabel('y');
    colorbar;
    title('Compare ZZ1 Global');
   subplot(1,3,3)
    plotElemField(x,ZZ2_Value);
    xlabel('x');
    ylabel('y');
    colorbar;
    title('Compare ZZ2');
figure('Name','tk_int');
   subplot(1,2,1)
    plotElemField(x,tk(:,1));
    xlabel('x');
    ylabel('y');
    colorbar;
    title('tk_int(x)');
   subplot(1,2,2)
    plotElemField(x,tk(:,2));
    xlabel('x');
    ylabel('y');
    colorbar;
    title('tk_int(y)');

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