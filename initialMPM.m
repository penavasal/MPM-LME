

function [ZV,Zvpold,Zupold,Zapold]=initialMPM(mat,tinic,V0,vp0,up0,ap0)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% FILE PREPARATION
    %Abrimos y cerramos ficheros para guardar resultados numericos
    
    %Primary variables at nodes
    Vfile = fopen('velocidades.asc','w');
    
    % Particles
    upfile = fopen('up.asc','w');   %desplazamiento particulas
    vpfile = fopen('vp.asc','w');   %velodicad particulas
    
    fclose(Vfile);fclose(upfile);fclose(vpfile);
    
    %Secundary variables at gauss points
    Cauchyfile=fopen('Cauchy.asc','w');
    fclose(Cauchyfile);
    Deformationgradientfile = fopen('Deformationgradient.asc','w');
    fclose(Deformationgradientfile ); 
    Deformationfile = fopen('Epsilon.asc','w');
    fclose(Deformationfile );

    %Other variables
    tiempofile=fopen('tiempo.asc','w');
    fclose(tiempofile);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% INITIAL CONDITIONS
    
    %Condiciones iniciales de las Primary and secondary variables 
    print_file('velocidades.asc', V0');

    print_file('up.asc', up0);

    print_file('vp.asc', vp0);
    
    %print_file_tens('Cauchy.asc', 'Deformationgradient.asc',mat);
    print_file_tens('Cauchy.asc', 'Epsilon.asc',mat);
    
    print_file('tiempo.asc', tinic);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Variables auxiliares para actualizar valores con el tiempo
    %las variables Z(:,1) ó que terminen en old se refieren al tiempo n, mientras que las variables
    %Z(:,2) ó que terminen en new se refieren al tiempo n+1.
    
    %ZU=zeros(size(U0,1),2);increU=zeros(2*nnodes,2);
    %ZU(:,1)=U0;
    
    ZV=zeros(size(V0,1),2);
    ZV(:,1)=V0;
    
    %ZA=zeros(size(A0,1),2);
    %ZA(:,1)=A0;
    
    Zvpold=vp0;
    %Zvpnew=Zvpold;
    Zupold=up0;
    %Zupnew=Zupold;
    Zapold=ap0;
    %Zapnew=Zapold;






end