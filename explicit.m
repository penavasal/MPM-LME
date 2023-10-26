function explicit(ti,tfin,tplot,ht,coordenadas,MP,DATA,mat,IDu,ContourU,...
    V0,vp0,up0,ap0)

[ZV,Zvpold,Zupold,~]=initial(mat,ti,V0,vp0,up0,ap0);

nnodes = DATA.nnodes;
nptgauss=size(MP,2);

nodal=State(coordenadas);

NinterpolanLME=zeros(nnodes,nptgauss);
GradNinterpolanLME=zeros(2,nnodes,nptgauss);

while ti<tfin

    [NinterpolanLME,GradNinterpolanLME,MP,nearpoint]=...
        ShapeandGradShapeLME(coordenadas,MP,DATA,Zupold,NinterpolanLME,GradNinterpolanLME);
    nodal=setMass(nodal,NinterpolanLME,nearpoint,mat);
     
    % Balance de FUERZAS y MOMENTO LINEAL nodal.
    %                -----------------------------------------------
    % Pasamos información de las partículas a los nodos. Siempre que hay transferencia 
    % empleamos ensamblaje, para poder contabilizar la participación en un nodo concreto 
    % de la información de TODAS las partículas que están alrededor suya

    nodalmomentum=setNodalMomentum(nnodes,nearpoint,NinterpolanLME,mat,IDu,Zvpold);
    nodalexternalf=setExternalForces(nnodes,nearpoint,NinterpolanLME,mat,IDu);
    nodalinternalf=setInternalForces(nnodes,nearpoint,GradNinterpolanLME,mat,IDu);

    %Calculamos momento lineal para el tiempo n+1 en nodos 
    nodalmomentum(:,2) = nodalmomentum(:,1) + ht*(nodalexternalf(:,1)+nodalinternalf(:,1));
    %aplicamos condiciones dirichlet en nodos
    nodalmomentum(IDu(1,ContourU{1}),2)=0;
    nodalmomentum(IDu(2,ContourU{2}),2)=0;

   
    % Obtenemos VELOCIDADES en tiempo n+1 a partir del momento lineal en tiempo n+1 
    %                 ----------------------
    % No hay ensamblaje. Simplemente empleamos las estructuras de almacenamiento IDu. 
    % Idea: Como es velocidades en nodos a partir de momento en nodos no es necesario ensamblaje

    ZV=NodalVelo(nodal,nodalmomentum,IDu,ZV);
    
    % Con las velocidades nodales en n y n+1 obtenemos por interpolación la nueva 
    % velocidad en particulas. FLIP y PIC

    Zvpnew=PVelo(ZV, Zvpold, NinterpolanLME, mat, MP,IDu);

    % Calculamos nuevamente el MOMENTO LINEAL en tiempo n+1.
    %                                    -----------------------------
    % Con estas nuevas velocidades en particulas volvemos a transferir información a los nodos. 
    % Como lo hacemos por ensamblaje, hay que borrar lo que teníamos previamente

    nodalmomentum(:,2)=setNodalMomentum(nnodes,nearpoint,NinterpolanLME,mat,IDu,Zvpnew);
   
    % Aplicamos condiciones dirichlet
    nodalmomentum(IDu(1,ContourU{1}),2)=0;
    nodalmomentum(IDu(2,ContourU{2}),2)=0;

    % Obtenemos VELOCIDADES en tiempo n+1 a partir del nuevo momento lineal en tiempo n+1 
    %                 ----------------------
    % No hay ensamblaje. Simplemente empleamos las estructuras de almacenamiento IDu. 
    % Idea: Como es velocidades en nodos a partir de momento en nodos no es necesario ensamblaje

    ZV=NodalVelo(nodal,nodalmomentum,IDu,ZV);
    
    % Con la información de los desplazamientos nodales a tiempo n+1
    % obtenemos TENSION (1ºPK) y DEFORMACION (F) a tiempo n+1
    %               ---------------             --------------------------  
    % Transferimos información de los desplazamientos nodales a desplazamientos de partículas

    [mat,Zupnew] = MatState(mat,ZV,Zupold,GradNinterpolanLME,NinterpolanLME,MP,ht,IDu);
   
  
    % Update and Save files
    
    ti=ti+ht;
    
    %ZU(:,1)=ZU(:,2);
    ZV(:,1)=ZV(:,2);
    Zvpold(:,:)=Zvpnew(:,:);
    [mat, MP]=update(mat, MP, Zupnew, Zupold(:,:));
    Zupold(:,:)=Zupnew(:,:);
    

    if mod(ti,tplot)<ht
        print_file('tiempo.asc',ti);
        %print_file('desp.asc', ZU(:,2)')
        print_file('velocidades.asc', ZV(:,2)')
        print_file('up.asc', Zupnew(:,:));
        print_file('vp.asc', Zvpnew(:,:));
        %print_file_tens('Cauchy.asc', 'Deformationgradient.asc',mat);
        print_file_tens('Cauchy.asc', 'Epsilon.asc',mat);
        
        fprintf('%i\t%i\n',ti,ht);
    end
  
end


end