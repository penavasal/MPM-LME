function postprocess(str2)

    load DataProblem.mat

    load tiempo.asc
    %load desp.asc
    load up.asc
    load vp.asc
    load Epsilon.asc
    %load Deformationgradient.asc
    %load Piola1.asc
    load Cauchy.asc

    
    nptgauss=DATA.PG;
    nnodes=DATA.nnodes;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %A partir de aqui es todo postproceso dentro de Matlab
    
    VonMisesptog=zeros(length(tiempo),nptgauss);
    VonMisesnodes=zeros(length(tiempo),nnodes);
    Meanstresptog=zeros(length(tiempo),nptgauss);
    Meanstresnodes=zeros(length(tiempo),nnodes);
    Cauchyptog11=zeros(length(tiempo),nptgauss);
    Cauchy_vec=zeros(nptgauss,length(tiempo)*4);
    Cauchynodes11=zeros(length(tiempo),nnodes);
    %Jptog=zeros(length(tiempo),nptgauss);
    Jnodes=zeros(length(tiempo),nnodes);
    Masa=zeros(nnodes,nnodes);
    RHSrecu=zeros(nnodes,1);
    for n=1:length(tiempo)
        for npto=1:nptgauss
            %despparticu=up((1:2)+(n-1)*2,npto);
            %FF=Deformationgradient((1:3)+(n-1)*3,(1:3)+(npto-1)*3);
            %JJ=det(FF);
            %PP=Piola1((1:3)+(n-1)*3,(1:3)+(npto-1)*3);
            %Cauchy=(1/JJ)*PP*FF';
            CC=Cauchy((1:3)+(n-1)*3,(1:3)+(npto-1)*3);
            Cauchyptog11(n,npto)=CC(1,1);
            dev=CC-(1/3)*trace(CC)*eye(3);
            %Jptog(n,npto)=JJ;
            VonMisesptog(n,npto)=sqrt(3/2)*sqrt(trace(dev'*dev));
            Meanstresptog(n,npto)=(1/3)*trace(CC);
            Cauchy_vec(npto,(n-1)*4+1)=CC(1,1);
            Cauchy_vec(npto,(n-1)*4+2)=CC(2,2);
            Cauchy_vec(npto,(n-1)*4+3)=CC(3,3);
            Cauchy_vec(npto,(n-1)*4+4)=CC(1,2);
        end
    end
    

    figure
    %hold on
    %plot(tiempo,desp(:,IDu(2,conectividad(4,1))))
    %Zupnew(:,(1:pgausselem)+(i-1)*pgausselem)
    plot(tiempo,vp(1+((1:length(tiempo))-1)*2, 1),'r')
    %plot(tiempo,up(2+((1:length(tiempo))-1)*2, 363),'r')
    %hold off
    %min(desp(:,IDu(2,1)))
    
    %%% PLOT MATLAB
    %plotMATLAB(coordenadas,conectividad,DATA,MP,tiempo,up);
    
    %%% PLOT VTK
    meshname='mesh.vtk';
    CELL_vtk(coordenadas,conectividad,meshname)
    for tn=1:size(tiempo,1)
        %str2='block';
        mpFileName=(['VTK/' str2 '_' num2str(tn) '.vtk']);
        
        uppost=zeros(2,nptgauss);
        U=zeros(2,nptgauss);
        V=zeros(2,nptgauss);
        S=zeros(nptgauss,4);
        S(:,1:4)=Cauchy_vec(:,(tn-1)*4+1:tn*4);
        for i=1:nptgauss
            uppost(:,i)=MP(i).coords+up((1:2)+(tn-1)*2,i);
            U(:,i)=up((1:2)+(tn-1)*2,i);
            V(:,i)=vp((1:2)+(tn-1)*2,i);
        end
        writeVTK(DATA,uppost,mpFileName,U,V,S);
    end

end

function writeVTK(DATA,coords,mpFileName,U,V,S)
    
    %VTK output file generation: material point data
    %--------------------------------------------------------------------------
    % Author: Pedro Navas (after Mian Xie(after William Coombs))
    % Date:   31/08/2023
    % Description:
    % Function to generate a VTK file containing the material point data.
    nmp = DATA.PG;
    nD  = DATA.ndim;

    fid=fopen(mpFileName,'wt');
    fprintf(fid,'# vtk DataFile Version 2.0\n');
    fprintf(fid,'MATLAB generated vtk file, MX\n');
    fprintf(fid,'ASCII\n');
    fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
    fprintf(fid,'POINTS %i double\n',nmp);
    if nD==3
        for i=1:nmp
            fprintf(fid,'%f %f %f \n',coords(:,i)');%mpD.x(i,:));
        end
    elseif nD==2
        for i=1:nmp
            fprintf(fid,'%f %f %f\n',coords(:,i)',0);
        end
    elseif nD==1
        for i=1:nmp
            fprintf(fid,'%f %f %f\n',coords(:,i)',0,0);
        end
    end
    fprintf(fid,'\n');

    fprintf(fid,'POINT_DATA %i\n',nmp);
        %du
    fprintf(fid, 'VECTORS d_u float \n');
    for j=1:nmp
        fprintf(fid,'%f %f %f\n',U(:,i)',0);
    end
    %Vu
    fprintf(fid, 'VECTORS v_u float \n');
    for j=1:nmp
        fprintf(fid,'%f %f %f\n',V(:,i)',0);
    end
    
    %STRESS
    fprintf(fid,'SCALARS sigma_xx FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    for i=1:nmp
        fprintf(fid,'%f\n',S(i,1));
    end
    fprintf(fid,'\n');

    fprintf(fid,'SCALARS sigma_yy FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    for i=1:nmp
        fprintf(fid,'%f\n',S(i,2));
    end
    fprintf(fid,'\n');

    fprintf(fid,'SCALARS sigma_zz FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    for i=1:nmp
        fprintf(fid,'%f\n',S(i,3));
    end
    fprintf(fid,'\n');

    fprintf(fid,'SCALARS sigma_xy FLOAT %i\n',1);
    fprintf(fid,'LOOKUP_TABLE default\n');
    for i=1:nmp
        fprintf(fid,'%f\n',S(i,4));
    end
    fprintf(fid,'\n');
    fclose('all');

end

function CELL_vtk(coord,etpl,meshName)
    
    %VTK output file generation: mesh data
    %--------------------------------------------------------------------------
    % Author: William Coombs
    % Date:   15/01/2019
    % Description:
    % Function to generate a VTK file containing the background mesh data.
    %
    %--------------------------------------------------------------------------
    % MAKEVTK(coord,etpl,meshName)
    %--------------------------------------------------------------------------
    % Input(s):
    % coord    - coordinates of the grid nodes (nodes,nD)
    % etpl     - element topology (nels,nen) 
    % meshName - VTK file name, for example 'mesh.vtk'  
    %--------------------------------------------------------------------------
    % See also:
    % 
    %--------------------------------------------------------------------------

    [nD,nodes]=size(coord);
    [nen,nels]=size(etpl);

    %% FEM etpl to VTK format
    if nD ==3
        if nen==20
            tvtk=[1 7 19 13 3 5 17 15 8 12 20 9 4 11 16 10 2 6 18 14];
            elemId=25;
            elemFormat='%i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n';
        elseif nen==8
            tvtk=[1 4 8 5 2 3 7 6];
            elemId=12;
            elemFormat='%i %i %i %i %i %i %i %i %i\n';
        elseif nen==10
            tvtk=[1 2 3 4 5 6 7 8 10 9];
            elemId=24;
            elemFormat='%i %i %i %i %i %i %i %i %i %i %i\n';
        elseif nen==4
            tvtk=[1 3 2 4];
            elemId=10;
            elemFormat='%i %i %i %i %i\n';
        elseif nen==9
            tvtk=[3 1 7 5 2 8 6 4 9];
            elemId=10;
            elemFormat='%i %i %i %i %i %i %i %i %i %i\n';
        end
    elseif nD==2
        if nen==3
            tvtk=[1 3 2];
            elemId=5;
            elemFormat='%i %i %i %i\n';
        elseif nen==4
            tvtk=[1 4 2 3];
            elemId=8;
            elemFormat='%i %i %i %i %i\n';
        elseif nen==8
            tvtk=[1 7 5 3 8 6 4 2];
            elemId=23;
            elemFormat='%i %i %i %i %i %i %i %i %i\n';
        end
    end

    %% Generation of vtk file
    fid=fopen(meshName,'wt');
    fprintf(fid,'# vtk DataFile Version 2.0\n');
    fprintf(fid,'MATLAB generated vtk file, WMC\n');
    fprintf(fid,'ASCII\n');
    fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
    fprintf(fid,'POINTS %i double\n',nodes);
    if nD==3
        for i=1:nodes
            fprintf(fid,'%f %f %f \n',coord(:,i));
        end
    elseif nD==2
        for i=1:nodes
            fprintf(fid,'%f %f %f \n',coord(:,i),0);
        end
    end
    fprintf(fid,'\n');
    fprintf(fid,'CELLS %i %i\n',nels,(nen+1)*nels);
    for i=1:nels
        fprintf(fid,elemFormat,nen,(etpl(tvtk,i)-1));       
    end
    fprintf(fid,'\n');
    fprintf(fid,'CELL_TYPES %i\n',nels);
    for i=1:nels
        fprintf(fid,'%i\n',elemId);       
    end
    fprintf(fid,'\n');
    
end

function plotMATLAB(coordenadas,conectividad,DATA,MP,tiempo,up)

    %geometría y malla
    ax=min(coordenadas(1,:));bx=max(coordenadas(1,:));
    ay=min(coordenadas(2,:));by=max(coordenadas(2,:));
    
    X = zeros(DATA.NNE,DATA.ne) ;UX = zeros(DATA.NNE,DATA.ne) ;
    Y = zeros(DATA.NNE,DATA.ne) ;UY = zeros(DATA.NNE,DATA.ne) ;
    nptgauss=DATA.PG;
    uppost=zeros(2,nptgauss);
    
    figure 
    for tn=1:size(tiempo,1)
        if abs(rem(tn,1)-0)<1e-5%activar si hay muchos pasos de tiempo
            for i=1:nptgauss
                uppost(:,i)=MP(i).coords+up((1:2)+(tn-1)*2,i);
            end
            %ux = desp(tn,IDu(1,:));
            %uy = desp(tn,IDu(2,:));
        
            for i=1:size(conectividad,2)
            
                %Si es de 4 nodos no hace falta ordenar los nodos
                nd=conectividad(:,i);
                X(:,i)=coordenadas(1,nd)';    
                Y(:,i)=coordenadas(2,nd)';
               % UX(:,i) = ux(nd) ;
               % UY(:,i) = uy(nd) ;
               % profile(:,i) = VonMisesnodes(tn,nd) ;  
               % profile(:,i) = Meanstresnodes(tn,nd);  
               %profile(:,i) = Jnodes(tn,nd) ; 
               % profile(:,i) = uy(nd) ; 
                
            end
        
           %defoX = X+UX;
           %defoY = Y+UY;
          
           fill(X,Y,'w')
           hold on
           % plot(uppost(1,:),uppost(2,:),'KO','MarkerSize',5,'MarkerFaceColor','k');
           scatter(uppost(1,:),uppost(2,:),80,up(2+(tn-1)*2,:),'filled');
           hold off
           axis equal;
           axis([ax-ax/10 bx+bx/10 ay-ay/10 by+by/10]);
           colormap(jet)
           colorbar
           caxis([min(min(up(2+((1:length(tiempo))-1)*2,:))) max(max(up(2+((1:length(tiempo))-1)*2,:)))])
           drawnow
        end
    end

end


    %for n=1:length(tiempo)
    %for iele=1:ne
       
    %Masa(conectividad(:,iele),conectividad(:,iele))=Masa(conectividad(:,iele),conectividad(:,iele))+...
    %    stimamasa(Fhat(:,:,(1:pgausselem)+(iele-1)*pgausselem),Ninterpolantref,pgausselem,...
    %    Pesosptog,nnodeselem);  
    
    %RHSrecu(conectividad(:,iele),1)=RHSrecu(conectividad(:,iele),1)+stimaRHS(Fhat(:,:,(1:pgausselem)+(iele-1)*pgausselem),...
    %    Ninterpolantref,pgausselem,Pesosptog,nnodeselem,VonMisesptog(n,(1:pgausselem)+(iele-1)*pgausselem));       
    
    %end
    %VonMisesnodes(n,:)=(Masa\RHSrecu)';
    %Meanstresnodes(n,:)=(Masa\RHSrecu)';
    %Jnodes(n,:)=(Masa\RHSrecu)';
    %end
    
    %for n=1:length(tiempo)
    %for i=1:ne
     
    %VonMisesnodes(n,conectividad(:,i))=...
    %    linearextrapolation(Ninterpolantref,pgausselem,nnodeselem,VonMisesptog(n,(1:pgausselem)+(i-1)*pgausselem));  
    %        
    %Meanstresnodes(n,conectividad(:,i))=...
    %    linearextrapolation(Ninterpolantref,pgausselem,nnodeselem,Meanstresptog(n,(1:pgausselem)+(i-1)*pgausselem));  
    
    %Cauchynodes11(n,conectividad(:,i))=...
    %    linearextrapolation(Ninterpolantref,pgausselem,nnodeselem,Cauchyptog11(n,(2:9)+(i-1)*pgausselem));  
    %Jnodes(n,conectividad(:,i))=...
    %    linearextrapolation(Ninterpolantref,pgausselem,nnodeselem,Jptog(n,(1:pgausselem)+(i-1)*pgausselem));  
    % ponemos (2:9) por que estamos empleando elementos de 9 puntos de gauss y no tenemos en cuenta el pto central
    %para extrapolar valores a los nodos. Para el de 4 ptos de gauss en elem de
    %4 nodos usamos (1:pgausselem)
      
           
    %end
    %end
