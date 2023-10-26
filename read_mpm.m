function [coordenadas,conectividad,DATA,MP,NODE_LIST]=read_mpm(filename)

    fid=fopen(filename,'rt');
    for i=1:5
        aux=fgets(fid);
    end
    aux=fgets(fid,7);  nnodes=fscanf(fid,'%d \n',[1,1]); %número de nodos en la malla
    aux=fgets(fid,7); nnodeselem=fscanf(fid,'%d \n',[1,1]); %número de nodos por elemento.
    aux=fgets(fid,7); ndim=fscanf(fid,'%d \n',[1,1]); %dimensiones del problema
    aux=fgets(fid,7);  ne=fscanf(fid,'%d \n',[1,1]); %número de elementos en la malla
    aux=fgets(fid,7);  nmats=fscanf(fid,'%d \n',[1,1]); %número de elementos en la malla
    aux=fgets(fid,7);  pgausselem=fscanf(fid,'%d \n',[1,1]); %número de puntos de gauss por elemento
                                                              %En los cuadrilateros pueden ser 4 o 9.
    for i=1:3
        aux=fgets(fid);
    end    
    t=14;
    %filas de conectividad= contenido numeración global/posición numeración local 
    %columnas de conectividad= elementos empleados en la discretización
    A=fscanf(fid,'%f \n',[nnodeselem+2,ne]);
    conectividad=A(3:(nnodeselem+2),:);
    material=A(2,:);
    %filas de coordenadas= coordenas x e y de cada nodo
    %columnas de coordenadas= nodos de la discretización (numeración global) 
    aux=fgets(fid);aux=fgets(fid);
    t=t+ne+2;
    A=fscanf(fid,'%f \n',[ndim+1,nnodes]);
    coordenadas=A(2:ndim+1,:);
    t=t+nnodes+4;

    % Condiciones de contorno
    formato = '%s %s %s %s %s %s %s %s %s %s %s'; % formato de cada línea 

    tsargs = {...
        'HeaderLines',0,...
        'HeaderColumns',0,...
        'ReturnOnError',false,...
        'EmptyValue',0,...
        'CollectOutput',true,...
        'EndOfLine','\r\n'};

    % Collect info
    res  = textscan(fid,formato,-1,tsargs{:});
    data=res{1};
    
    fclose(fid);
    [l,~]=size(data);

    bcs=0;
    pls=0;
    lls=0;
    vls=0;
    abcs=0;
    rbs=0;
    rbs_nds=[];
    BC={};
    PL={};
    ABC={};
    LL={};
    VL={};
    RB={};

    fin=0;

    t=4;
    while fin==0
        t=t+1;
        s1=data{t,1};
        switch s1
            case 'BC_SET'
                num=str2double(data{t,3});
                if num~=0
                    bcs=bcs+1;
                    if bcs==str2double(data{t,2})
                        list=zeros(num,1);
                        for i=1:num
                            t=t+1;
                            list(i)=str2double(data{t,1});
                        end
                        BC{bcs}=list;
                        clear num list
                    else
                        fprintf('Error, bad numbering of set of BC!!\n')
                        stop
                    end
                end
                continue
            case 'RIGID_BODY'
                num=str2double(data{t,3});
                if num~=0
                    rbs=rbs+1;
                    if rbs==str2double(data{t,2})
                        list=zeros(num,3);
                        for i=1:num
                            t=t+1;
                            list(i,1)=x(str2double(data{t,1}),1);
                            list(i,2)=x(str2double(data{t,1}),2);
                            list(i,3)=x(str2double(data{t,2}),1);
                            list(i,4)=x(str2double(data{t,2}),2);
                        end
                        rbs_nds=union(rbs_nds,[str2double(data{t,1}) str2double(data{t,2})]);
                        RB{rbs}=list;
                        clear num list
                    else
                        fprintf('Error, bad numbering of set of LOADs!!\n')
                        stop
                    end
                end
                continue
            case 'ABSORBING_BC'
                num=str2double(data{t,3});
                if num~=0
                    abcs=abcs+1;
                    if abcs==str2double(data{t,2})
                        list=zeros(num,3);
                        for i=1:num
                            t=t+1;
                            list(i,1)=str2double(data{t,1});
                            list(i,2)=str2double(data{t,2});
                            if NNE==8 || NNE==6
                              list(i,3)=str2double(data{t,3});  
                            end
                        end
                        ABC{abcs}=list;
                        clear num list
                    else
                        fprintf('Error, bad numbering of set of LOADs!!\n')
                        stop
                    end
                end
                continue
            case 'POINT_LOAD'
                num=str2double(data{t,3});
                if num~=0
                    pls=pls+1;
                    if pls==str2double(data{t,2})
                        list=zeros(num,1);
                        for i=1:num
                            t=t+1;
                            list(i)=str2double(data{t,1});
                        end
                        PL{pls}=list;
                        clear num list
                    else
                        fprintf('Error, bad numbering of set of LOADs!!\n')
                        stop
                    end
                end
                continue
            case 'LINE_LOAD'
                num=str2double(data{t,3});
                if num~=0
                    lls=lls+1;
                    if lls==str2double(data{t,2})
                        list=zeros(num,3);
                        for i=1:num
                            t=t+1;
                            list(i,1)=str2double(data{t,1});
                            list(i,2)=str2double(data{t,2});
                            if NNE==8 || NNE==6
                              list(i,3)=str2double(data{t,3});  
                            end
                        end
                        LL{lls}=list;
                        clear num list
                    else
                        fprintf('Error, bad numbering of set of LOADs!!\n')
                        stop
                    end
                end
                continue
            case 'VOLUME_LOAD'
                num=str2double(data{t,3});
                if num~=0
                    vls=vls+1;
                    if vls==str2double(data{t,2})
                        list=zeros(num,1);
                        for i=1:num
                            t=t+1;
                            list(i)=str2double(data{t,1});
                        end
                        VL{vls}=list;
                        clear num list
                    else
                        fprintf('Error, bad numbering of set of LOADs!!\n')
                        stop
                    end 
                end
                continue
            case 'END_BC'
                t=t+2;
                if t>l
                    fin=1;
                elseif strcmp(data{t,1},'LOADS')==0
                    fprintf('Error, bad reading of geometry!!\n')
                    stop
                end
                continue
            case 'END_LDS'
                fin=1;
                continue
            otherwise
                fprintf('Error, bad reading of geometry!!\n')
                stop
        end
    end
    
    if isempty(rbs_nds)==0
       [a,b]=size(elem);
       for i=length(rbs_nds):-1:1
           x=[x(1:rbs_nds(i)-1,:);x(rbs_nds(i)+1:end,:)];
           for j=1:a
               for k=1:b
                   if elem(j,k)>rbs_nds(i)
                       elem(j,k)=elem(j,k)-1;
                   end
               end
           end
           %BCS
           for j=1:bcs
               list=BC{j};
               for k=1:length(list)
                   if list(k)>rbs_nds(i)
                        list(k)=list(k)-1;
                   end
               end 
               BC{j}=list;
           end
           %PLS
           for j=1:pls
               list=PL{j};
               for k=1:length(list)
                   if list(k)>rbs_nds(i)
                        list(k)=list(k)-1;
                   end
               end
               PL{j}=list;
           end
           %VLS
           for j=1:vls
               list=VL{j};
               for k=1:length(list)
                   if list(k)>rbs_nds(i)
                        list(k)=list(k)-1;
                   end
               end  
               VL{j}=list;
           end
           %LLS
           for j=1:lls
               list=LL{j};
               [a1,b1]=size(list);
               for k=1:a1
                   for l=1:b1
                        if list(k,l)>rbs_nds(i)
                            list(k,l)=list(k,l)-1;
                        end
                   end
               end
               LL{j}=list;
           end
       end
    end

    NODE_LIST=struct('bcs',bcs,'pls',pls,'lls',lls,'vls',vls,'rbs',rbs,...
        'abcs',abcs,'BC',{BC},'PL',{PL},'LL',{LL},'VL',{VL},'RB',{RB},...
        'ABC',{ABC});


    %aux=fgets(fid);aux=fgets(fid);aux=fgets(fid);aux=fgets(fid);
    %BC{1}=[];
    %aux=fgets(fid,9);  bcs=fscanf(fid,'%d \n',[1,1]);
    %A=fscanf(fid,'%f \n',[1,bcs]);
    %BC{1}=A;
    %fclose(fid);       
    
    % Puntos de gauss, tantos bucles como materiales.
    p=0;
    for i=1:nmats
        h=0;
        for e=1:ne
            if material(e)==i
               h=h+1;
               conect1(:,h)=conectividad(:,e);
               elemPG(h)=e;
               Area(h)=area(nnodeselem,conect1(:,h),coordenadas');
            end
        end

        coordptog=ptogauss(coordenadas,conect1,h,pgausselem);%coordenadas de los ptos de gauss

        for e=1:h
            for j=1:pgausselem
                p=p+1;
                MP(p).mat=i;
                MP(p).el=elemPG(e);
                MP(p).coords=coordptog(:,p);
                MP(p).area=Area(e)/pgausselem;
                if nnodeselem==4
                    MP(p).spacing = sqrt(MP(p).area);
                else
                    disp('elemento no identificado')
                    stop
                end
                MP(p).upc=[];
            end
        end
        
        clear elemPG conect1

    end


    DATA.ne=ne;
    DATA.nnodes=nnodes;
    DATA.ndim=ndim;
    DATA.NNE=nnodeselem;
    DATA.PGE=pgausselem;
    DATA.PG=p;

end

function Area=area(NNE,elem,x_a)

        xn=zeros(NNE+1,1);
        yn=zeros(NNE+1,1);

        for i=1:NNE
            nd=elem(i);
            xn(i)=x_a(nd,1);
            yn(i)=x_a(nd,2);
        end
        nd=elem(1);
        xn(i+1)=x_a(nd,1);
        yn(i+1)=x_a(nd,2);

        a1 = 0;
        a2 = 0;

        for i=1:NNE
            a1 = a1 + xn(i)*yn(i+1);
            a2 = a2 + yn(i)*xn(i+1);
        end
        Area=abs(a1-a2)/2;

end

function punto=ptogauss(coord,nodelement,ne,pgausselem)

%Empleamos elementos isoparamétricos para establecer la posición de los
%puntos de gauss X(ptog)=Nhat_alpha(ptogr)*X_alpha donde 
%X(ptog) son las coordenadas del punto de gauss en el elemento físico en la configuración de referencia
%Nhat_alpha(ptog) son las funciones de forma en el elemento de referencia(isoparamétrico) evaluadas en las posiciones de los puntos de gauss en ese elemento de referencia.
%X_alpha son las coordenadas de los nodos de cada elemento físico en la configuración de referencia

punto=zeros(2,pgausselem*ne);
%filas de punto= coordenas x e y de cada punto de gauss
%columnas de punto= puntos de gauss en toda la malla. Se van rellenando
%elemento a elemento. 
%Primer elmento -> todos los puntos de gauss
%Segundo elemento -> todos los puntos de gauss
%y así sucesivamente

if abs(size(nodelement,1)-6)<1e-10
%puntos de gauss para elementos de 6 nodos
    
    if pgausselem==3 %tres puntos de gauss por elemento

    
        Gpd2con3=[1/6 1/6
          2/3 1/6
          1/6 2/3];
   
    
        for i=1:ne
            for pto=1:pgausselem 
                % X(ptog)=Nhat_alpha(ptogr)*X_alpha    
                punto(1:2,pto+(i-1)*pgausselem)=(1-Gpd2con3(pto,1)-Gpd2con3(pto,2))*(2*(1-Gpd2con3(pto,1)-Gpd2con3(pto,2))-1)*coord(:,nodelement(1,i))+...
                                    Gpd2con3(pto,1)*(2*Gpd2con3(pto,1)-1)*coord(:,nodelement(2,i))+...
                                    Gpd2con3(pto,2)*(2*Gpd2con3(pto,2)-1)*coord(:,nodelement(3,i))+...
                                    4*Gpd2con3(pto,1)*(1-Gpd2con3(pto,1)-Gpd2con3(pto,2))*coord(:,nodelement(4,i))+...
                                    4*Gpd2con3(pto,1)*Gpd2con3(pto,2)*coord(:,nodelement(5,i))+...
                                    4*Gpd2con3(pto,2)*(1-Gpd2con3(pto,1)-Gpd2con3(pto,2))*coord(:,nodelement(6,i));
                               
         
            end  
        end

    elseif pgausselem==7

        a=(6+sqrt(15))/21;
        b=(4/7)-a;

        Gpd2con7=[1/3 1/3
         a a
         1-2*a a
         a 1-2*a
         b b
         1-2*b b
         b 1-2*b];
        for i=1:ne
            for pto=1:pgausselem 
     % X(ptog)=Nhat_alpha(ptogr)*X_alpha      
                punto(1:2,pto+(i-1)*pgausselem)=(1-Gpd2con7(pto,1)-Gpd2con7(pto,2))*(2*(1-Gpd2con7(pto,1)-Gpd2con7(pto,2))-1)*coord(:,nodelement(1,i))+...
                                    Gpd2con7(pto,1)*(2*Gpd2con7(pto,1)-1)*coord(:,nodelement(2,i))+...
                                    Gpd2con7(pto,2)*(2*Gpd2con7(pto,2)-1)*coord(:,nodelement(3,i))+...
                                    4*Gpd2con7(pto,1)*(1-Gpd2con7(pto,1)-Gpd2con7(pto,2))*coord(:,nodelement(4,i))+...
                                    4*Gpd2con7(pto,1)*Gpd2con7(pto,2)*coord(:,nodelement(5,i))+...
                                    4*Gpd2con7(pto,2)*(1-Gpd2con7(pto,1)-Gpd2con7(pto,2))*coord(:,nodelement(6,i));
             
            end
        end
    end

elseif abs(size(nodelement,1)-4)<1e-10

    if pgausselem==4
        Gpd2con4=[-1/sqrt(3) -1/sqrt(3)
          1/sqrt(3) -1/sqrt(3)
          1/sqrt(3) 1/sqrt(3)
          -1/sqrt(3) 1/sqrt(3)];
   
    
        for i=1:ne
            for pto=1:pgausselem 
                % X(ptog)=Nhat_alpha(ptogr)*X_alpha      
                punto(1:2,pto+(i-1)*pgausselem)=((1-Gpd2con4(pto,1))*(1-Gpd2con4(pto,2))/4)*coord(:,nodelement(1,i))+...
                                    ((1+Gpd2con4(pto,1))*(1-Gpd2con4(pto,2))/4)*coord(:,nodelement(2,i))+...
                                    ((1+Gpd2con4(pto,1))*(1+Gpd2con4(pto,2))/4)*coord(:,nodelement(3,i))+...
                                    ((1-Gpd2con4(pto,1))*(1+Gpd2con4(pto,2))/4)*coord(:,nodelement(4,i));
                                    
                                
                            
         
            end  
        end

    elseif pgausselem==9

    

        Gpd2con9=[0 0
           -sqrt(3/5) -sqrt(3/5)
              0 -sqrt(3/5)
           sqrt(3/5) -sqrt(3/5)
           sqrt(3/5) 0
           sqrt(3/5) sqrt(3/5)
           0 sqrt(3/5)
           -sqrt(3/5) sqrt(3/5)
           -sqrt(3/5) 0];
       
        for i=1:ne
            for pto=1:pgausselem 
                % X(ptog)=Nhat_alpha(ptogr)*X_alpha      
                punto(1:2,pto+(i-1)*pgausselem)=((1-Gpd2con9(pto,1))*(1-Gpd2con9(pto,2))/4)*coord(:,nodelement(1,i))+...
                                    ((1+Gpd2con9(pto,1))*(1-Gpd2con9(pto,2))/4)*coord(:,nodelement(2,i))+...
                                    ((1+Gpd2con9(pto,1))*(1+Gpd2con9(pto,2))/4)*coord(:,nodelement(3,i))+...
                                    ((1-Gpd2con9(pto,1))*(1+Gpd2con9(pto,2))/4)*coord(:,nodelement(4,i));
             
            end
        end
    end 
 
    elseif abs(size(nodelement,1)-8)<1e-10

        if pgausselem==4
            Gpd2con4=[-1/sqrt(3) -1/sqrt(3)
          1/sqrt(3) -1/sqrt(3)
          1/sqrt(3) 1/sqrt(3)
          -1/sqrt(3) 1/sqrt(3)];
   
    
for i=1:ne
    for pto=1:pgausselem 
     % X(ptog)=Nhat_alpha(ptogr)*X_alpha      
    punto(1:2,pto+(i-1)*pgausselem)=(-(1-Gpd2con4(pto,1))*(1-Gpd2con4(pto,2))*(1+Gpd2con4(pto,1)+Gpd2con4(pto,2))/4)*coord(:,nodelement(1,i))+...
                                    (-(1+Gpd2con4(pto,1))*(1-Gpd2con4(pto,2))*(1-Gpd2con4(pto,1)+Gpd2con4(pto,2))/4)*coord(:,nodelement(2,i))+...
                                    (-(1+Gpd2con4(pto,1))*(1+Gpd2con4(pto,2))*(1-Gpd2con4(pto,1)-Gpd2con4(pto,2))/4)*coord(:,nodelement(3,i))+...
                                    (-(1-Gpd2con4(pto,1))*(1+Gpd2con4(pto,2))*(1+Gpd2con4(pto,1)-Gpd2con4(pto,2))/4)*coord(:,nodelement(4,i))+...
                                    ((1-Gpd2con4(pto,1)^2)*(1-Gpd2con4(pto,2))/2)*coord(:,nodelement(5,i))+...
                                    ((1+Gpd2con4(pto,1))*(1-Gpd2con4(pto,2)^2)/2)*coord(:,nodelement(6,i))+...
                                    ((1-Gpd2con4(pto,1)^2)*(1+Gpd2con4(pto,2))/2)*coord(:,nodelement(7,i))+...
                                    ((1-Gpd2con4(pto,1))*(1-Gpd2con4(pto,2)^2)/2)*coord(:,nodelement(8,i));
                               
         
    end  
end

elseif pgausselem==9

    

    Gpd2con9=[0 0
          -sqrt(3/5) -sqrt(3/5)
              0 -sqrt(3/5)
           sqrt(3/5) -sqrt(3/5)
           sqrt(3/5) 0
           sqrt(3/5) sqrt(3/5)
           0 sqrt(3/5)
           -sqrt(3/5) sqrt(3/5)
           -sqrt(3/5) 0];
       
   for i=1:ne
    for pto=1:pgausselem 
     % X(ptog)=Nhat_alpha(ptogr)*X_alpha      
    punto(1:2,pto+(i-1)*pgausselem)=(-(1-Gpd2con9(pto,1))*(1-Gpd2con9(pto,2))*(1+Gpd2con9(pto,1)+Gpd2con9(pto,2))/4)*coord(:,nodelement(1,i))+...
                                    (-(1+Gpd2con9(pto,1))*(1-Gpd2con9(pto,2))*(1-Gpd2con9(pto,1)+Gpd2con9(pto,2))/4)*coord(:,nodelement(2,i))+...
                                    (-(1+Gpd2con9(pto,1))*(1+Gpd2con9(pto,2))*(1-Gpd2con9(pto,1)-Gpd2con9(pto,2))/4)*coord(:,nodelement(3,i))+...
                                    (-(1-Gpd2con9(pto,1))*(1+Gpd2con9(pto,2))*(1+Gpd2con9(pto,1)-Gpd2con9(pto,2))/4)*coord(:,nodelement(4,i))+...
                                    ((1-Gpd2con9(pto,1)^2)*(1-Gpd2con9(pto,2))/2)*coord(:,nodelement(5,i))+...
                                    ((1+Gpd2con9(pto,1))*(1-Gpd2con9(pto,2)^2)/2)*coord(:,nodelement(6,i))+...
                                    ((1-Gpd2con9(pto,1)^2)*(1+Gpd2con9(pto,2))/2)*coord(:,nodelement(7,i))+...
                                    ((1-Gpd2con9(pto,1))*(1-Gpd2con9(pto,2)^2)/2)*coord(:,nodelement(8,i));
             
    end
   end
 end   
elseif abs(size(nodelement,1)-9)<1e-10

 if pgausselem==4
 Gpd2con4=[-1/sqrt(3) -1/sqrt(3)
          1/sqrt(3) -1/sqrt(3)
          1/sqrt(3) 1/sqrt(3)
          -1/sqrt(3) 1/sqrt(3)];
   
    
for i=1:ne
    for pto=1:pgausselem 
     % X(ptog)=Nhat_alpha(ptogr)*X_alpha      
    punto(1:2,pto+(i-1)*pgausselem)=((1/4)*(Gpd2con4(pto,1)^2-Gpd2con4(pto,1))*(Gpd2con4(pto,2)^2-Gpd2con4(pto,2)))*coord(:,nodelement(1,i))+...
                                    ((1/4)*(Gpd2con4(pto,1)^2+Gpd2con4(pto,1))*(Gpd2con4(pto,2)^2-Gpd2con4(pto,2)))*coord(:,nodelement(2,i))+...
                                    ((1/4)*(Gpd2con4(pto,1)^2+Gpd2con4(pto,1))*(Gpd2con4(pto,2)^2+Gpd2con4(pto,2)))*coord(:,nodelement(3,i))+...
                                    ((1/4)*(Gpd2con4(pto,1)^2-Gpd2con4(pto,1))*(Gpd2con4(pto,2)^2+Gpd2con4(pto,2)))*coord(:,nodelement(4,i))+...
                                    ((1/2)*(Gpd2con4(pto,2)^2-Gpd2con4(pto,2))*(1-Gpd2con4(pto,1)^2))*coord(:,nodelement(5,i))+...
                                    ((1/2)*(Gpd2con4(pto,1)^2+Gpd2con4(pto,1))*(1-Gpd2con4(pto,2)^2))*coord(:,nodelement(6,i))+...
                                    ((1/2)*(Gpd2con4(pto,2)^2+Gpd2con4(pto,2))*(1-Gpd2con4(pto,1)^2))*coord(:,nodelement(7,i))+...
                                    ((1/2)*(Gpd2con4(pto,1)^2-Gpd2con4(pto,1))*(1-Gpd2con4(pto,2)^2))*coord(:,nodelement(8,i))+...
                                    (((1-Gpd2con4(pto,1)^2)*(1-Gpd2con4(pto,2)^2)))*coord(:,nodelement(9,i));
                               
                     
                                
                                
                                
    end  
end

elseif pgausselem==9

    

    Gpd2con9=[0 0
          -sqrt(3/5) -sqrt(3/5)
              0 -sqrt(3/5)
           sqrt(3/5) -sqrt(3/5)
           sqrt(3/5) 0
           sqrt(3/5) sqrt(3/5)
           0 sqrt(3/5)
           -sqrt(3/5) sqrt(3/5)
           -sqrt(3/5) 0];
       
   for i=1:ne
    for pto=1:pgausselem 
     % X(ptog)=Nhat_alpha(ptogr)*X_alpha      
    punto(1:2,pto+(i-1)*pgausselem)=((1/4)*(Gpd2con9(pto,1)^2-Gpd2con9(pto,1))*(Gpd2con9(pto,2)^2-Gpd2con9(pto,2)))*coord(:,nodelement(1,i))+...                                     
                                    ((1/4)*(Gpd2con9(pto,1)^2+Gpd2con9(pto,1))*(Gpd2con9(pto,2)^2-Gpd2con9(pto,2)))*coord(:,nodelement(2,i))+...                                     
                                    ((1/4)*(Gpd2con9(pto,1)^2+Gpd2con9(pto,1))*(Gpd2con9(pto,2)^2+Gpd2con9(pto,2)))*coord(:,nodelement(3,i))+...                                     
                                    ((1/4)*(Gpd2con9(pto,1)^2-Gpd2con9(pto,1))*(Gpd2con9(pto,2)^2+Gpd2con9(pto,2)))*coord(:,nodelement(4,i))+...                                   
                                    ((1/2)*(Gpd2con9(pto,2)^2-Gpd2con9(pto,2))*(1-Gpd2con9(pto,1)^2))*coord(:,nodelement(5,i))+...                                     
                                    ((1/2)*(Gpd2con9(pto,1)^2+Gpd2con9(pto,1))*(1-Gpd2con9(pto,2)^2))*coord(:,nodelement(6,i))+...                                    
                                    ((1/2)*(Gpd2con9(pto,2)^2+Gpd2con9(pto,2))*(1-Gpd2con9(pto,1)^2))*coord(:,nodelement(7,i))+...                                     
                                    ((1/2)*(Gpd2con9(pto,1)^2-Gpd2con9(pto,1))*(1-Gpd2con9(pto,2)^2))*coord(:,nodelement(8,i))+...                                     
                                    (((1-Gpd2con9(pto,1)^2)*(1-Gpd2con9(pto,2)^2)))*coord(:,nodelement(9,i));
                                  
                                
                                 
          

                
          
             
    end
   end
 end    
 
end
end


