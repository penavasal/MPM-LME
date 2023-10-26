classdef Material
    properties
        type
        
        E
        rhos
        nu
        gravitymodulus

        area
        ZCauchyold=[0,0,0;0,0,0;0,0,0];
        ZCauchynew=[0,0,0;0,0,0;0,0,0];

        ZDeformgradold=[1,0,0;0,1,0;0,0,1];
        ZDeformgradnew=[1,0,0;0,1,0;0,0,1];
        
        ZEpsilonOld=[0,0,0;0,0,0;0,0,0];
        ZEpsilonNew=[0,0,0;0,0,0;0,0,0];

        ZUold=[0;0];
        ZUnew=[0;0];

        %MPM PIC y FLIP

        %alpha=0 es PIC puro. Es muy disipativo pero es estable 
        %alpha=1 es FLIP puro combination. No disipa, pero es más inestable que PIC. 
        alpha=0.99;
    end
    properties (Dependent)
        Mass
        lambda
        mu
    end
    methods
        function obj=Material(MP,const)
            if nargin > 0
                obj = Material.empty(size(MP,2),0);
                for i=1:size(MP,2)
                     obj(i).area=MP(i).area;
                     mat=MP(i).mat;
                     obj(i).type=const(mat).type;
                     obj(i).gravitymodulus=const(mat).gravitymodulus;
                     obj(i).rhos=const(mat).rhos;

                     if const(mat).type=="Linear Elastic"
                         obj(i).E=const(mat).E;
                         obj(i).nu=const(mat).nu;
                     end

                end
            end
        end
        function a = get.lambda(obj)
            a=(obj.E*obj.nu)/((1+obj.nu)*(1-2*obj.nu));%N/m^2
        end
        function a = get.mu(obj)
            a=obj.E/(2*(1+obj.nu));%N/m^2
        end
        function a = get.Mass(obj)
            a=obj.area*obj.rhos;
        end
        function Courant=courant(obj,ht)
            Courant=ht*sqrt(obj.E/obj.rhos)/sqrt(obj.area);
        end

        function [obj, MP]=update(obj, MP, Unew,Uold)
            for a=1:size(obj,2) 
                %obj(a).ZPiola1old=obj(a).ZPiola1new;
                obj(a).ZCauchyold=obj(a).ZCauchynew;
                obj(a).ZDeformgradold=obj(a).ZDeformgradnew;
                obj(a).ZEpsilonOld=obj(a).ZEpsilonNew;

                MP(a).coords = MP(a).coords + (Unew(:,a) - Uold(:,a));
            end
        end

        function nodalinternalf=setInternalForces(nnodes,nearpoint,GradNinterpolanLME,mat,IDu)
            nodalinternalf=zeros(2*nnodes,1);        
            for a=1:nnodes
                lennodes=size(nearpoint{a},2);
                if lennodes
                    for i=1:lennodes
                        pto=nearpoint{a}(i);        
                        nodalinternalf(IDu(1:2,a),1) = nodalinternalf(IDu(1:2,a),1) - ...     %El signo - en nodalinternalf viene de la forma habitual de trabajar con el MPM.
                            mat(pto).ZCauchynew(1:2,1:2)*GradNinterpolanLME(:,a,pto)*mat(pto).area;
                            %mat(pto).ZPiola1old(1:2,1:2)*GradNinterpolanLME(:,a,pto)*mat(pto).area;
                    end   
                end
            end
        end

        function nodalexternalf=setExternalForces(nnodes,nearpoint,NinterpolanLME,mat,IDu,ti)

            nodalexternalf=zeros(2*nnodes,1);
            for a=1:nnodes
                lennodes=size(nearpoint{a},2);
                if lennodes
                    for i=1:lennodes
                        pto=nearpoint{a}(i);
                        g=mat(pto).gravitymodulus;

                        if ti<10
                            ff=0.5*g*(1+sin(2*ti*pi/20-pi/2));
                        else
                            ff=g;
                        end
       
                        nodalexternalf(IDu(1:2,a),1) = nodalexternalf(IDu(1:2,a),1)+...
                            NinterpolanLME(a,pto)*mat(pto).Mass*[0;-ff];
                            %NinterpolanLME(a,pto)*mat(pto).Mass*[0;-mat(pto).gravitymodulus];
                    end   
                end
            end
        end

        function nodalmomentum=setNodalMomentum(nnodes,nearpoint,NinterpolanLME,mat,IDu,ZV)

            nodalmomentum=zeros(2*nnodes,1);
            for a=1:nnodes
                lennodes=size(nearpoint{a},2);
                if lennodes
                    for i=1:lennodes
                        pto=nearpoint{a}(i);
                        nodalmomentum(IDu(1:2,a),1) = nodalmomentum(IDu(1:2,a),1)+...
                            NinterpolanLME(a,pto)*mat(pto).Mass*ZV(1:2,pto);

                    end   
                end
            end
        end

        function Zvpnew=PVelo(ZV, Zvpold,NinterpolanLME, mat, MP,IDu)

        % PIC interpola sólo la velocidad nodal en n+1
        % FLIP interpola la diferencia de velocidades nodales (VIn+1)-(VIn) 
            nptgauss=size(mat,2);
            PIC=zeros(2,nptgauss);
            FLIP=zeros(2,nptgauss);
            Zvpnew=zeros(2,nptgauss);
            for pto=1:nptgauss
                   for a=1:size(MP(pto).nears,2)
                       PIC(1:2,pto) = PIC(1:2,pto) + NinterpolanLME(MP(pto).nears(a),pto) * ...
                           ZV(IDu(1:2,MP(pto).nears(a)),2);
                       FLIP(1:2,pto) = FLIP(1:2,pto) + NinterpolanLME(MP(pto).nears(a),pto) * ...
                           (ZV(IDu(1:2,MP(pto).nears(a)),2)-ZV(IDu(1:2,MP(pto).nears(a)),1));
                   end
                   Zvpnew(:,pto)=mat(pto).alpha*(Zvpold(:,pto)+FLIP(1:2,pto))+(1-mat(pto).alpha)*PIC(1:2,pto);
            end
        end

        function [Zvpnew,Zapnew,Zupnew]=PVelo_PC(ZV, Zvpold,Zupold,NinterpolanLME, mat, MP,IDu,ht)

        % PIC interpola sólo la velocidad nodal en n+1
        % FLIP interpola la diferencia de velocidades nodales (VIn+1)-(VIn) 
            nptgauss=size(mat,2);
            PIC=zeros(2,nptgauss);
            FLIP=zeros(2,nptgauss);
            aux=zeros(2,nptgauss);
            Zvpnew=zeros(2,nptgauss);
            Zupnew=zeros(2,nptgauss);
            Zapnew=zeros(2,nptgauss);
            for pto=1:nptgauss
                   for a=1:size(MP(pto).nears,2)
                       PIC(1:2,pto) = PIC(1:2,pto) + NinterpolanLME(MP(pto).nears(a),pto) * ...
                           ZV(IDu(1:2,MP(pto).nears(a)),2);
                       FLIP(1:2,pto) = FLIP(1:2,pto) + NinterpolanLME(MP(pto).nears(a),pto) * ...
                           (ZV(IDu(1:2,MP(pto).nears(a)),2)-ZV(IDu(1:2,MP(pto).nears(a)),1));
                   end
                       
                   Zvpnew(:,pto)=mat(pto).alpha*(Zvpold(:,pto)+FLIP(:,pto))+(1-mat(pto).alpha)*PIC(:,pto);
                   aux(:,pto)=(Zvpnew(:,pto) - Zvpold(:,pto))/ht;
                   Zapnew(:,pto)=aux(:,pto);
                   Zupnew(:,pto)=Zupold(:,pto)+ht*Zvpnew(:,pto)+0.5*ht*ht*aux(:,pto);

            end
        end

        function [obj,Zupnew] = MatState(obj,ZV,Zupold,GradNinterpolanLME,NinterpolanLME,MP,ht,IDu)
            Zupnew=zeros(2,size(obj,2));
 
            for pto=1:size(obj,2)
                Fpunto=zeros(2,2);
                Ee=zeros(3);
                F=eye(3);

                aux=zeros(2,1);
                for a=1:size(MP(pto).nears,2)
                    aux=aux+ht*NinterpolanLME(MP(pto).nears(a),pto)*ZV(IDu(1:2,MP(pto).nears(a)),2); 
                    Fpunto(1:2,1:2)=Fpunto(1:2,1:2)+ZV(IDu(1:2,MP(pto).nears(a)),2)*GradNinterpolanLME(:,MP(pto).nears(a),pto)';
                end
                
                Zupnew(:,pto)=Zupold(:,pto)+aux;           
                % Total lagrangian
                %F(1:2,1:2)=obj(pto).ZDeformgradold(1:2,1:2)+ht*Fpunto(1:2,1:2);
                %
                % 
                % Updated  
                F(1:2,1:2)=(eye(2)+ ht*Fpunto(1:2,1:2)) * obj(pto).ZDeformgradold(1:2,1:2);
                obj(pto).ZDeformgradnew=F;
                
                if const(mat).type=="Linear Elastic"
                    % Small Strain
                    H(1:2,1:2)=F(1:2,1:2)-eye(2);
                    Ee(1:2,1:2)=0.5*(H(1:2,1:2)+H(1:2,1:2)');
                    obj(pto).ZEpsilonNew = Ee;

                    Eeo=obj(pto).ZEpsilonOld;

                    dEe=Ee - Eeo;
                    obj(pto).area=obj(pto).area*(1+trace(dEe));
                else
                
                    Jpto=det(F);
                    B=F*F';
                
                    obj(pto).area = MP(pto).area*Jpto;
                end

                %%%% CONSTITUTIVE MODELLING %%%%%%%

                if const(mat).type=="Linear Elastic"
                    I1=trace(Ee);
                    obj(pto).ZCauchynew = obj(pto).lambda*I1*eye(3)+2*obj(pto).mu*Ee;
                elseif const(mat).type=="Neo-Hookean"
                    %obj(pto).ZPiola1new=(obj(pto).lambda/2)*(Jpto^2-1)*eye(3)\F'+obj(pto).mu*(F-inv(F'));
                    obj(pto).ZCauchynew=(obj(pto).lambda/2/Jpto)*(Jpto^2-1)*eye(3) + obj(pto).mu/Jpto*(B-eye(3));
                end

                
                %%%%%%%
                   
                
             end
        end

        function obj = MatStatePC(obj,ZV,GradNinterpolanLME,MP,ht,IDu)

            for pto=1:size(obj,2)
                Fpunto=zeros(2,2);
                Ee=zeros(3);
                F=eye(3);
                for a=1:size(MP(pto).nears,2)
                    Fpunto(1:2,1:2)=Fpunto(1:2,1:2)+ZV(IDu(1:2,MP(pto).nears(a)),1)*GradNinterpolanLME(:,MP(pto).nears(a),pto)';
                end
                            
                % Total lagrangian
                %F(1:2,1:2)=obj(pto).ZDeformgradold(1:2,1:2)+ht*Fpunto(1:2,1:2);
                %
                % 
                % Updated  
                F(1:2,1:2)=(eye(2)+ ht*Fpunto(1:2,1:2)) * obj(pto).ZDeformgradold(1:2,1:2);
                obj(pto).ZDeformgradnew=F;
                    

                if obj(pto).type == "Linear Elastic"
                    % Small Strain 
                    H(1:2,1:2)=F(1:2,1:2)-eye(2);
                    Ee(1:2,1:2)=0.5*(H(1:2,1:2)+H(1:2,1:2)');
                    obj(pto).ZEpsilonNew = Ee;

                    Eeo=obj(pto).ZEpsilonOld;

                    dEe=Ee - Eeo;
                    obj(pto).area=obj(pto).area*(1+trace(dEe));

                else
                
                    Jpto=det(F);
                    B=F*F';
                
                    obj(pto).area = MP(pto).area*Jpto;
                end

                %%%% CONSTITUTIVE MODELLING %%%%%%%

                if obj(pto).type == "Linear Elastic"
                    I1=trace(Ee);
                    obj(pto).ZCauchynew = obj(pto).lambda*I1*eye(3)+2*obj(pto).mu*Ee;
                else

                    %obj(pto).ZPiola1new=(obj(pto).lambda/2)*(Jpto^2-1)*eye(3)\F'+obj(pto).mu*(F-inv(F'));
                    obj(pto).ZCauchynew=(obj(pto).lambda/2/Jpto)*(Jpto^2-1)*eye(3) + obj(pto).mu/Jpto*(B-eye(3));

                end

                
                %%%%%%%
                   
                
             end
        end

    end
end








