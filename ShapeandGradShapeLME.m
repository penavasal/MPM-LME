function [Shape,DShape,MP,nearpoint]=ShapeandGradShapeLME(coordenadas,MP,DATA,Up,Shape,DShape)

    global LME

   Na=DATA.nnodes;
   allnodes=1:Na;
   dimension=DATA.ndim;
   
   Np=DATA.PG;
   allpoints=1:Np;
   
   
   % LME parameters
   %gamma=7.0;
   %target_zero=1.0e-5; 
   %TOL=1e-6;
   %node_spacing=0.0015;% measure of the nodal spacing
   %node_spacing=0.01;% measure of the nodal spacing


   gamma=LME.gamma;
   target_zero=LME.target_zero; 
   TOL=LME.TOL;


   Nap=zeros(Np,1);
   for i=1:Np

        node_spacing(i)=MP(i).spacing;% measure of the nodal spacing
        
        Beta(i)=gamma/node_spacing(i)^2;
        range1=node_spacing(i)*sqrt(-1/gamma*log(target_zero));
        range_lme(i)=max(range1, 2*node_spacing(i));

        coordptog=MP(i).coords;
        aux=[];
        for a=1:Na       
            if norm(coordptog-coordenadas(:,a),2)<range_lme(i)
                aux=[aux,allnodes(a)];
            end
        end
        nearnode(i)={aux};
        MP(i).nears=aux;
        Nap(i,1)=length(nearnode{i});
   end
   
   for a=1:Na
        aux=[];
        for i=1:Np
            if ismember(a,nearnode{i})
                aux=[aux,allpoints(i)];
            end
        end
        nearpoint(a)={aux};
   end
   
   
   for p=1:Np
        si=0;
        upc=MP(p).upc;
        upg=Up(:,p);
        
        if isempty(upc)
            si=1;
        else
            up=sqrt((upc(1)-upg(1))^2+(upc(2)-upg(2))^2)/sqrt(MP(p).area);
            if up>0.1
                si=1;
            else
                si=0;
            end
        end
        
        if si
            Shape(:,p)=0;
            DShape(:,:,p)=0;
            coordptog=MP(p).coords;        
            lambda=zeros(dimension,1);
            for ite=1:50
                Z=0;
                for a=1:Nap(p,1)
                    Z=Z+exp(-Beta(i)*(coordptog-coordenadas(:,nearnode{p}(a)))'*(coordptog-coordenadas(:,nearnode{p}(a)))+...
                        lambda'*(coordptog-coordenadas(:,nearnode{p}(a))));
                end

                pa=zeros(Nap(p,1),1);
                for a=1:Nap(p,1)
                    pa(a,1)=(1/Z)*exp(-Beta(i)*(coordptog-coordenadas(:,nearnode{p}(a)))'*(coordptog-coordenadas(:,nearnode{p}(a)))+...
                        lambda'*(coordptog-coordenadas(:,nearnode{p}(a))));
                end

                R=zeros(2,1);
                for a=1:Nap(p,1)
                    R(:,1)=R(:,1)+(coordptog-coordenadas(:,nearnode{p}(a)))*pa(a,1);  
                end

                aux=zeros(2,2);
                for a=1:Nap(p,1)
                    aux=aux+pa(a,1)*((coordptog-coordenadas(:,nearnode{p}(a)))*(coordptog-coordenadas(:,nearnode{p}(a)))');        
                end 

                DR=aux-R(:,1)*R(:,1)';
                if abs(rcond(DR))<1e-8
                    rcond(DR);
                    %disp('Newton Failed, near to singular matrix')
                end 

                dpa=zeros(2,Nap(p,1));
                for a=1:Nap(p,1)
                    dpa(:,a)=-DR\(pa(a,1)*(coordptog-coordenadas(:,nearnode{p}(a))));
                end

                %fprintf('%i\t%i\t%i\n',p,ite,norm(R,2));
                if norm(R,2)<TOL 

                    Shape(nearnode{p},p)=pa(:,1);
                    DShape(:,nearnode{p},p)=dpa(:,:);

                    break
                end

                lambda=lambda-DR\R(:,1);
                %lambda=lambda-(DR+norm(R(:,1),2)*eye(2))\R(:,1);
            end

            if abs(ite-50)<1e-5
                error('max ite alcanzado');
            end
            
            MP(p).upc=Up(:,p);
        end
      
   end

   
