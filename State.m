classdef State
    properties
        mass
        coords
    end
%     properties (Dependent)
%         Mass
%     end
    methods
        function obj=State(coordenadas)
            if nargin > 0
                obj = State.empty(size(coordenadas,2),0);
                for i=1:size(coordenadas,2)
                     obj(i).coords=coordenadas(:,i);
                end
            end
        end

        function obj=setMass(obj,NinterpolanLME,nearpoint,mat)
            for a=1:size(obj,2)
                obj(a).mass=0;
                for pto=1:size(nearpoint{a},2)                                                                
                    obj(a).mass=obj(a).mass+NinterpolanLME(a,nearpoint{a}(pto))*mat(pto).Mass;         
                end                            
            end
        end

        function ZV=NodalVelo(nodal,nodalmomentum,IDu,ZV)

            for a=1:size(nodal,2)
                if nodal(a).mass~=0
                    ZV(IDu(1:2,a),2)=nodalmomentum(IDu(1:2,a),2)/nodal(a).mass;  
                end
            end
        end

        function ZV=PredVelo(nodal,nnodes,nearpoint,NinterpolanLME,mat,IDu,Zvpold,Zapold,ht,gamma,ContourU)
            
            nodalmomentum=setNodalMomentum(nnodes,nearpoint,NinterpolanLME,mat,IDu,Zvpold);
            nodalInertial=setNodalMomentum(nnodes,nearpoint,NinterpolanLME,mat,IDu,Zapold);

            nodalmomentum(IDu(1,ContourU{1}),1)=0;
            nodalmomentum(IDu(2,ContourU{2}),1)=0;
            nodalInertial(IDu(1,ContourU{1}),1)=0;
            nodalInertial(IDu(2,ContourU{2}),1)=0;
            
            ZV=zeros(nnodes*2,1);
            for a=1:size(nodal,2)
                if nodal(a).mass~=0
                    ZV(IDu(1:2,a),1)=(nodalmomentum(IDu(1:2,a),1)+(1-gamma)*ht*nodalInertial(IDu(1:2,a),1))/nodal(a).mass;  
                end
            end
        end

    end
end




