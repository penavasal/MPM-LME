
function [Ms,As]=MP_load(x_a,nd,elem,pge,MP)

    ndl=length(nd);
    [a,b]=size(elem);

    [pairs,les]=Cpairs(x_a,nd);
    
    El=zeros(ndl-1,1);
    MPs=zeros(ndl-1,pge);
    for j=1:ndl-1
        par=pairs(j,:);
        for e=1:a
            t=0;
            for i=1:b
                if ismember(elem(e,i),par)
                    t=t+1;
                end
            end
            if t==2
                El(j)=e;
                break;
            end
        end 
    end
    
    for j=1:ndl-1
        e=El(j);
        t=0;
        for i=1:size(MP,2)
            if MP(i).el==e
                t=t+1;
                MPs(j,t)=i;
            end
        end
    end
    for j=1:ndl-1
        [Ms,As]=distMPs(MPs,pairs,les,x_a,MP);
    end

end


function [d1,d2,i1,i2]=dist(x_a,nd,j)
    x1=x_a(nd(j),1);
    y1=x_a(nd(j),2);
    d=zeros(length(nd),1);
    for i=1:length(nd)
        if(i~=j)
            d(i)=sqrt((x_a(nd(i),1)-x1)^2+(x_a(nd(i),2)-y1)^2);
%             if strcmp(r,'X')
%                 if (x_a(nd(i),1)-x1)<0
%                     d(i)=-d(i);
%                 end
%             end
        else
            d(i)=1.0e32;
        end
    end
    [~,i]=min(abs(d));
    d1=d(i);
    i1=i;
    d(i)=1.0e32;
    [d2,i2]=min(abs(d));   
end


function [pairs,les]=Cpairs(x_a,nd)

    evalua=zeros(length(nd),4);
    for j=1:length(nd)
        [d1,d2,i1,i2]=dist(x_a,nd,j);
        evalua(j,:)=[i1,i2,d1,d2];
    end

    times=zeros(length(nd),1);
    for j=1:length(nd)
        for i=1:2
            times(evalua(j,i))=times(evalua(j,i))+1;
        end
    end

    [~,ini]=min(times);
    pairs=zeros(length(nd)-1,2);
    les=zeros(length(nd)-1,1);
    
    % First
    j=1;
    sig=evalua(ini,1);
    pairs(j,:)=[nd(ini),nd(sig)];
    les(j)=evalua(ini,3);
    for i=1:2
        if evalua(sig,i)==ini
           evalua(sig,i)=1e32;
           evalua(sig,i+2)=1e32;
        end
    end
    ini=sig;     
    % Rest
    for j=2:length(nd)-1
        sig=min(evalua(ini,1:2));
        pairs(j,:)=[nd(ini),nd(sig)];
        les(j)=min(evalua(ini,3:4));
        for i=1:2
            if evalua(sig,i)==ini
                evalua(sig,i)=1e32;
                evalua(sig,i+2)=1e32;
            end
        end
        ini=sig;
    end

end



function [Ms,As]=distMPs(MPs,pairs,les,x_a,MP)
    
    [a,b]=size(MPs);

    x1=x_a(pairs(1),1);
    y1=x_a(pairs(1),2);
    x2=x_a(pairs(2),1);
    y2=x_a(pairs(2),2);
    
    pb=sqrt(b);
    aa=les/pb;
    t=0;
    Ms=zeros(pb*a,1);
    As=zeros(pb*a,1);
    for i=1:a
        d=zeros(b,1);
        for j=1:b
            mp=MPs(i,j);
            d1=sqrt((MP(mp).coords(1)-x1)^2+(MP(mp).coords(2)-y1)^2);
            d2=sqrt((MP(mp).coords(1)-x2)^2+(MP(mp).coords(2)-y2)^2);
            d(j)=min(d1,d2);
        end
        for j=1:pb
            [~,m]=min(abs(d));
            t=t+1;
            Ms(t)=MPs(i,m);
            As(t)=aa(i);
            d(m)=1.0e32;
        end
    end
end