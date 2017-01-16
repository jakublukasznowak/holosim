% zmiana w matlabie
function c55 = mymesh (a,b)

Na=size(a,1); Ma=size(a,2);
Nb=size(b,1); Mb=size(b,2);
N=max([Na Nb]);

if Ma>1 || Mb>1
    a=num2cell(a,2);
    b=num2cell(b,2);
    c3=cell([Na ceil(Na/Nb)*Nb]);
else
    c3=nan([Na ceil(Na/Nb)*Nb]);
end


%if Na>Nb
    c1=cat(3,repmat(a,1,ceil(Na/Nb)*Nb),repmat(b',Na,ceil(Na/Nb)));
    for i=1:N
        temp=circshift(c1,i-1+ceil(Na/Nb)*Nb,2);
        try
            c3(:,i,1)=diag(temp(:,:,1));
            c3(:,i,2)=diag(temp(:,:,2));
        catch
            for j=1:Na
                c3(j,i,1)=temp(j,j,1);
                c3(j,i,2)=temp(j,j,2);
            end
        end
    end
% %else
%     c1=cat(3,repmat(a,1,Nb),repmat(b',Na,1));
%     for i=1:N
%         temp=circshift(c1,i-1+N,2);
%         try
%             c3(:,i,1)=diag(temp(:,:,1));
%             c3(:,i,2)=diag(temp(:,:,2));
%         catch
%             for j=1:Na
%                 c3(j,i,1)=temp(j,j,1);
%                 c3(j,i,2)=temp(j,j,2);
%             end
%         end
%     end
% end

c55=reshape(c3(:,1:Nb,:),[Na*Nb,2]);

end


