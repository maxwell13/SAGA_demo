function [D] = update_d(P,X, mask, lam_3)

if ~isempty(mask)
% D=(P+lam_3*mask.*X)./(eye(size(mask))+lam_3*mask);
D=P;
%vecotrize X and P and remove missing values 
x=X(:);
x(mask)=[];
p=P(:);
p(mask)=[];
d=(p+lam_3*x)*inv(1+lam_3);
D(setdiff(1:end,mask))=d;

else
D=X;
end


end

