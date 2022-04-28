
function [objs,U,A]=SAGA_ortho_l1_Phi_not_ortho(X,L,Phi,opts);
    
%% parameters
K=opts.K;


lambda_1=opts.lambda_1;
lambda_0=opts.lambda_0;
rho_1=opts.rho_1;
rho_0=opts.rho_0;
MAX_ITER=opts.iter;


%% initialization
[n,t]=size(X);

PhiPhiT=Phi*Phi';
Phi2_eye=eye(size(PhiPhiT));


[p,t]=size(Phi);

U=eye([K,p]);

A=eye(n,K);

Z=A;
V=U; 

Gamma_1=0;
Gamma_0=0;


obj_old = 0;
objs =[];


% 
%for z update
k_eye=eye(K,K);
S_A_eye=eye(size(A));


if issparse(L)
    
[QA, DA] = eigs( opts.lambda_2*2*L,length(L));
dA = diag(DA);
else
[QA, dA] = eig( opts.lambda_2*2*L, 'vector');

end



D=X;
for i = 1:MAX_ITER
    
%   A-update
     B=U*Phi;

  %% A with L
     A = pre_compute_sylvester(QA,dA,(B*B'+k_eye*rho_0),(D*B'+rho_0*Z-Gamma_0));   

    
    %% posetive
      A = (A>0).*A;     
      
    %% orthnormal
    [U_A,S_A,V_A]=svd(A);
     A=U_A*S_A_eye*V_A';

      %% U-update

     % if A is orthnormal    
     U=(2*A'*D*Phi'+rho_1*V-Gamma_1)*inv(2*PhiPhiT+Phi2_eye*rho_1);  



    
%     %V update
    h_1=U-Gamma_1/rho_1;
    V = sign(h_1).*max(abs(h_1)-lambda_1/rho_1,0);
    
    
%     Z update
     h_2=A-Gamma_0/rho_0;
     Z = sign(h_2).*max(abs(h_2)-lambda_0/rho_0,0);

    if ~isempty(opts.mask)
    P=A*U*Phi;
    
    D=update_d(P,X,opts.mask,opts.lambda_3);
    else 
        D=X;   
    end


    % update Gamma
     Gamma_1 = Gamma_1 + rho_1*(V-U);
     Gamma_0 = Gamma_0 + rho_0*(Z-A);

    rho_1=min(rho_1*1.1, 1e5);

    
    % stop condition
    % ++++++++++++++++++++++++++++ stop condition check
    if (0==mod(i,100))
        [obj_new,term0,term1,term2,term3]  = getObj(D,X,A,U,Phi,L,opts);
        
        
        objs.total(i/10)=obj_new;
        objs.term0(i/10)=term0;
        objs.term1(i/10)=term1;
        objs.term2(i/10)=term2;
        objs.term3(i/10)=term3;


        residual = abs(obj_old-obj_new); %/obj_new;
        disp(['iteration=',num2str(i),', obj=',num2str(obj_new),...
            ',','residual-',num2str(i),'=',num2str(residual),' , term0=', ...
            num2str(term0),', term1=',num2str(term1),', term2=',num2str(term2),', term3=',num2str(term3)]);
        fprintf('\n')

        if residual <  1e-5
            break;
        else
            obj_old = obj_new;
        end
    end
    % update
    
    
end
% save result

end


function [obj_new,term0,term1,term2,term3] = getObj(D,X,A,U,Phi,L,opts)


term0=norm(D-A*U*Phi);
term1=opts.lambda_0*norm(A,1);
term2=opts.lambda_1*norm(U,1);
term3=opts.lambda_2*trace(A'*L*A);

obj_new = [term1+term2+term3+term0];




end