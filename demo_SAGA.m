              



load('data.mat')
        SAGA_K = 30;
        lam0 =.1;
        lam1 =.1;
        lam2 =.1;
            
        train_data=disagg_Signal;
                  
        parm_AGS.K = SAGA_K;
        parm_AGS.lambda_0 = lam0;
        parm_AGS.lambda_1 = lam1;
        parm_AGS.lambda_2 = lam2;
        parm_AGS.lambda_3 =1;  
        parm_AGS.rho_1= parm_AGS.lambda_1;
        parm_AGS.rho_0=parm_AGS.lambda_0;
        parm_AGS.iter=500;
        parm_AGS.mask = [];
        parm_AGS.A_hard=1;
        

        time_SAGA=[];
        tic;
        [objs_DFT,U_DFT,A_DFT]=SAGA_ortho_l1_Phi_ortho(train_data,L,PhiDFT,parm_AGS);
        time_SAGA=[time_SAGA,toc]
  



        time_SAGA=[];
        tic;
        [objs_REMI,U_REMI,A_REMI]=SAGA_ortho_l1_Phi_not_ortho(train_data,L,PhiREMI,parm_AGS);
        time_SAGA=[time_SAGA,toc]
        
        
        time_SAGA=[];
        tic;
        [objs_SPLIN,U_SPLIN,A_SPLIN]=SAGA_ortho_l1_Phi_not_ortho(train_data,L,PhiSPLIN,parm_AGS);
        time_SAGA=[time_SAGA,toc]



              
           

            
