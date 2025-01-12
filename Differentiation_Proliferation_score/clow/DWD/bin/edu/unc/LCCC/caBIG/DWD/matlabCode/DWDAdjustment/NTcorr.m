%%************************************************************************
%% NTcorr: corrector step for the NT direction. 
%%
%% SDPT3: version 3.1 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 15 Sep 2004
%%************************************************************************

  function [dX,dy,dZ] = NTcorr(blk,A,par,rp,Rd,sigmu,hRd,...
            dX,dZ,coeff,L,X,Z); 

    global matlabversion printlevel iter
    global matfct_options solve_ok 
%%
    [rhs,EinvRc]  = NTrhsfun(blk,A,par,X,Z,rp,Rd,sigmu,hRd,dX,dZ);
    m = length(rp); ncolU = size(coeff.mat12,2); 
    rhs = [rhs; zeros(m+ncolU-length(rhs),1)]; 
%%
    solve_ok = 1; resnrm = norm(rhs);
    if strcmp(matfct_options,'chol') | strcmp(matfct_options,'spchol')
       [xx,resnrm,solve_ok] = symqmr(coeff,rhs,L);
       if (solve_ok<=0) & (printlevel)
          fprintf('\n  warning: symqmr fails: %3.1f.',solve_ok); 
       end
    else
       [xx,resnrm,solve_ok] = mybicgstab(coeff,rhs,L);       
       if (solve_ok<=0) & (printlevel)
          fprintf('\n  warning: mybicgstab fails: %3.1f.',solve_ok); 
       end
    end
    if (printlevel>=3); fprintf(' %2.0d',length(resnrm)-1); end
%%   
    if (any(isnan(xx)) | any(isinf(xx)))
       solve_ok = 0;
       fprintf('\n  NTcorr: dy contains NaN or inf.');
    end
    [dX,dy,dZ] = NTdirfun(blk,A,par,Rd,EinvRc,xx); 
%%************************************************************************


