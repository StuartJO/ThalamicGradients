function DisplayCorr(X,Y)

[RHO,P,CIL,CIU] = corrcoef(X,Y,'Rows','complete');
    
disp(['Pearson''s r(',num2str(length(X)-2),') = ',num2str(RHO(1,2)),', p = ',num2str(P(1,2)),', CI = [',num2str(CIL(1,2)),', ',num2str(CIU(1,2)),'], two-tailed'])
       