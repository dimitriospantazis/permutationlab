%% check if test is exact (error control should be 0.05 under null hypothesis)

numpermutation = 1000;
falsepositive = zeros(numpermutation,1);
clear data1 data2 fp1 fp2 fp3 fp4;
for i = 1:numpermutation
    disp(i)
    
    %create random data
    for k = 1:10
        data1{k} = conv2(randn(100,10),gausswin(5)*gausswin(5)');
        data2{k} = conv2(randn(100,10),gausswin(5)*gausswin(5)');
    end
    
    stat1 = pl_permtest(data1,'statistic','mean');
    fp1(i) = nnz(stat1.FWER.criticalmap(:))>0;
    
    stat2 = pl_permtestcluster(data1,'statistic','mean');
    fp2(i) = nnz(stat2.criticalmap(:))>0;
    
    stat3 = pl_permtest2(data1,data2,'statistic','mean');
    fp3(i) = nnz(stat3.FWER.criticalmap(:))>0;

    stat4 = pl_permtestcluster2(data1,data2,'statistic','mean');
    fp4(i) = nnz(stat4.criticalmap(:))>0;
end
nnz(fp1)/numpermutation
nnz(fp2)/numpermutation
nnz(fp3)/numpermutation
nnz(fp4)/numpermutation