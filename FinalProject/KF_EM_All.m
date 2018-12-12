% KF_EM_All.m

files = dir('*-MLP.mat');

data = cell(1,1);chainLen = 1e4;
heldOutData = cell(1,1);

dataCount = 0;heldOutCount = 0;
for ii=1:length(files)
    load(files(ii).name);
    numIter = floor(N/chainLen)-1;
    
    if numIter>0
        for jj=1:numIter
            dataCount = dataCount+1;
            data{dataCount} = [pupilDiameter(1+(jj-1)*chainLen:chainLen*jj)';...
                pupilRotation(1+(jj-1)*chainLen:chainLen*jj,:)'];
        end
        heldOutCount = heldOutCount+1;
        heldOutData{heldOutCount} = [pupilDiameter(end-chainLen/2:end)';...
                pupilRotation(end-chainLen/2:end,:)'];
    end
end

[A,C,Gamma,Sigma,z,mu0,V0] = KalmanFilterEM_MultiChain(data,heldOutData);

for ii=1:5
    option = unidrnd(dataCount);
    tmp = data{option};
    estTmp = C*z{option};
    
    figure;
    for jj=1:3
        subplot(3,1,jj);plot(tmp(jj,:));hold on;plot(estTmp(jj,:));
    end
end
