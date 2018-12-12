%% Final Project
% Visualization of code to simulate a continuous-space Hidden Markov Model
%  (HMM) [also known as a linear dynamical system] and fit model parameters 
%  using the Expectation-Maximization (EM) algorithm. This is also known as 
%  a Kalman filter. I use the algorithm for an eye-tracking problem. 
%  We have data tracking the eye movements of mice (pupil diameter, along 
%  with 2 dimensions for eye position [altitude and azimuth]). I fit the 
%  model to get a smoothed estimate of those variables. There is no ground
%  truth dataset in this case for comparison.

%% Simulate Hidden Markov Model, Fit with EM
N = 5000;d = 3;
[x,z,Gamma,Sigma,A,C,mu0,V0] = SimulateKalmanFilterData(N,d);

figure;
plot(z(1,1:100),'c');hold on;plot(z(2,1:100),'m');
xlabel('Time (AU)');ylabel('Magnitude');
title('HMM Hidden States');legend('Dim 1','Dim 2');

% x is the "observed" data, fit the model using EM, find the sequence of
%  hidden states z-hat
[Ahat,Chat,Gammahat,Sigmahat,zhat,mu0hat,V0hat] = KalmanFilterEM(x);

% check dimensions against each other
fprintf('Model Predictions versus True Values\n\n');
figure;
plot(z(1,1:100));hold on;plot(zhat(1,1:100));
title(sprintf('Estimated HMM State Space, Dimension-%d',1));
xlabel('Time (AU)');ylabel('Magnitude');legend('True','Est');

fprintf('State-Space Prediction R-Squared\n');
for ii=1:d
    [r,~] = corrcoef(z(ii,:),zhat(ii,:));
    fprintf('Dimension-%d: %3.2f\n',ii,r(1,2).^2);
end

fprintf('Autoregressive Transition Matrix A\n');
fprintf('True Matrix\n');
disp(A);
fprintf('Estimated Matrix\n');
disp(Ahat);

fprintf('Observation Covariance Matrix Sigma\n');
fprintf('True Matrix\n');
disp(Sigma);
fprintf('Estimated Matrix\n');
disp(Sigmahat);

fprintf('Hidden-To-Observed Transformation Matrix C\n');
fprintf('True Matrix\n');
disp(C);
fprintf('Estimated Matrix\n');
disp(Chat);

%% Simulate Projectile Trajectory, Fit HMM with EM

% use equations of motion to simulate projectile
t = 0:1000;N = length(t);dt = 0.01;z = zeros(N,2);
vel = zeros(N,2);vel(1,:) = [68,68];a = [0,-9.8];
for ii=2:N
    ahat = a-0.005*(vel(ii-1,:).^2).*sign(vel(ii-1,:));
    for jj=1:2
        z(ii,jj) = z(ii-1,jj)+dt*vel(ii-1,jj);
        vel(ii,jj) = vel(ii-1,jj)+dt*ahat(jj);
    end
end
x = z;Sigma = [20,0.05;0.05,20]; % add noise for measurements
for ii=2:N
    x(ii,:) = z(ii,:)+mvnrnd([0,0],Sigma);
end

figure;subplot(2,1,1);plot(z(:,1),z(:,2));xlabel('Distance (AU)');
ylabel('Height (AU)');title('True Projectile Trajectory (State-Space)');
axis([0 max(x(:,1)) 0 max(x(:,2))+10]);
subplot(2,1,2);plot(x(:,1),x(:,2));xlabel('Distance (AU)');
ylabel('Height (AU)');title('Observed Projectile Trajectory (Noisy Measurements)');
axis([0 max(x(:,1)) 0 max(x(:,2))+10]);

% fit HMM
[Ahat,Chat,Gammahat,Sigmahat,zhat,mu0hat,V0hat] = KalmanFilterEM(x');

estTraj = Chat*zhat;estTraj = estTraj';
figure;plot(z(:,1),z(:,2));hold on;plot(estTraj(:,1),estTraj(:,2));
axis([0 max(x(:,1)) 0 max(x(:,2))+10]);
xlabel('Distance (AU)');ylabel('Height (AU)');
title('True and Estimated Projectile Trajectories');
legend('True','Est');

fprintf('Residual Variances\n');
dimension = {'Distance (X)','Height (Y)'};
for ii=1:2
    fprintf('Dimension: %s\n',dimension{ii});
    fprintf('Original Observations: %3.2f\n',var(x(:,ii)-z(:,ii)));
    fprintf('Model Estimates: %3.2f\n\n',var(estTraj(:,ii)-z(:,ii)));
end

% v = VideoWriter('Projectile_Trajectory.avi');
% v.FrameRate = 60;
% open(v);
% for ii=2:1000
%     plot(x(ii,1),x(ii,2),'.b');xlabel('Distance (AU)');
%     ylabel('Height (AU)');title('Observed and Estimated Projectile Trajectories');
%     axis([0 max(x(:,1)) 0 max(x(:,2))+10]);
%     hold on;
%     plot(estTraj(ii-1:ii,1),estTraj(ii-1:ii,2),'-r');F = getframe;
%     writeVideo(v,F.cdata);
% end;
% close(v);

%% Simulation Results
% We see from the simulations that the EM algorithm accurately 
%  recovers model parameters and the hidden state-space variables for 
%  these test cases. In the first case, an HMM is simulated exactly
%  and so we can check the results against the true parameters. With 5000
%  observations, the algorithm captures about 40-60% of the explained variance 
%  (R^2) for the hidden state-space variables, even in the presence of 
%  considerable observation noise. It does an impressive job recovering the
%  parameters used to simulate the data as well. It could perhaps be aided
%  by some kind of regularization, e.g. C is set to the identity matrix but
%  the algorithm estimates C as having non-zero off-diagonal elements.
%  These non-zero elements are quite small though, with the largest rarely
%  exceeding 0.1. With regards to the second simulation, the HMM does an
%  excellent job recovering the trajectory of a simulated projectile. The
%  trajectory was simulated using the equations of motion, and observations 
%  were created by adding correlated Gaussian noise to the parabola. 
%  The HMM estimate of the trajectory is substantially better than using 
%  the observations alone, with residual variance decreasing from about 20 
%  to 0.7, a 30-fold improvement with 1000 observations.

%% Eye-Tracking Data
filename = 'RetinoCall_20180922-4016323';
figure;
load('ExampleIm.mat');
imshow(im);title('Example Eye-Tracker Frame');

load(strcat(filename,'-Init.mat'),'minX','maxX','minY','maxY');
load(strcat(filename,'-MLP.mat'),'pupilDiameter','pupilTranslation','pupilRotation','Fs');
imCoords = pupilTranslation(1,:)+pupilRotation(1,:);
pupDiam = pupilDiameter(1);

figure;
tmp = mean(im(minY:maxY,minX:maxX,:),3);
imshow(tmp);caxis([40 60]);
title('Example Eye-Tracker Frame, Zoom');

for ii=1:size(tmp,1)
    for jj=1:size(tmp,2)
        R = sqrt((imCoords(2)-ii).^2+(imCoords(1)-jj).^2);
        if abs(R-pupDiam/2)<1
            tmp(ii,jj) = 60;
        end
    end
end
figure;
imshow(tmp);caxis([40 60]);
title('Example Eye-Tracker Frame, with Pupil Estimate');

tmpN = 1000;
time = linspace(0,tmpN/Fs,tmpN);
figure;plot(time,pupilDiameter(1:tmpN));
xlabel('Time (seconds)');ylabel('Pupil Diameter (pixels)');
title('Estimated Pupil Diameter');

figure;plot(time,pupilRotation(1:tmpN,1));
hold on;plot(time,pupilRotation(1:tmpN,2));
xlabel('Time (seconds)');
ylabel('Pupil Position (pixels)');
title('Estimated Pupil Position');
legend('Azimuth','Altitude');


% run EM algorithm with about 24,000 frames
data = [pupilDiameter';pupilRotation'];
data = data(:,1:5000); % use subset of data for these purposes
[Ahat,Chat,Gammahat,Sigmahat,zhat,mu0hat,V0hat] = KalmanFilterEM(data);

estData = Chat*zhat;

figure;plot(time,pupilDiameter(1:tmpN));
hold on;plot(time,estData(1,1:tmpN));
xlabel('Time (seconds)');ylabel('Pupil Diameter (pixels)');
title('Estimated Pupil Diameter');
legend('True','Est');

figure;plot(time,pupilRotation(1:tmpN,1));
hold on;plot(time,estData(2,1:tmpN));
xlabel('Time (seconds)');
ylabel('Pupil Position (pixels)');
title('Estimated Pupil Azimuth');
legend('True','Est');

figure;plot(time,pupilRotation(1:tmpN,2));
hold on;plot(time,estData(3,1:tmpN));
xlabel('Time (seconds)');
ylabel('Pupil Position (pixels)');
title('Estimated Pupil Altitude');
legend('True','Est');

% origFrames = cell(1000,1);
% for ff=1:1000
%     tmp = frames{ff};
%     imCoords = pupilTranslation(ff,:)+pupilRotation(ff,:);
%     pupDiam = pupilDiameter(ff);
%     imCoords(1) = imCoords(1)+minX;imCoords(2) = imCoords(2)+minY;
%     tmp2 = tmp;
%     for ii=1:size(tmp,1)
%         for jj=1:size(tmp,2)
%             R = sqrt((imCoords(2)-ii).^2+(imCoords(1)-jj).^2);
%             if abs(R-pupDiam/2)<1
%                 tmp2(ii,jj) = 60;
%             end
%         end
%     end
%     origFrames{ff} = tmp2(30:175,60:225);
% end
% 
% kalmanFrames = cell(1000,1);
% % load('RetinoCall_20180922-4016323-Kalman.mat')
% for ff=1:1000
%     tmp = frames{ff};
%     imCoords = pupilTranslation(ff,:)+estData(2:3,ff)';
%     pupDiam = estData(1,ff);
%     imCoords(1) = imCoords(1)+minX;imCoords(2) = imCoords(2)+minY;
%     tmp2 = tmp;
%     for ii=1:size(tmp,1)
%         for jj=1:size(tmp,2)
%             R = sqrt((imCoords(2)-ii).^2+(imCoords(1)-jj).^2);
%             if abs(R-pupDiam/2)<1
%                 tmp2(ii,jj) = 60;
%             end
%         end
%     end
%     kalmanFrames{ff} = tmp2(30:175,60:225);
% end
% 
% v = VideoWriter('PostHMMTracking.avi');
% v.FrameRate = 50;
% open(v);
% for ii=1:1000
%     figure(1);
%     subplot(1,2,1);imshow(origFrames{ii});caxis([40 60]);subplot(1,2,2);imshow(kalmanFrames{ii});caxis([40 60]);
%     F = getframe(figure(1));
%     writeVideo(v,F.cdata);
% end
% close(v);


%% Eye-Tracker Results
% The Kalman filter provides smoothed estimates of the three eye-tracking
% variables (pupil diameter, pupil azimuth, and pupil altitude). Based on
% the covariance matrices, it is clear that the algorithm "trusts" two of
% the observations. The Sigma matrix (observation covariance) has much smaller
% diagonal values than the Gamma matrix (latent covariance), except for the
% estimate of the altitude. This is clear from the plots as well: the
% diameter and azimuth estimates closely track the observed data, while the
% altitude estimate is more smoothed relative to the observed data. Though 
% not shown, it is clear from viewing the videos that the Kalman estimates
% are better than the original observations. We don't have any ground-truth
% data for the eye-tracking, though, so we can't prove improvement. The
% biggest concerns with tracking the pupil are twofold: 1) the edges of the
% eye sometimes look a bit like a pupil, so the initial estimates of the 
% pupil position occasionally jump. Ideally, the Kalman filter would "notice"
% these jumps and weight the estimate away from a bad observation. However,
% there is also the second concern: 2) the mouse makes periodic saccades, 
% which means the autoregressive process underlying the hidden space probably 
% has two components, one that occurs 99% of the time while the eyes are 
% resting or drifting, and the other that occurs during saccades. The 
% saccades present a problem for any algorithm, because in some sense they 
% look like the noisy jumps that I described above. 