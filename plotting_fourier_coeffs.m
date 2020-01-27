load('DisplacementInOdor.mat')
%%

% Convert to trials (row) x time (col)
displacement_ = permute(displacement,[2,1,3]);
displacement_ = reshape( displacement_ , size(displacement_,1), size(displacement_,2)*size(displacement_,3) );
displacement_ = displacement_';
displacement_ = displacement_(:,[1:10:end]);
%%

figure; subplot(2,1,1); plot( squeeze(displacement(1,:,1)) );
subplot(2,1,2); plot( displacement_(1,:) )

figure; subplot(2,1,1); plot( squeeze(displacement(1,:,2)) );
subplot(2,1,2); plot( displacement_(11,:) )

%%

t0=80000;
tend=160000;

T = [t0/10:tend/10];

figure('color','w'); 
subplot(1,2,1);
plot( displacement_(10, T) )
x = displacement_(1,T);
L=length(x);
NFFT= 2^nextpow2( L );
X=fft(x,NFFT);
fs=1000;
Px=X.*conj(X)/(NFFT*L); %Power of each freq components
fVals=fs*(0:NFFT/2-1)/NFFT;
plot(fVals,Px(1:NFFT/2),'b','LineSmoothing','on','LineWidth',1);
title('One Sided Power Spectral Density');
xlabel('Frequency (Hz)')
ylabel('PSD');
xlim([0,10])
subplot(1,2,2); plot(x);


%%

% First ten: CS+
% Second ten: CS-
% Third ten: O3

% Odor: 80 000 to 160 000
% Pre-odor: from 65 000 (80 000 - 15 000) to 80 000
t0 = 8000;
tend = 9000;

figure; for i = 1:30;
    plot( displacement_(i,:) ); patch([t0,t0,tend,tend,t0],[20,-4,-4,20,20],'r','facealpha',0.5)
    title(i);
    pause()
end;

%%


figure('color','w'); 
mycolor = [1,0,0];

%% Smoothing out displacements

smooth_displacement_ = arrayfun( @(x) smooth( displacement_(x,:), 10 ), [1:30], 'uniformoutput', false );
smooth_displacement_ = cell2mat( smooth_displacement_ )';

%% Unsmoothed
trials = [4,5,7,8,10,12,14,15,16]';
tend = [12903,14516,11967,13157,12788,16294,15162,16853,16456]';
%tend = min(tend,12000);

mytable = table( trials, tend );

output = rowfun( @(trial,tend) myfourier( displacement_(trial,:), 8000, tend,'k' ), mytable(:,[1,2]) );
output.Var2([1:5]) = 1;
output.Var2([6:9]) = 2;

figure; boxplot(output.Var1,output.Var2); hold on;
scatter(output.Var2,output.Var1,'jitter','on','jitteramount',0.2)
output.Var3 = trials;a
output.Var4 = tend;

%% Smoothed

mytable = table( trials, tend );

output = rowfun( @(trial,tend) myfourier( smooth_displacement_(trial,:), 8000, tend,'k' ), mytable(:,[1,2]) );
output.Var2([1:5]) = 1;
output.Var2([6:9]) = 2;

figure; boxplot(output.Var1,output.Var2); hold on;
scatter(output.Var2,output.Var1,'jitter','on','jitteramount',0.2)
output.Var3 = trials;
output.Var4 = tend;


%% Starts here
