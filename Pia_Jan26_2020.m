%% Jan 26 2020

% spectrogram(x,window,noverlap,f,fs)
t = 0:0.001:10;
x = chirp(t,1,1,20,'quadratic');
figure; spectrogram(x,128,50,[],1e3);
xlim([0,100])

%%

N = 10000; % 10 seconds of samples
n = 0:N-1;

w0 = pi/20; % 1000*pi / 20 = 50 half periods, 25 full
x0 = sin(w0*n);

w1 = pi/100;
x1 = sin(w1*n);

y = [x0,x1];

figure; spectrogram(y,512,127,1024,1e3);

%%

csplus_table_coeff_disp = [csplus_table,table(csplus_displacements( include_csplus,:))];
csminus_table_coeff_disp = [csminus_table,table(csminus_displacements( include_csminus,:))];

csplus_table_coeff_disp.Properties.VariableNames = {'Index','Mouse','Trial','csplus_max','csplus_coeffs','Displacement'};
csminus_table_coeff_disp.Properties.VariableNames = {'Index','Mouse','Trial','csplus_max','csplus_coeffs','Displacement'};

