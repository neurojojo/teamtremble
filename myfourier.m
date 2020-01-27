function output = myfourier( x, t0, tend, varargin )


vis =0;

% Set parameters
T = [t0:tend];
x = x(T);
L=length(x);
NFFT= 2^nextpow2( L );

X=fft(x,NFFT);
fs=1000;
Px=X.*conj(X)/(NFFT*L); %Power of each freq components
fVals=fs*(0:NFFT/2-1)/NFFT;

    if vis==1
        figure;
        plot(fVals,Px(1:NFFT/2),'b','LineSmoothing','on','LineWidth',1,'color','k');
        title('One Sided Power Spectral Density');
        xlabel('Frequency (Hz)')
        ylabel('PSD');
        xlim([0,10]);
    end

    output = sum( Px([5:12]) ); % For a 1000 timepoint vector this is 4 Hz to 10 Hz

end