fringes_forward = table2array(readtable("C:\Users\korzen\Documents\praca\ted_calc\labview\pig_eye + rubber\p1for.csv"));
fringes_backward = table2array(readtable("C:\Users\korzen\Documents\praca\ted_calc\labview\pig_eye + rubber\p1bac.csv"));
fringes_forward_u = fringes_forward;
fringes_backward_u =fringes_backward;
OCT_forw_uncor=20*log10(abs((fft(fringes_forward_u,size(fringes_forward,2),2))));

OCT_back_uncor=20*log10(abs((fft(fringes_backward_u,size(fringes_forward,2),2))));

uenvelopf=zeros(size(fringes_forward,1),size(fringes_forward,2));
uenvelopb=zeros(size(fringes_forward,1),size(fringes_forward,2));
sigphasef=zeros(size(fringes_forward,1),size(fringes_forward,2));
sigphaseb=zeros(size(fringes_forward,1),size(fringes_forward,2));
corphasf=zeros(size(fringes_forward,1),size(fringes_forward,2));
corphasb=zeros(size(fringes_forward,1),size(fringes_forward,2));
xx=linspace(1,size(corphasf,2),size(corphasf,2));
corr_fring_forward=fringes_forward;
corr_fring_backward=fringes_backward;
xt = linspace(1,size(fringes_forward,2),size(fringes_forward,2));
% OCT_forw=20*log10(abs((fft(corf,size(fringes_forward,2),2))));
% OCT_back=20*log10(abs((fft(corb,size(fringes_forward,2),2))));
% imagesc(fftshift(OCT_forw_uncor',1))
% figure (5)
% imagesc(rot90(OCT_forw,3))
ascuc=fringes_forward(200,:);
ascbuc=fringes_backward(200,:);
BPF = zeros(1752,1)';
BPF2 = zeros(1752,1)';
BPF(1:862) = 1;
BPF(892:end) = 1;
BPF2(601:end) = 1;
BPF2(end-601:end) = 0;
BP = BPF.*BPF2;
for j=1:size(fringes_forward,1)
    sas = fft(fringes_forward(j,:)); %jeden a-scan
    signal = abs(fftshift(sas));
    outc = signal.*BP;
    signalc=ifft(ifftshift(outc));
%     signalc = abs(signalc);
%     signalc=signalc-mean(signalc);
    re=signalc;
    re(1:10)=0;
    re(end-10:end)=0;
    fringes_forward(j,:)=re;
    sab = fft(fringes_backward(j,:));
    signal = abs(fftshift(sab));
    outc = signal.*BP;
    signalc=ifft(ifftshift(outc));
%     signalc = abs(signalc);
%     signalc=signalc-mean(signalc);
    re=signalc;
    %res = abs(signalc);
    %re=signalc-mean(signalc);
    re(1:10)=0;
    re(end-10:end)=0;
    fringes_backward(j,:) = re; 
end

OCT_forw_u=20*log10(abs((fft(fringes_forward,size(fringes_forward,2),2))));
OCT_back_u=20*log10(abs((fft(fringes_backward,size(fringes_forward,2),2))));
% imagesc(fftshift(OCT_forw_u',1))

for iii=1:size(fringes_forward,1)

    uenvelopf(iii,:) = envelope(fringes_forward(iii,:));
    hh = hilbert(fringes_forward(iii,:));
    sigphasef(iii,:) = (unwrap(angle(hh)));
    fitphf = fit(xx',sigphasef(iii,:)','poly1');
    corphasf(iii,:)=xx.*fitphf.p1+fitphf.p2;
    %corr_fring_forward(iii,:)=real((hilbert(fringes_forward(iii,:).*exp(1i.*corphasf(iii,:)))));
    %corr_fring_forward(iii,:)=real(real(hilbert(fringes_forward(iii,:))).*exp(1i.*corphasf(iii,:)));
    %corr_fring_forward(iii,:)=real(uenvelopf(iii,:).*exp(1i.*corphasf(iii,:)));
    corr_fring_forward(iii,:)=real((hilbert(fringes_forward(iii,:))).*exp(1i.*corphasf(iii,:)));
    corr_fring_backward(iii,:)=real(abs(fringes_forward(iii,:)).*exp(1i.*corphasf(iii,:)));
    %corr_fring_forward(iii,:)=uenvelopf(iii,:).*exp(-1i.*corphasf(iii,:));
    %corr_fring_forward(iii,:) = abs(hh).*exp(1i.*corphasf(iii,:)); raczej
    %nie
    uenvelopb(iii,:) = envelope(fringes_backward(iii,:));
    hhh = hilbert(fringes_backward(iii,:));
    sigphaseb(iii,:) = (unwrap(angle(hhh)));
    fitphb = fit(xx',sigphaseb(iii,:)','poly1');
    corphasb(iii,:)=xx.*fitphb.p1+fitphb.p2;
    %plot(abs(fftshift(fft(real(hilbert(s12).*exp(-1i*phi_corr)),length(s12)*N_pad))))
    %corr_fring_backward(iii,:)=real(real(hilbert(fringes_backward(iii,:))).*exp(1i.*corphasb(iii,:)));
    %corr_fring_backward(iii,:)=real(uenvelopb(iii,:).*exp(1i.*corphasb(iii,:)));
    %corr_fring_backward(iii,:)=real((hilbert(fringes_backward(iii,:).*exp(1i.*corphasb(iii,:)))));
    %corr_fring_backward(iii,:)=uenvelopb(iii,:).*exp(-1i.*corphasb(iii,:));
    corr_fring_backward(iii,:)=real(abs(fringes_backward(iii,:)).*exp(1i.*corphasf(iii,:)));
    %corr_fring_backward(iii,:) = abs(hhh).*exp(1i.*corphasb(iii,:));
%     avgph(iii,:)=xx.*(fitphf.p1+fitphb.p1)/2;
%     
%     avgphfor(iii,:)=xx.*(((fitphf.p1+fitphb.p1)/2)+fitphf.p1)/2;
%     avgphbac(iii,:)=xx.*(((fitphf.p1+fitphb.p1)/2)+fitphb.p1)/2;
%     avg_fringes1(iii,:)=real((hilbert(fringes_backward(iii,:))).*exp(1i.*avgphfor(iii,:)));
%     avg_fringes2(iii,:)=real((hilbert(fringes_forward(iii,:))).*exp(1i.*avgphbac(iii,:)));
end


OCT_forw=20*log10(abs((fft(corr_fring_forward,size(fringes_forward,2),2))));
OCT_back=20*log10(abs((fft(corr_fring_backward,size(fringes_forward,2),2))));


figure (4)
imagesc(fftshift(OCT_forw_uncor',1))
figure (5)
imagesc(fftshift(OCT_forw',1))
% N = size(ascuc,2);
% Fs = 10 * 1000;
% t = linspace(0,N*1/Fs,N);
% Ts = 1;                                                             % Sampling Interval
% Fs = 1/Ts;                                                          % Sampling Frequency
% Fn = Fs/2;                                                          % Nyquist Frequency
% L = size(data,1);                                                   % Length Of ‘data’ Vector
% t = 1:L*Ts;  % Time Vector

% dF = Fs/N;
% f = (-Fs/2:dF:Fs/2-dF)';
% 
% 
% 
% test1 = fringes_forward_c(400,:);
% test2 = fringes_forward(400,:);
% 
% ftes1 = abs(fft(test1));
% ftes2 = abs(fft(test2));
iii=1;
g = 77777;
ix =0;
g2=0;
ix2=0;

%segmentacja forward
%for iii=size(OCT_forw,1):-1:1
for iii=1:size(OCT_forw,1)
   [mx in] = max(OCT_forw(iii,:));
   [mx2 in2] = max(OCT_forw(iii,:));
%    [mx in] = max(OCT_forw(iii,20:150));
%    [mx2 in2] = max(OCT_forw(iii,1455:end-15));
   %[mx in] = max(OCT_forw(iii,:));
   if g == 77777
       g = mx;
       ix = in;
       g2 = mx2;
       ix2 = in2;
   else
       g = [g mx]; 
       ix = [ix in];
       g2 = [g2 mx2];
       ix2 = [ix2 in2];
   end
end
g = 77777;
for iii=1:size(OCT_back,1)
   [mxb, inb] = max(OCT_back(iii,:));
   [mxb2, inb2] = max(OCT_back(iii,:));
%    [mxb, inb] = max(OCT_back(iii,20:150));
%    [mxb2, inb2] = max(OCT_back(iii,1455:end-15));
   %[mx in] = max(OCT_forw(iii,:));
   if g == 77777
       g = mxb;
       ixb = inb;
       g2 = mxb2;
       ixb2 = inb2;
   else
       g = [g mxb]; 
       ixb = [ixb inb];
       g2 = [g2 mxb2];
       ixb2 = [ixb2 inb2];
   end
end
figure (9)
plot(ix);hold on;plot(ixb);hold off;
sred=(ix+ixb)./2;
I = (fftshift(OCT_back',1));

medfor = medfilt1(ix);

medback = medfilt1(ix);

xq = linspace(1,length(ix),length(ix));
medix = medfilt1(ix); 
medixb = medfilt1(ixb); 
mixs = smooth(xq,medix,10,'rloess');
mixbs = smooth(xq,medixb,10,'rloess');
sreds=(mixs+mixbs)./2;
sred=(ix+ixb)./2;
plot(mixs);hold on;plot(mixbs);plot(sreds);hold off;

figure (5)
imagesc(fftshift(OCT_forw(:,876:end)',1))
hold on
plot(-mixs+420)
hold off
figure (6)
imagesc(fftshift(OCT_back(:,876:end)',1))
hold on
plot(-mixbs+420)
hold off