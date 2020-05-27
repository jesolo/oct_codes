fringes_real = table2array(readtable("C:\Users\korzen\Documents\praca\ted_calc\labview\pig_eye + rubber\real.csv"));
asr=fringes_real(200,:);

OCT_real=20*log10(abs((fft(fringes_real,size(fringes_real,2),2))));
imagesc(fftshift(OCT_real',1))
uenvelopf=zeros(size(fringes_real,1),size(fringes_real,2));

sigphasef=zeros(size(fringes_real,1),size(fringes_real,2));

corphasf=zeros(size(fringes_real,1),size(fringes_real,2));

xx=linspace(1,size(corphasf,2),size(corphasf,2));
corr_fring_real=fringes_real;


asre=fringes_real(200,:);

% BPF = zeros(size(fringes_real,2),1)';
% 
% BPF(1:(size(fringes_real,2)/2)-15) = 1;
% BPF((size(fringes_real,2)/2)+15:end) = 1;
% 
% for j=1:size(fringes_real,1)
%     sas = fft(fringes_real(j,:)); %jeden a-scan
%     signal = abs(fftshift(sas));
%     outc = signal.*BPF;
%     signalc=ifft(ifftshift(outc));
% %     signalc = abs(signalc);
% %     signalc=signalc-mean(signalc);
%     re=signalc;
%     re(1:10)=0;
%     re(end-10:end)=0;
%     fringes_real(j,:)=real(re);
%     fringes_real1(j,:)=abs(re);
% end
% 
% for iii=1:size(fringes_real,1)
% 
%     uenvelopf(iii,:) = envelope(fringes_real(iii,:));
%     hh = hilbert(fringes_real(iii,:));
%     sigphasef(iii,:) = (unwrap(angle(hh)));
%     fitphf = fit(xx',sigphasef(iii,:)','poly1');
%     corphasf(iii,:)=xx.*fitphf.p1+fitphf.p2;
%     %corr_fring_forward(iii,:)=real((hilbert(fringes_forward(iii,:).*exp(1i.*corphasf(iii,:)))));
%     %corr_fring_forward(iii,:)=real(real(hilbert(fringes_forward(iii,:))).*exp(1i.*corphasf(iii,:)));
%     corr_fring_real(iii,:)=real(uenvelopf(iii,:).*exp(1i.*corphasf(iii,:)));
%     %corr_fring_real(iii,:)=real((hilbert(fringes_real(iii,:))).*exp(1i.*corphasf(iii,:)));
%     %corr_fring_backward(iii,:)=real(abs(fringes_forward(iii,:)).*exp(1i.*corphasf(iii,:)));
%     %corr_fring_forward(iii,:)=uenvelopf(iii,:).*exp(-1i.*corphasf(iii,:));
%     %corr_fring_forward(iii,:) = abs(hh).*exp(1i.*corphasf(iii,:)); raczej
%     %nie
%     
% end

OCT_real=20*log10(abs((fft(corr_fring_real,size(fringes_real,2),2))));

figure (5)
imagesc(fftshift(OCT_real',1))
imagesc(OCT_real(:,1:190));

as=fringes_real(200,:);
as1=fringes_real(200,:);
iii=1;
g = 77777;
ix =0;
g2=0;
ix2=0;

mdwa=0;
idwa=0;
limval = 19;
varlim = 7;
%vl= 4;
%segmentacja forward
%for iii=size(OCT_forw,1):-1:1  
for iii=1:size(OCT_real,1)
    ccc=[];
    temp=0;
    temp2=0;
    [mxr inr] = max(OCT_real(iii,:));
    temp = OCT_real(iii,:);
    tresh = 100;
    %tresh = mxr-mxr*0.1;
    while(isempty(ccc)==1)
        
        ccc = find(temp>tresh);
        for jj=1:length(ccc)
            dolim = limval - varlim;
            golim = limval + varlim;
            if (dolim<17)
                dolim = 17;
            end 
            if (golim > 170)
                golim = 170;
            end 
            
            if ccc(jj)<dolim || ccc(jj)>golim
                ccc(jj)=0;
                %            inr = ccc(jj);
                %            break;
            end
        end
        ccc = ccc(ccc~=0);
        if length(ccc)>1
            for jj=1:length(ccc)
                temp2(jj) = temp(ccc(jj));
            end
        else
            temp2=temp(ccc);
        end
        [mdwa idwa] = max(temp2);
        inr = ccc(idwa);
        tresh = tresh - 1;
         
    end
    limval = inr;
    [mx2r in2r] = max(OCT_real(iii,300:end-10));
    %    [mx in] = max(OCT_real(iii,10:160));
%    [mx2 in2] = max(OCT_real(iii,1455:end-15));
   %[mx in] = max(OCT_forw(iii,:));
   if g == 77777
       g = mxr;
       ix = inr;
       g2 = mx2r;
       ix2 = in2r;
   else
%        if inr < 10 || mrx < 100 
%            [mxr inr] = max(OCT_real(iii,(ix(iii-1))+10:150); 
%        end
        itx = 20;
%        while (inr<20)
%           [mxr inr] = max(OCT_real(iii,itx:150));
%           itx=itx+1;
%        end
       g = [g mxr]; 
       ix = [ix inr];
       g2 = [g2 mx2r];
       ix2 = [ix2 in2r];
   end
end


II = rot90(OCT_real(:,1:190));
figure (2)
imagesc(II)
hold on 
plot(-ix+199)
hold off
figure (4)
imagesc(II)

xqq = linspace(1,800,4000);
tmp = [ix zeros(1,1800)];
vq1 = interp1(xqq,tmp,xq);
vq1 = circshift(vq1,50);
vq1 =vq1/2.7;
vq1 =vq1+20;
medv = medfilt1(vq1);

figure (7)
plot(sreds)
hold on
plot(medv)
hold off

