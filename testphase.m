% fringes_forward = table2array(readtable("C:\Users\korzen\Documents\praca\ted_calc\labview\24-11-2019\ruber1for.csv"));
% fringes_backward = table2array(readtable("C:\Users\korzen\Documents\praca\ted_calc\labview\24-11-2019\ruber1bac.csv"));

fringes_forward = table2array(readtable("C:\Users\korzen\Documents\praca\ted_calc\realdata\rubberv1closer2zopd\forward.csv"));
fringes_backward = table2array(readtable("C:\Users\korzen\Documents\praca\ted_calc\realdata\rubberv1closer2zopd\backward.csv"));

% fringes_forward = table2array(readtable("C:\Users\korzen\Documents\praca\ted_calc\labview\airpuffrobber\forrobb.csv"));
% fringes_backward = table2array(readtable("C:\Users\korzen\Documents\praca\ted_calc\labview\airpuffrobber\bacrobb.csv"));


OCT_forw_uncor=20*log10(abs((fft(fringes_forward,size(fringes_forward,2),2))));

OCT_back_uncor=20*log10(abs((fft(fringes_backward,size(fringes_forward,2),2))));
asc=fringes_forward(255,:);
ascb=fringes_backward(255,:);
uenvelopf=zeros(size(fringes_forward,1),size(fringes_forward,2));
uenvelopb=zeros(size(fringes_forward,1),size(fringes_forward,2));
sigphasef=zeros(size(fringes_forward,1),size(fringes_forward,2));
sigphaseb=zeros(size(fringes_forward,1),size(fringes_forward,2));
corphasf=zeros(size(fringes_forward,1),size(fringes_forward,2));
corphasb=zeros(size(fringes_forward,1),size(fringes_forward,2));
xx=linspace(1,size(corphasf,2),size(corphasf,2));
corr_fring_forward=fringes_forward;
corr_fring_backward=fringes_backward;
xq=linspace(0,1658,322);

for iii=1:size(fringes_forward,1)
    fringes_forward(iii,:)=fringes_forward(iii,:);%-3.6e4;
    fringes_backward(iii,:)=fringes_backward(iii,:);%-3.6e4;
    uenvelopf(iii,:) = envelope(fringes_forward(iii,:));
    uenvelopb(iii,:) = envelope(fringes_backward(iii,:));

    hh = hilbert(fringes_forward(iii,:));
    sigphasef(iii,:) = (unwrap(angle(hh)));
    intfu= interp1(xx,sigphasef(iii,:),xq);
    intfu2= interp1(xq,intfu,xx);
    [temp1,edges]=histcounts(gradient(intfu2),130);
    [mx,ind]=max(temp1);
    corphasf(iii,:)=xx.*edges(ind);
    
    corr_fring_forward(iii,:)=real(uenvelopf(iii,:).*exp(1i.*corphasf(iii,:)));
    hhh = hilbert(fringes_backward(iii,:));
    sigphaseb(iii,:) = (unwrap(angle(hhh)));
    intfu= interp1(xx,sigphaseb(iii,:),xq);
    intfu2= interp1(xq,intfu,xx);
    
    [temp2,edges2]=histcounts(gradient(intfu2),200);
    [mx2,ind2]=max(temp2);
    corphasb(iii,:)=xx.*edges2(ind2);
    corr_fring_backward(iii,:)=real(uenvelopb(iii,:).*exp(1i.*corphasb(iii,:)));
    avgph(iii,:)=xx.*(edges(ind)+edges2(ind2))/2;
    avg_fringes1(iii,:)=uenvelopf(iii,:).*exp(-1i.*avgph(iii,:));
    avg_fringes2(iii,:)=uenvelopb(iii,:).*exp(1i.*avgph(iii,:));
end

% for iii=1:size(fringes_forward,1)
%     fringes_forward(iii,:)=fringes_forward(iii,:);%-3.6e4;
%     fringes_backward(iii,:)=fringes_backward(iii,:);%-3.6e4;
%     uenvelopf(iii,:) = envelope(fringes_forward(iii,:));
%     uenvelopb(iii,:) = envelope(fringes_backward(iii,:));
% 
%     hh = hilbert(fringes_forward(iii,:));
%     sigphasef(iii,:) = (unwrap(angle(hh)));
%     fitphf = fit(xx',sigphasef(iii,:)','poly1');
%     corphasf(iii,:)=xx.*fitphf.p1+fitphf.p2;
%     corr_fring_forward(iii,:)=real(uenvelopf(iii,:).*exp(-1i.*corphasf(iii,:)));
%     hhh = hilbert(fringes_backward(iii,:));
%     sigphaseb(iii,:) = (unwrap(angle(hhh)));
%     fitphb = fit(xx',sigphaseb(iii,:)','poly1');
%     corphasb(iii,:)=xx.*fitphb.p1+fitphb.p2;
%     corr_fring_backward(iii,:)=real(uenvelopb(iii,:).*exp(-1i.*corphasb(iii,:)));
%     avgph(iii,:)=xx.*(fitphf.p1+fitphb.p1)/2+(fitphf.p2+fitphb.p2)/2;
%     avg_fringes1(iii,:)=uenvelopf(iii,:).*exp(-1i.*avgph(iii,:));
%     avg_fringes2(iii,:)=uenvelopb(iii,:).*exp(1i.*avgph(iii,:));
% end


i=300;
figure (5)
plot(avgph(i,:));hold on;
plot(corphasf(i,:));plot(corphasb(i,:));
hold off;
%asc=fringes_forward(577,:);
OCT_forw=20*log10(abs((fft(corr_fring_forward,size(fringes_forward,2),2))));


OCT_back=20*log10(abs((fft(corr_fring_backward,size(fringes_backward,2),2))));


figure (6)
imagesc(fftshift(OCT_back',1))

OCT_forw_avg=20*log10(abs((fft(avg_fringes1,size(fringes_forward,2),2))));
    
OCT_back_avg=20*log10(abs((fft(avg_fringes2,size(fringes_backward,2),2))));
figure
imagesc(fftshift(OCT_back_avg',1))
figure
imagesc(fftshift(OCT_forw_avg',1))

i=170;
figure (7)
plot(OCT_forw(i,:));hold on;
plot(OCT_forw_uncor(i,:));
hold off;

figure (8)
plot(corphasf(i,:));hold on;
plot(sigphasef(i,:));
hold off;

