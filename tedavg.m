load('Displacement1_7mm.mat')

x = linspace(1,1000,1000);
%OCT_forw=20*log10(abs(fft(OCT_fringe_forward,1000,2)));

%OCT_back=20*log10(abs(fft(OCT_fringe_backward,1000,2)));

OCT_forw=20*log10(abs(fftshift(fft(OCT_fringe_forward,1000,2))));

OCT_back=20*log10(abs(fftshift(fft(OCT_fringe_backward,1000,2))));
%gauss_envelope = gauss_distribution(x,size(OCT_forw,2)/2,size(OCT_forw,2)/4);
oct_for_cut=OCT_forw(:,500:end);
oct_bac_cut=OCT_back(:,500:end);

fringeFor = OCT_fringe_forward(1,:);

fringeBack = OCT_fringe_backward(1,:);
avgpos=0;
dz=-1;
resvec=zeros(1,500);

OCT_fim=rot90(OCT_forw,-1);
OCT_bim=rot90(OCT_back,-1);
for i=size(OCT_fim,2):-1:1
    %peakFor=OCT_forw(i,:);
    peakFor=OCT_fim(:,i);
    %[mx,pos]=max(peakFor(1:495));
    [mx,pos]=max(peakFor);
    peakBack=OCT_bim(:,i);
    [mxb,posb]=max(peakBack);
    z=abs((pos+posb)/2);
    if avgpos==0
        avgpos=z;
    else
        avgpos = [avgpos z];
    end
    resvec(1,i)=((pos-posb)/2);
%     if dz==-1
%         dz=((pos-posb)/2);
%     else
%         dz = [dz((pos-posb)/2)];
%     end
end
%figure


%resvec=fliplr(resvec);

path=cumtrapz(fliplr(resvec));
temp2=fftshift(path);
temp2(1,251:end)=127;
pathn=temp2./max(temp2);


