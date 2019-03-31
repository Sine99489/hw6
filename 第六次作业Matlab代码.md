第六次作业Matlab代码

~~~matlab
function im_out=im_filter(im_in,type,M,a)
%空域滤波器
switch type
    %算术均值滤波
    case 'amean'
        h=fspecial('average',[M,M]);
        im_out=imfilter(im_in,h);
    %几何均值滤波
    case 'gmean'
        h=fspecial('average',[M,M]);
        im_out=real(exp(imfilter(log(im_in+0.1),h,'replicate')))-0.1;
    %谐波均值滤波
    case 'hmean'
        h=fspecial('average',[M,M]);
        im_out=1./imfilter(1./(im_in+0.1),h,'replicate')-0.1;
    %逆谐波均值滤波
    case 'chmean'
        Q=a;
        h=ones(M,M);
        im_out=imfilter(im_in.^(Q+1),h,'replicate')./(eps+imfilter(im_in.^Q,h,'replicate'));
    %中值滤波
    case 'median'
        im_out=medfilt2(im_in,[M,M],'symmetric');
    %最大值滤波
    case 'max'
        im_out=ordfilt2(im_in,M*M,ones(M,M),'symmetric');
    %最小值滤波
    case 'min'
        im_out=ordfilt2(im_in,1,ones(M,M),'symmetric');
    %中点滤波
    case 'middle'
        im_out=0.5*(ordfilt2(im_in,M*M,ones(M,M),'symmetric')+ordfilt2(im_in,1,ones(M,M),'symmetric'));
    %修正的阿尔法均值滤波
    case 'alpha'
        d=a;
        im_out=imfilter(im_in,ones(M,M),'symmetric');
        for i=1:d/2
            im_out=im_out-ordfilt2(im_in,i,ones(M,M),'symmetric')-ordfilt2(im_in,M*M+1-i,ones(M,M),'symmetric');
        end
        im_out=im_out/(M*M-d);
    %自适应局部降噪滤波
    case 'adplf'
        nv=a;
        h=fspecial('average',[M,M]);
        im_lm=imfilter(im_in,h,'symmetric');
        im_lv=imfilter((im_in-im_lm).^2,h,'symmetric');
        im_lv(im_lv>nv)=nv;
        im_out=im_in-nv*(im_in-im_lm)./im_lv;
    %自适应中值滤波
    case 'adpmed'
        alP=false(size(im_in));
        im_out=zeros(size(im_in));
        for i=3:2:M
            zmin=ordfilt2(im_in,1,ones(i,i),'symmetric');
            zmax=ordfilt2(im_in,i*i,ones(i,i),'symmetric');
            zmed=medfilt2(im_in,[i,i],'symmetric');
            pB=(zmed>zmin)&(zmax>zmed)&(~alP);
            zB=(im_in>zmin)&(zmax>im_in);
            opzxy=pB&zB;
            opzmed=pB&(~zB);
            im_out(opzxy)=im_in(opzxy);
            im_out(opzmed)=zmed(opzmed);
            alP=alP|pB;
            if all(alP(:))
                break;
            end
        end
    otherwise
        error('Unknow type!');
end
~~~

~~~matlab
%p6_1.m
%在测试图像上产生高斯噪声lena图-需能指定均值和方差
%并用多种滤波器恢复图像
lena=im2double(imread('lena.bmp'));
[M,N]=size(lena);
subplot(2,4,1)
imshow(lena);
title('原图像');

%用高斯噪声污染图像
m=0;    %设定噪声均值
sd=0.1; %设定标准差
imn=m+sd*randn(M,N);
lena_n=lena+imn;
lena_n(lena_n<0)=0;
lena_n(lena_n>1)=1;
subplot(2,4,2)
imshow(lena_n);
title('加入高斯噪声后的图像');

%用多种滤波器恢复图像
L=5;  %设定空域滤波模板大小
%算术均值滤波
lena_am=im_filter(lena_n,'amean',L);
subplot(2,4,3)
imshow(lena_am);
title('使用算术均值滤波复原图像');
%几何均值滤波
lena_gm=im_filter(lena_n,'gmean',L);
subplot(2,4,4)
imshow(lena_gm);
title('使用几何均值滤波复原图像');
%谐波均值滤波
lena_hm=im_filter(lena_n,'hmean',L);
subplot(2,4,5)
imshow(lena_hm);
title('使用谐波均值滤波复原图像');
%中点滤波
lena_mid=im_filter(lena_n,'middle',L);
subplot(2,4,6)
imshow(lena_mid);
title('使用中点滤波复原图像');
%修正的阿尔法均值滤波
lena_al=im_filter(lena_n,'alpha',L,L);
subplot(2,4,7)
imshow(lena_al);
title('使用修正的阿尔法均值滤波复原图像');
%自适应局部降噪滤波
lena_lf=im_filter(lena_n,'adplf',L,sd^2);
subplot(2,4,8)
imshow(lena_lf);
title('使用局部降噪滤波复原图像');
~~~

~~~matlab
%p6_2.m
%在测试图像lena图加入椒盐噪声(椒和盐噪声密度均是0.1)
%用学过的滤波器恢复图像
%在使用反谐波分析Q大于0和小于0的作用
lena=im2double(imread('lena.bmp'));
[M,N]=size(lena);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
subplot(1,5,1)
imshow(lena);
title('原图像');
%使用椒盐噪声污染图像
r=rand(M,N);
lena_sp=lena;
lena_sp(r<=0.1)=0;
lena_sp(r>=0.9)=1;
subplot(1,5,2)
imshow(lena_sp);
title('加入椒盐噪声后的图像');

%用多种滤波器恢复图像
L=5;  %设定空域滤波模板大小
%中值滤波
lena_med=im_filter(lena_sp,'median',L);
subplot(1,5,3)
imshow(lena_med);
title('使用中值滤波复原图像');
%修正的阿尔法均值滤波
lena_al=im_filter(lena_sp,'alpha',L,2*L);
subplot(1,5,4)
imshow(lena_al);
title('使用修正的阿尔法均值滤波复原图像');
%自适应中值滤波
lena_adm=im_filter(lena_sp,'adpmed',7);
subplot(1,5,5)
imshow(lena_adm);
title('使用自适应中值滤波复原图像');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
subplot(1,4,1)
imshow(lena);
title('原图像');
%仅使用胡椒噪声污染图像
r=rand(M,N);
lena_p=lena;
lena_p(r<=0.1)=0;
subplot(1,4,2)
imshow(lena_p);
title('加入胡椒噪声后的图像');

%用多种滤波器恢复图像
L=5;  %设定空域滤波模板大小
%逆谐波均值滤波(Q=1.5)
lena_chm1=im_filter(lena_p,'chmean',L,1.5);
subplot(1,4,3)
imshow(lena_chm1);
title('使用逆谐波均值滤波(Q=1.5)复原图像');
%最大值滤波
lena_max=im_filter(lena_p,'max',L);
subplot(1,4,4)
imshow(lena_max);
title('使用最大值滤波复原图像');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
subplot(1,5,1)
imshow(lena);
title('原图像');
%仅使用盐粒噪声污染图像
r=rand(M,N);
lena_s=lena;
lena_s(r<=0.1)=1;
subplot(1,5,2)
imshow(lena_s);
title('加入盐粒噪声后的图像');

%用多种滤波器恢复图像
L=5;  %设定空域滤波模板大小
%谐波均值滤波
lena_hm=im_filter(lena_s,'hmean',L);
subplot(1,5,3)
imshow(lena_hm);
title('使用谐波均值滤波复原图像');
%逆谐波均值滤波(Q=-1.5)
lena_chm2=im_filter(lena_s,'chmean',L,-1.5);
subplot(1,5,4)
imshow(lena_chm2);
title('使用逆谐波均值滤波(Q=-1.5)复原图像');
%最小值滤波
lena_min=im_filter(lena_s,'min',L);
subplot(1,5,5)
imshow(lena_min);
title('使用最小值滤波复原图像');
~~~

~~~matlab
function  [ff,snr]=wiener(f,G,H,K)
%维纳滤波器
[M,N]=size(f);
FF=(conj(H)./(H.*conj(H)+K)).*G;
ff=real(ifft2(FF));
ff=ff(1:M,1:N);
snr=sum(ff(:).^2)/sum((f(:)-ff(:)).^2);
~~~

~~~matlab
function [ff,gama,r1]=clsf(G,H,m,sd)
%约束最小二乘方滤波器
[M,N]=size(G);
M=M/2;
N=N/2;
n=M*N*(m^2+sd^2);
p=[0,-1,0;-1,4,-1;0,-1,0];
P=fft2(p,2*M,2*N);
gama=0.2;
a=0.01;
dg=0.001;
while 1
    FF=conj(H)./(conj(H).*H+gama*conj(P).*P).*G;
    R=G-H.*FF;
    r=real(ifft2(R));
    r=r(1:M,1:N);
    r1=sum(r(:).^2);
    if abs(r1-n)<=a
        break;
    elseif r1<n-a
        gama=gama+dg;
    else
        gama=gama-dg;
    end
end
ff=real(ifft2(FF));
ff=ff(1:M,1:N);
~~~

~~~matlab
%p6_3.m
%推导维纳滤波器并实现下边要求；
%(a) 实现模糊滤波器如方程Eq. (5.6-11).
%(b) 模糊lena图像：45度方向，T=1；
%(c) 再模糊的lena图像中增加高斯噪声，均值= 0 ，方差=10 pixels 以产生模糊图像；
%(d) 分别利用方程 Eq. (5.8-6)和(5.9-4)，恢复图像；并分析算法的优缺点.
lena=im2double(imread('lena.bmp'));
[M,N]=size(lena);
figure(1);
subplot(1,3,1)
imshow(lena);
title('原图像');

%均匀线性运动模糊
F=fft2(lena,2*M,2*N);
%F=fftshift(F);
H=zeros(2*M,2*N);
T=1;
a=0.1;
b=0.1;
for u=1:2*M
    for v=1:2*N
        if ((u-M)*a+(v-N)*b)~=0
            H(u,v)=T/(pi*((u-M)*a+(v-N)*b))*sin(pi*((u-M)*a+(v-N)*b))*exp(-1i*pi*((u-M)*a+(v-N)*b));
        else
            H(u,v)=T*exp(-1i*pi*((u-M)*a+(v-N)*b));
        end
    end
end
H=ifftshift(H);
G=H.*F;
%G=ifftshift(G);
g=ifft2(G);
lena_m=real(g(1:M,1:N));
subplot(1,3,2)
imshow(lena_m);
title('运动模糊后的图像');

%加入高斯噪声
m=0;    %设定噪声均值
sd=sqrt(10)/255; %设定标准差
imn=m+sd*randn(2*M,2*N);
GG=G+fft2(imn);
gg=ifft2(GG);
lena_n=real(gg(1:M,1:N));
% lena_n(lena_n<0)=0;
% lena_n(lena_n>1)=1;
subplot(1,3,3)
imshow(lena_n);
title('加入高斯噪声后的图像');

%维纳滤波
for i=1:200
    k(i)=i*0.0001;
    [~,snr(i)]=wiener(lena,GG,H,k(i));
end
figure(2)
plot(k,snr)
xlabel('K')
ylabel('SNR')
[~,j]=max(snr);
K=k(j);
[lena_w,snr(i)]=wiener(lena,GG,H,K);
figure(3)
imshow(lena_w)
title('维纳滤波处理后的图像');

%约束最小二乘滤波
[lena_c,gama,r1]=clsf(GG,H,m,sd);
figure(4)
imshow(lena_c)
title('约束最小二乘滤波处理后的图像');
~~~

