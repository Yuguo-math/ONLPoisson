% Script reproducing the results in Table 1 of the LNLA2009 paper [1].
%
% References:
% [1] Mäkitalo, M., and A. Foi, "On the inversion of the Anscombe transformation in low-count Poisson image denoising", Proc. Int. Workshop on Local and Non-Local Approx. in Image Process., LNLA 2009, Tuusula, Finland, pp. 26-32, August 2009.
% [2] Mäkitalo, M., and A. Foi, "Optimal inversion of the Anscombe transformation in low-count Poisson image denoising", submitted, October 2009.
% [3] Zhang B., J.M. Fadili, and J-L. Starck, "Wavelets, ridgelets, and curvelets for Poisson noise removal", IEEE Trans. Image Process., vol. 17, no. 7, pp. 1093-1108, July 2008.
%
%  Alessandro Foi and Markku Mäkitalo - Tampere University of Technology - 2009
% -------------------------------------------------------------------------------

clear all
close all

% mex oxf_poisson004.c
load images.mat   % Load the five images used for the experiments reported in [3] and [1,2]. These images are kindly provided by the authors of [3].
% The image 'Galaxy' is copyright of Commissariat ?l'Énergie Atomique (CEA) / Jean-Luc Starck, www.cea.fr, included here with permission.
% The image 'Cells' is originally from the ImageJ package http://rsb.info.nih.gov/ij (see http://rsb.info.nih.gov/ij/disclaimer.html).
y{1}=spots;
y{2}=galaxy;
y{3}=Ridges;
y{4}=Barbara;
y{5}=cells;
image_name{1}='Spots';
image_name{2}='Galaxy';
image_name{3}='Ridges';
image_name{4}='Barbara';
image_name{5}='Cells';
m1=1;m2=10;
s1=1;s2=10;
m0=9;
s0=6
disp('   ')
disp('----------------------------------------------------------------')

data_nmise=[];
data_psnr =[];
NMISE_nlmsum=zeros(1,5);
for iii=1:5
    for jjj=1:5
        if jjj==1
            m0=9
            s0=6
            J=5
            FPR=0.01
            N_MAX=5
            J_B3=5
            FPR_B3=0.01
            N_MAX_B3=5
            wid=5
            si=1
            nlm0=9
            nls0=6
            nlwid=5
            nlsi=1
            nlh=0.15
        elseif jjj==2
            m0=7
            s0=2
            J=5
            FPR=0.0001
            N_MAX=5
            J_B3=3
            FPR_B3=0.0001
            N_MAX_B3=10
            wid=5
            si=1
            nlm0=6
            nls0=1
            nlwid=5
            nlsi=1
            nlh=0.2
        elseif jjj==3
            m0=4
            s0=9
            J=5
            FPR=0.001
            N_MAX=5
            J_B3=3
            FPR_B3=0.00001
            N_MAX_B3=10
            wid=7
            si=2
            nlm0=4
            nls0=10
            nlwid=7
            nlsi=2
            nlh=0.2
        elseif  jjj==4
            m0=7
            s0=10
            J=4
            FPR=0.001
            N_MAX=5
            J_B3=5
            FPR_B3=0.001
            N_MAX_B3=5
            wid=1
            si=0
            nlm0=7
            nls0=10
            nlwid=1
            nlsi=0
            nlh=0.1
        else
            m0=5
            s0=8
            J=5
            FPR=0.0001
            N_MAX=5
            J_B3=5
            FPR_B3=0.001
            N_MAX_B3=10
            wid=3
            si=0.6
            nlm0=3
            nls0=6
            nlwid=3
            nlsi=1
            nlh=0.2
        end
        
        %      nlw=nlm_poisson004(z,nlm0,nls0,nlh);
        %jjj=3
        
        NMSIE5s0=[m0 0 s1:s2];
        
        
        %jjj=1
        %randn('seed',0);  rand('seed',0);   %% fixes pseudo-random noise
        jj=0;
        NMISEs0=[m0 jjj];
        
        
        z=poissrnd(y{jjj});
        imshow(z,[])
        %     [y_hat, PSNR_y_hat, NMISE_y_hat] = Poisson_denoising_Anscombe_exact_unbiased_inverse(z, y{jjj});   %%  denoise
        jj=jj+1;%% generates Poisson-distributed observations (noisy data)
        %     Idenoised1 = msvst79Denoise (z, FPR, J, 1, 0, 0, 1, 0, N_MAX);
        %     Idenoised2 = msvstB3Denoise (z, 0.001,5, 1, 0, 0, 1, 0, 10, 0, 0);
        
        
        %     timestart = cputime;
        %     w=oxf_poisson004(z,m0,s0);
        %
        %     timefini= cputime;
        %     %NMISE_o_hat=fnmise(y{jjj},w)
        %     if  wid>1
        %         h = fspecial('gaussian', wid,si);
        %         wh{jjj}=imfilter(w,h,'symmetric');
        %     else
        %         wh{jjj}=w;
        %     end
        %wh{jjj}=oxf_poisson004(z,m0,s0);
        %     matchtime = timefini-timestart
        
        %nlw=nlm_poisson004(z,nlm0,nls0,O.1);
        nlw=nlm_poisson004(z,nlm0,nls0,nlh);
        %figure;imshow(z);
        %figure;imshow(nlw);
        if nlwid>1
            h = fspecial('gaussian', nlwid,nlsi);
            nlwh{jjj}=imfilter(nlw,h,'symmetric');
        else
            nlwh{jjj}=nlw;
        end
        %figure;imshow( nlwh{jjj});
        %figure;imshow(wh{jjj});
        %figure;image(wh{jjj});
        %     NMISE_bo79d=fnmise(y{jjj},Idenoised1)
        %     NMISE_bob3d=fnmise(y{jjj},Idenoised2)
        %
        %     NMISE_opf=fnmise(y{jjj},wh{jjj})
        NMISE_nlm(jjj)=fnmise(y{jjj},nlwh{jjj})
        %     NMISE_o_hat=fnmise(y{jjj},wh{jjj})
        %     NMISEs0(2+jj)=NMISE_o_hat;
        
        %data_psnr =[];
        %     data_nmise_0 = [jjj, NMISE_nlm, NMISE_opf,NMISE_y_hat,NMISE_bo79d,NMISE_bob3d]
        %     data_nmise   = cat(1,data_nmise,data_nmise_0)
        %data_psnr_0  = [jjj, psnr(y{jjj},wh{jjj}),psnr(y{jjj},y_hat),psnr(y{jjj},wh{jjj}),psnr(y{jjj},wh{jjj})]
        
        %     disp('----------------------------------------------------------------')
        %     disp('   ')
        %
        %     original = y{jjj};
        %     max_or   = max(original(:));
        %     min_or   = min(original(:));
        %     original = standard_image(original,max_or,min_or);
        %     original_name=[image_name{jjj},'_original' ]
        %     imwrite(uint8(original),[original_name,'.bmp'],'bmp');
        %     imwrite(uint8(original),[original_name,'.tiff'],'tiff');
        %
        %     noisy = standard_image(z,max_or,min_or);
        %     noisy_name=[image_name{jjj},'_noisy' ]
        %     imwrite(uint8(noisy),[noisy_name,'.bmp'],'bmp');
        %     imwrite(uint8(noisy),[noisy_name,'.tiff'],'tiff');
        %
        %
        %     outopfimage= [image_name{jjj},'_owf' ]
        %     wh1=wh{jjj};
        %     wh1=standard_image(wh1,max_or,min_or);
        %     imwrite(uint8(wh1),[outopfimage,'.bmp'],'bmp');
        %     imwrite(uint8(wh1),[outopfimage,'.tiff'],'tiff');
        %
        %     outnlmimage= [image_name{jjj},'_nlm' ]
        %     nlwh1=nlwh{jjj};
        %     nlwh1=standard_image(nlwh1,max_or,min_or);
        %     imwrite(uint8(nlwh1),[outnlmimage,'.bmp'],'bmp');
        %     imwrite(uint8(nlwh1),[outnlmimage,'.tiff'],'tiff');
        %
        %     outbm3d = [image_name{jjj},'_bm3d' ]
        %     wh1=standard_image(y_hat,max_or,min_or);
        %     imwrite(uint8(wh1),[outbm3d,'.bmp'],'bmp');
        %     imwrite(uint8(wh1),[outbm3d,'.tiff'],'tiff');
        %
        %     outbm3d = [image_name{jjj},'_bo79d' ]
        %     wh1=standard_image(Idenoised1,max_or,min_or);
        %     imwrite(uint8(wh1),[outbm3d,'.bmp'],'bmp');
        %     imwrite(uint8(wh1),[outbm3d,'.tiff'],'tiff');
        %
        %     outbm3d = [image_name{jjj},'_bob3d' ]
        %     wh1=standard_image(Idenoised2,max_or,min_or);
        %     imwrite(uint8(wh1),[outbm3d,'.bmp'],'bmp');
        %     imwrite(uint8(wh1),[outbm3d,'.tiff'],'tiff');
        
    end
    NMISE_nlmsum = NMISE_nlmsum+NMISE_nlm
end
NMISE_nlmMeans = NMISE_nlmsum./5
