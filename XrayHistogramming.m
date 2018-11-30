clear all;
addpath('bin');

patient_folder = dir('DOI_CLEAN');
p1 = [4.6098,0.4234,0.8597];
p2 = [7.0903 0.3686 0.6471];
for j = 1:3
    tic;
    patient_id = patient_folder(j).name;
    if patient_id(1) == 'L'
        CT_path = strcat('DOI_CLEAN/',patient_id,'/CT/');
        Xray_path = strcat('DOI_CLEAN/',patient_id,'/Xray/');
        load(strcat(CT_path,'Proj.mat')); % proj
        
        load(strcat(Xray_path,'Xray1.mat'));
        xray1 = double(img)/2^14;
        load(strcat(Xray_path,'Xray2.mat'));
        xray2 = double(img)/2^14;
        tmp1 = proj(:,:,1)/10;
        tmp2 = proj(:,:,2)/15;
        
%         p1 = Extraction_Fit(tmp1(:),xray1(:),[5,0.4,0.8]);
%         p2 = Extraction_Fit(tmp2(:),xray2(:),[5,0.4,1]);
        proj_processed = zeros(size(proj),'single');
        
        proj_processed(:,:,1) = 1./(1+exp(-(p1(1).*(tmp1-p1(2))))).*p1(3);
        proj_processed(:,:,1) = imsharpen(proj_processed(:,:,1),'Radius',10,'Amount',3);
        
        proj_processed(:,:,2) = 1./(1+exp(-(p2(1).*(tmp2-p2(2))))).*p2(3);
        proj_processed(:,:,2) = imsharpen(proj_processed(:,:,2),'Radius',10,'Amount',3);
        
        figure(1);
 %       subplot(1,3,1); imagesc(max(tmp1(:,:,1)',0)); axis off; axis equal; colormap gray; colorbar; title('CT projection')
        subplot(1,2,1); imagesc(max(proj_processed(:,:,1)',0),[0 1]); axis off; axis equal; colormap gray; colorbar; title('Transferred')
        subplot(1,2,2); imagesc(max(xray1,0),[0 1]); axis off; axis equal; colormap gray; colorbar; title('X-ray')
        
        figure(2);
        subplot(1,3,1); imagesc(max(tmp2(:,:,1)',0)); axis off; axis equal; colormap gray; colorbar; title('CT projection')
        subplot(1,3,2); imagesc(max(proj_processed(:,:,2)',0),[0 1]); axis off; axis equal; colormap gray; colorbar; title('Transferred')
        subplot(1,3,3); imagesc(max(xray2,0),[0 1]); axis off; axis equal; colormap gray; colorbar; title('X-ray')
        saveas(1,'proj_hist.png');
    toc;    
    end
end

%% test;
% p1
% p2
% x = 0:0.01:1;
% y = 1./(1+exp(-(p1(1).*(x-p1(2))))).*p1(3);
% figure(91); plot(x,y); title('Transfer');



