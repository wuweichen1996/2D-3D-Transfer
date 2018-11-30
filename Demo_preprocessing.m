%This is the demo for the whole data-preprocessing process

clear all;
addpath('bin');
WaterAtt = 0.0269; % 50 kVp, 'SOMATOM Definition Flash'
p1 = [4.6098,0.4234,0.8597];
p2 = [7.0903 0.3686 0.6471]; 

%% CT Image generation:from .dcm to .dat
patient_folder = dir('DOI_CLEAN');
for j = 1:length(patient_folder)
    tic;
    patient_id = patient_folder(j).name;
    if patient_id(1) == 'L'
        CT_path = strcat('DOI_CLEAN/',patient_id,'/CT/');
        File = dir(strcat(CT_path,'*.dcm'));

        info = dicominfo(strcat(CT_path,File(1).name));

        dx = info.PixelSpacing(1);
        dy = info.PixelSpacing(1);
        dz = info.SliceThickness;

        % usually nx, ny have the same number
        nx = info.Rows;
        ny = info.Columns;
        nz = length(File);

        img = zeros(nx, ny, nz);

        MaxSliceLocation = -10000;
        MinSliceLocation = 10000;
        % HU value
        for i = 1:nz
            info = dicominfo(strcat(CT_path,File(i).name));
            tmp = info.SliceLocation;
            if tmp<MinSliceLocation
                MinSliceLocation = tmp;
            end
    
            if tmp>MaxSliceLocation
                MaxSliceLocation = tmp;
            end
    
        %     img(:,:,i) = dicomread(File(i).name)*info.RescaleSlope+info.RescaleIntercept;
        end

        for i = 1:nz
            info = dicominfo(strcat(CT_path,File(i).name));    
            iz = (info.SliceLocation-MinSliceLocation)/info.SliceThickness;
            iz = int16(iz);
            img(:,:,end-iz) = dicomread(strcat(CT_path,File(i).name))*info.RescaleSlope+info.RescaleIntercept;
        end
    
        % return;
        % Attenuation value convert
        
        Attimg = max(single(img)/1000*WaterAtt + WaterAtt,0);

        figure(1);
        subplot(1,2,1); imagesc(img(:,:,round(end/2))); axis off; axis equal; colorbar; title('HU');
        subplot(1,2,2); imagesc(Attimg(:,:,round(end/2))); axis off; axis equal; colorbar; title('mm^{-1}');
        
        CT_image_path = strcat(CT_path,'CT_image.mat');
        save(CT_image_path,'Attimg','dx','dy','dz');
    end
    toc;
    fprintf('%d/%d patients completed .mat transformation\n',j,length(patient_folder));
end



%% CT Image Resampling: from 512*512*N to 256*256*256
for j = 1:length(patient_folder)
    tic;
    patient_id = patient_folder(j).name;
    if patient_id(1) == 'L'
        CT_path = strcat('DOI_CLEAN/',patient_id,'/CT/');
        % atlas configuration
        Anx = 256;
        Any = 256;
        Anz = 256;
        Adx = 2; % mm
        Ady = 2;
        Adz = 2;

        resolution = [Adx, Ady, Adz];
        % Load CT image

        load(strcat(CT_path,'CT_image.mat')); % Attimg, dx, dy, dz

        [nx, ny, nz] = size(Attimg);
        % Resampling and store

        ox = single([0:nx-1] - (nx-1.0)/2.0)*dx;
        oy = single([0:ny-1] - (ny-1.0)/2.0)*dy;
        oz = single([0:nz-1] - (nz-1.0)/2.0)*dz;

        [oxx, oyy, ozz] = meshgrid(oy, ox, oz);

        x = single([0:Anx-1] - (Anx-1.0)/2.0)*Adx;
        y = single([0:Any-1] - (Any-1.0)/2.0)*Ady;
        z = single([0:Anz-1] - (Anz-1.0)/2.0)*Adz;
        [xx, yy, zz] = meshgrid(x, y, z);

        img_att = interp3(oxx, oyy, ozz, single(Attimg), xx, yy, zz);
        img_att(isnan(img_att)) = 0;
        
        CT_resample_path = strcat(CT_path,'CT_image_Resample.mat');
        save(CT_resample_path,'img_att','resolution');

        figure(2); 
        subplot(131); imagesc(img_att(:,:,round(end/2))); axis off; axis equal; colormap gray;
        subplot(132); imagesc(squeeze(img_att(:,round(end/2),:))'); axis off; axis equal; colormap gray;
        subplot(133); imagesc(squeeze(img_att(round(end/2),:,:))'); axis off; axis equal; colormap gray;
        
    end
    toc;
    fprintf('%d/%d patients completed resampling\n',j,length(patient_folder));
end

%% Projection of Resampled CT image

%Parameter setting in general
WaterAtt = 0.0269; % 50 keV, 'SOMATOM Definition Flash'
ParamSetting;

for j = 1:length(patient_folder)
    tic;
    patient_id = patient_folder(j).name;
    if patient_id(1) == 'L'
        CT_path = strcat('DOI_CLEAN/',patient_id,'/CT/');
        load(strcat(CT_path,'CT_image_Resample.mat'));
        Orgimg = img_att;
        proj = CTprojectionC(Orgimg,param);
    
        figure(3); 
        subplot(1,2,1); imagesc(max(proj(:,:,1)',0)); axis off; axis equal; colormap gray; colorbar;
        subplot(1,2,2); imagesc(max(proj(:,:,2)',0)); axis off; axis equal; colormap gray; colorbar;
        proj_path = strcat(CT_path,'Proj.mat');
        save(proj_path,'proj');
    end
    toc;
    fprintf('%d/%d patients completed projection\n',j,length(patient_folder));
end

%% X-ray histogram

for j = 1:length(patient_folder)
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

        figure(4);
 %       subplot(1,3,1); imagesc(max(tmp1(:,:,1)',0)); axis off; axis equal; colormap gray; colorbar; title('CT projection')
        subplot(1,2,1); imagesc(max(proj_processed(:,:,1)',0),[0 1]); axis off; axis equal; colormap gray; colorbar; title('Transferred')
        subplot(1,2,2); imagesc(max(xray1,0),[0 1]); axis off; axis equal; colormap gray; colorbar; title('X-ray')
        saveas(4,strcat(CT_path,'proj_hist1'));
        
        figure(5);
 %       subplot(1,3,1); imagesc(max(tmp2(:,:,1)',0)); axis off; axis equal; colormap gray; colorbar; title('CT projection')
        subplot(1,2,1); imagesc(max(proj_processed(:,:,2)',0),[0 1]); axis off; axis equal; colormap gray; colorbar; title('Transferred')
        subplot(1,2,2); imagesc(max(xray2,0),[0 1]); axis off; axis equal; colormap gray; colorbar; title('X-ray')
        saveas(5,strcat(CT_path,'proj_hist2'));
        
        save(strcat(CT_path,'proj_hist.mat'),'proj_processed');
    end
    toc;
end





