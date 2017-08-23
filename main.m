%% empty everything
clear all;
clc;

%% change current path to where "main.m" is
mpath = strrep(which(mfilename),[mfilename '.m'],'');
cd(mpath);

template1=rgb2gray(imread('template/template1.jpg'));

%% load all file under "database" folder
data_num=dir('database');
cd database;

fileID = fopen('record.txt','w');

for k=3:length(data_num)

    %% read images from database one by one and build a same-size mask
    
    fname=data_num(k).name;   
    fprintf(fileID, '====================');
    fprintf(fileID,'%s',fname);
    fprintf(fileID, '====================\r\n');

    image_origin = imread('361.jpg');
    image_gray = rgb2gray(image_origin);
        
    %% hue and saturation clip 
    image_hsv = rgb2hsv(image_origin);
    image_hue=image_hsv(:,:,1);
    
    image_hue(image_hue<15/360) = 1;
    image_hue(image_hue>310/360)= 1;
    image_hue = floor(image_hue);
    
    image_sat = image_hsv(:,:,2);
    image_mask=image_hue.*im2bw(image_sat,0.05);

    %% use open method to eliminate some debris, then fill up holes                 
    afterFill = imfill(image_mask,8,'holes');
    se = strel('disk',23);
    afterCH = imopen(afterFill,se);
%     figure; imshow(afterCH,[]);
    afterCH= logical(afterCH);
    
    %% use 8-neighborhood connection to label parts
    stats = regionprops(afterCH, 'BoundingBox', 'Centroid','Area','Perimeter');     
    
    
    
    
    
    
    
    figure;
    imshow(image_origin);     title(fname);
    hold on;  
    ssimval_list=zeros(size(stats,2));
    ratio_list=zeros(size(stats,2));
    
    for object = 1:length(stats)          
        
        
        roundness = stats(object).Area/((stats(object).Perimeter)^2)*4*pi;
        if(roundness>1.2 || roundness<0.8)
            continue;
        end
        
        
        
        
        
        %% get masked texture 
        bb = stats(object).BoundingBox;
        bc = stats(object).Centroid;    
        bb = int16(floor(bb));
        
        if bb(1) == 0
            bb(1) = 1;
        end
        if bb(2) == 0
            bb(2) = 1;
        end
        i = bb(1) : bb(1)+bb(3)-1;
        j = bb(2) : bb(2)+bb(4)-1;

        
        image_masked = image_gray(j,i);    
%         figure;
%         imshow(image_masked);
%         
        %%SIFT
        %%resize template to get better features that might match target pic at same size.
        image_masked = imresize(image_masked, size(template1));     
        points_image = detectSURFFeatures(image_masked,'MetricThreshold',800,'NumOctaves',6,'NumScaleLevels',10,'ROI',[0.2*size(image_masked,2),0.2*size(image_masked,1),0.6*size(image_masked,2),0.6*size(image_masked,1)]);            
        template1_ratio = 0;
        template2_ratio = 0;
        ssimval=0;
        
        if points_image.Count ~= 0    
            points_template1 = detectSURFFeatures(template1,'MetricThreshold',800,'NumOctaves',6,'NumScaleLevels',10,'ROI',[0.2*size(template1,2),0.2*size(template1,1),0.6*size(template1,2),0.6*size(template1,1)]);         
            
%              figure;subplot(1,2,1);imshow(image_masked); 
%              hold on;
%              plot(points_image);subplot(1,2,2);imshow(template1); 
%              hold on;
%              plot(points_template1);

            [f1,vpts1] = extractFeatures(image_masked,points_image);
            [f2,vpts2] = extractFeatures(template1,points_template1);
            indexPairs = matchFeatures(f1,f2,'MatchThreshold',15,'MaxRatio',0.98,'Unique',1,'Metric','SAD') ;
            matchedPoints1 = vpts1(indexPairs(:,1));
            matchedPoints2 = vpts2(indexPairs(:,2));
             
%            figure;showMatchedFeatures(image_masked,template1,matchedPoints1,matchedPoints2,'montage');

            if (size(indexPairs,1)>3)
                [tform, inlierMaskPoints, inlierTemplatePoints] = estimateGeometricTransform(matchedPoints1, matchedPoints2, 'similarity','MaxNumTrials',1000);

%                 figure;showMatchedFeatures(image_masked,template1,matchedPoints1,matchedPoints2,'montage');
                
                outputView = imref2d(size(template1));
                Ir = imwarp(image_masked,tform,'OutputView',outputView);
                
                ir_col=round(0.2*size(Ir,1)):round(0.8*size(Ir,1));
                ir_row=round(0.2*size(Ir,2)):round(0.8*size(Ir,2));
                t_col=round(0.2*size(template1,1)):round(0.8*size(template1,1));
                t_row=round(0.2*size(template1,2)):round(0.8*size(template1,2));
                
                
                roi_ir = uint8(Ir(ir_col, ir_row));
                roi_template1 = uint8(template1(t_col, t_row));
                
                
                ssimval = ssim(roi_ir,roi_template1) ;
                ssimval_list(object)=ssimval;
                ratio=size(inlierMaskPoints,1)/size(vpts1,1);
                
                ratio_list(object)=ratio;
%                 figure;imshow(roi_ir);
            end
        end
   
        if (ssimval>0.45 && ratio > 0.01)
%             figure;showMatchedFeatures(image_masked,template1,matchedPoints1,matchedPoints2,'montage');
            targetPolygon = [bb(1), bb(2);...             % top-left
                             bb(1)+bb(3), bb(2);...       % top-right
                             bb(1)+bb(3), bb(2)+bb(4);... % bottom-right
                             bb(1),bb(2)+bb(4);...        % bottom-left
                             bb(1), bb(2)];               % top-left again to close the polygon                        
            line(targetPolygon(:, 1), targetPolygon(:, 2), 'Color', 'b','LineWidth',2);
            fprintf(fileID,'%f\t%f\r\n',ratio,ssimval);
        end         
    end  
    fprintf(fileID,'\r\n\r\n\r\n');
end
fclose(fileID);