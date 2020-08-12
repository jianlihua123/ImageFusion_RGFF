close all;
clear;
clc;

% source images path
IRDir = './SourceImages/IR_BMP/';
CCDDir = './SourceImages/VIS_BMP/';


output_path = './result/';


img_path_list = dir(strcat(IRDir, '*.bmp'));
img_num = length(img_path_list);

% parameter settings
lambda = 3; alpha = 1;
f = 6; scale = 3;
w1 = 27; sp1 = 30; ra1 = 300;
w2 = floor(w1/scale); sp2 = floor(sp1/scale); ra2 = floor(ra1/scale);
w3 = floor(w2/scale); sp3 = floor(sp2/scale); ra3 = floor(ra2/scale);
w4 = floor(w3/scale); sp4 = floor(sp3/scale); ra4 = floor(ra3/scale);

% main function
tic;
if img_num > 0   % from the third document
    for i = 1:img_num  % from the third document
        
        
        IR_image_name = img_path_list(i).name;% IR image name
        [filepath,name,ext] = fileparts(IR_image_name);
        index = regexp(name,'\d*\.?\d*','match');
        
        IR_image =  double(imread(strcat(IRDir,IR_image_name))); % IR image
        
        VS_image_name = strcat('VIS',num2str(index{1}),'.bmp');
        VS_image =  double(imread(strcat(CCDDir,VS_image_name))); % VIS image
        
        fprintf('%d %s\n',i,strcat(IR_image_name,'----',VS_image_name));
        
        [h ,w] = size(VS_image);
        I = zeros(h, w, 2);
        I(:,:,1)  = IR_image;
        I(:,:,2)  = VS_image;
        
        F = RGF_MultiScale_KIR_Operator(I, lambda, alpha, f,w1,sp1,ra1,w2,sp2,ra2,w3,sp3,ra3,w4,sp4,ra4);
        
        imwrite(F,strcat(output_path,strcat('F',num2str(index{1}),'.bmp')));        
        I = [];
    end
end
toc;

