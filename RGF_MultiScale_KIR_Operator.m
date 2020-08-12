function [F] = RGF_MultiScale_KIR_Operator(I, lambda, alpha,factor1,w1,sp1,ra1,w2,sp2,ra2,w3,sp3,ra3,w4,sp4,ra4)

%% Saliency map obtaining

M1 = kirschedge(I(:,:,1));
M2 = kirschedge(I(:,:,2));

H(:,:,1) = M1;
H(:,:,2) = M2;

% Saliecny map comparison

S = GauSaliency(H);
P = IWconstruct(S);

%% Weigth maps for detail layers

% Weight Optimization with Joint Bilateral Filtering

for i = 1:1
    L = I(:,:,i);
    L = L./max(L(:));
    
    W_D3(:,:,i) = jbfilter2(P(:,:,i),L,w1,sp1,ra1);
    W_D2(:,:,i) = jbfilter2(P(:,:,i),L,w2,sp2,ra2);
    W_D1(:,:,i) = jbfilter2(P(:,:,i),L,w3,sp3,ra3);
    W_D0(:,:,i) = jbfilter2(P(:,:,i),L,w4,sp4,ra4);
    
end

% Optimal Correction for weight maps of detail layers

W_D3(:,:,1) = Gammacorrect( W_D3(:,:,1) ,alpha );
W_D3(:,:,2) = 1- W_D3(:,:,1);
W_D2(:,:,i) = Gammacorrect( W_D2(:,:,1) ,alpha );
W_D2(:,:,2) = 1- W_D2(:,:,1);
W_D1(:,:,i) = Gammacorrect( W_D1(:,:,1) ,alpha );
W_D1(:,:,2) = 1- W_D1(:,:,1);
W_D0(:,:,i) = Gammacorrect( W_D0(:,:,1) ,alpha );
W_D0(:,:,2) = 1- W_D0(:,:,1);


%% Fuse image reconstruction

F = GuidFuseWithS_new2(I,W_D3,W_D2, W_D1, W_D0,lambda, factor1);

F = uint8(F*255);



function [ S ] = GauSaliency( H )
% Using the local average of the absolute value of H to construct the
% saliency maps
N = size(H,3);
S = zeros(size(H,1),size(H,2),N);
for i=1:N
    se = fspecial('gaussian',11,5);
    S(:,:,i) = imfilter(H(:,:,i),se,'replicate');
end
S = S + 1e-12; %avoids division by zero
S = S./repmat(sum(S,3),[1 1 N]);%Normalize the saliences in to [0-1]

function [P] = IWconstruct( S )
% construct the initial weight maps
[r ,c ,N] = size(S);
[X ,Labels] = max(S,[],3); % find the labels of the maximum
clear X
for i = 1:N
    mono = zeros(r,c);
    mono(Labels==i) = 1;
    P(:,:,i) = mono;
end


function [F] = GuidFuseWithS_new2(I,W_D3,W_D2, W_D1, W_D0, lambda, factor1)
I = double(I)/255;

se = fspecial('average', [31 31]);

if size(I,3) == 3
    [r,c,M,N] = size(I);
    
    F_D = zeros(r,c,M);
    for n = 1:N
        
        w_D3 = W_D3(:,:,n);
        w_D2 = W_D2(:,:,n);
        w_D1 = W_D1(:,:,n);
        w_D0 = W_D0(:,:,n);
        
        G = I(:,:,:,n);
        B = tsmooth(G,0.015,3);
        D = G - B;
        
        F_D = F_D + D.*repmat(w_D,[1 1 3]);
    end
else
    [r,c,N] = size(I);
    F_B = zeros(r,c);
    F_D3 = zeros(r,c);
    F_D2 = zeros(r,c);
    F_D1 = zeros(r,c);
    F_D0 = zeros(r,c);
    for n=1:N
        
        w_D3 = W_D3(:,:,n);
        w_D2 = W_D2(:,:,n);
        w_D1 = W_D1(:,:,n);
        w_D0 = W_D0(:,:,n);
        
        factor = factor1;
        B1(:,:,n) = RollingGuidanceFilter(I(:,:,n), lambda, 0.2 );
        B2(:,:,n) = RollingGuidanceFilter(I(:,:,n), factor*lambda,   0.2);
        B3(:,:,n) = RollingGuidanceFilter(I(:,:,n), factor*factor*lambda,  0.2);
        
        B4(:,:,n) = RollingGuidanceFilter(I(:,:,n), factor*factor*factor*lambda,  0.2);
        u0 = I(:,:,n);  u1 = B1(:,:,n);
        u2 = B2(:,:,n); u3 = B3(:,:,n);
        u4 = B4(:,:,n);
        
        d0 = u0 - u1;
        d1 = u1 - u2;
        d2 = u2 - u3;
        d3 = u3 - u4;
        
        
        
        %% Fusion for detail layers
        
        F_D3 = F_D3 + d3.*w_D3;
        F_D2 = F_D2 + d2.*w_D2;
        F_D1 = F_D1 + d1.*w_D1;
        F_D0 = F_D0 + d0.*w_D0;
        
    end
    
    %% Fusion for base layer
    x1 = B4(:,:,1);x2 = B4(:,:,2);V1 = var(x1(:));V2 = var(x2(:));wb = V1./(V1 + V2);
    F_B = wb.*B4(:,:,1) + (1 - wb).*B4(:,:,2);
    
end
F = F_B + F_D3 + F_D2 + F_D1 + F_D0;

