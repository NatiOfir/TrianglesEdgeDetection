clc;
clear;
close all;

param = 0;

disp('Creating Masks');
tic;
if ~exist('Files/l0.mat', 'file')
    createPatchMasks(param);
end

disp('Building Tree');

I = im2double(imread('Images/CC.png'));
fastRun = true;
if fastRun
    I = imresize(I,[65 65]);
end
sigma = 0.1;
[res,im] = runIm(I,false,param);
disp('Displaying Edges');
figure; 
subplot(1,2,1);imshow(res);
subplot(1,2,2);imshow(I);
toc;

