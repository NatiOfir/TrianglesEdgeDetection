clc;
clear;
close all;
profile on;
set(0,'DefaultFigureWindowStyle','docked');
warning('off', 'images:imshow:magnificationMustBeFitForDockedFigure');

param = 0;

disp('Creating Masks');
tic;
if ~exist('Files/l0.mat', 'file') | true
    createPatchMasks(param);
end

disp('Building Tree');
str = {'Sqr','CC','Sines'};
n = 2^3+1;
%[img gt]=LinesCosntructSinimage(n,2);
%img = circle(30,129);
%img = (im2double(img)-0.5)*0.4+0.5;
sigma = 0.1;

%for i = 1:length(str);

 %   img = imread(sprintf('%s.png',str{i}));
    %img = imresize(img,[n n]);
    img = zeros(n)+0.5;
    img(end/2-0.5:end,:) = 0.6;
    %img = img+sigma*randn(n);
    im = Image(img,param,sigma);
    im = im.buildTree(true);
    %im = im.detectEdges(false);
    im = im.detectEdgesPlusPlus();

    disp('Displaying Edges');
    figure; subplot(1,2,1);imshow(im.resIgray);
    subplot(1,2,2);imshow(im.I);
    toc;
    
 %   imwrite(im.resIgray,sprintf('%s_res.png',str{i}),'PNG');
  %  break;
%end
profile viewer;
%p = profile('info');
%profsave(p,sprintf('profile_results_%s',datestr(clock)));
profile off;
