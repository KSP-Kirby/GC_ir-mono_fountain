%load('flow.mat')
constrast = 5;
offset = .5;
imgs = cast(uvi_ir{1},'double')/256;
imgs1 = imgs(:,:,1) - offset;
imgs2 = imgs1*constrast;
imgs2(:,:,2) = imgs(:,:,2)*constrast;
imgs2(:,:,3) = imgs(:,:,3)*constrast;

imtool(imgs2)