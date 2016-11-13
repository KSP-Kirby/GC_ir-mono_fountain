% computes flow using brox
% 24 Aug 2016 data, Dual flea images

addpath 'C:\Users\Richard\Documents\MATLAB\pami2010Matlab';
%         C:\Users\Richard\Documents\MATLAB\pami2010Matlab
addpath 'C:\Users\Richard\Documents\MATLAB\flow-code-matlab';

h = waitbar(0,'Processing Data')
images = [1, 2, 3, 4, 5];
%images = [4, 5, 6];

%b = 55.84;
fl_f = [ 1285.77473   1284.27991 ] * .006;  % for front camera use y
fl_b = [ 973.17464   972.30769  ] * .006;   % for back camera use x

%for jj = 1:length(method)
vl = 'm';
ir = 'ir';

framesBetweenFlowComp = 4;

    for j = 1:length(images)-framesBetweenFlowComp
        i = images(j);
        i1 = images(j+framesBetweenFlowComp);
        if i < 10
            vl_1=imread(strcat(vl,'0',num2str(i),'_rect.tif'));
            ir_1=imread(strcat(ir,'0',num2str(i),'_rect.tif'));
            titlef1 = strcat(vl,'0',num2str(i));
            titleb1 = strcat(ir,'0',num2str(i));
        else
            vl=imread(strcat(vl,num2str(i),'_rect.tif'));
            ir=imread(strcat(ir,num2str(i),'_rect.tif'));
            titlef1 = strcat(vl,num2str(i));
            titleb1 = strcat(ir,num2str(i));
        end

        if i1 < 10
            vl_2=imread(strcat(vl,'0',num2str(i1),'_rect.tif'));
            ir_2=imread(strcat(ir,'0',num2str(i1),'_rect.tif'));
            titlef2 = strcat(vl,'0',num2str(i1));
            titleb2 = strcat(ir,'0',num2str(i1));
        else
            vl_2=imread(strcat(vl,num2str(i1),'_rect.tif'));
            ir_2=imread(strcat(ir,num2str(i1),'_rect.tif'));
            titlef2 = strcat(vl,num2str(i1));
            titleb2 = strcat(ir,num2str(i1));
        end
        
        vl_1 = convert1Dto3D(vl_1);
        vl_2 = convert1Dto3D(vl_2);
        ir_1 = convert1Dto3D(ir_1);
        ir_2 = convert1Dto3D(ir_2);
        
        
         [m, n, p] = size(vl_1);
         upperLeftCornerM = (m - 600)/2;
         upperLeftCornerN = (n - 800)/2;

         ir_1 = imresize(ir_1,[round(m*(1000/700)),n]);
         ir_2 = imresize(ir_2,[round(m*(1000/700)),n]);
         
         m1 = size(ir_1,1);
         m2 = round((m1 - m)/2);
         ir_1 = imcrop(ir_1,[0,m2,n,m-1]);
         ir_2 = imcrop(ir_2,[0,m2,n,m-1]);
         
         vl_1 = imcrop(vl_1,[upperLeftCornerN, upperLeftCornerM, 799, 599]);
         vl_1 = imresize(vl_1, .8);
         
         vl_2 = imcrop(vl_2,[upperLeftCornerN, upperLeftCornerM, 799, 599]);
         vl_2 = imresize(vl_2, .8);         
         
         ir_1 = imcrop(ir_1,[upperLeftCornerN, upperLeftCornerM, 799, 599]);
         ir_1 = imresize(ir_1, .8);
         
         ir_2 = imcrop(ir_2,[upperLeftCornerN, upperLeftCornerM, 799, 599]);
         ir_2 = imresize(ir_2, .8);
         

        plotTitle = {strcat('From 20 Aug 2016 Images:', titlef1,'-',titlef2,'<<->>', titleb1, '-', titleb2);...
            strcat('method: Brox')};
        
        trimX = 80;
        trimY = 80;
        horzLine = (480 - 2*trimX)/2;

%         vl_1 = mirrorHorz(vl_1);
%         vl_2 = mirrorHorz(vl_2);
        
        %uv_vl{i} = estimate_flow_interface(vl(100:end-100,100:end-100,:), vl2(100:end-100,100:end-100,:),method{jj});
        %uv_vl{j} = mex_LDOF(double(vl(trimX:end-trimX,trimY:end-trimY,:)), double(vl2(trimX:end-trimX,trimY:end-trimY,:)));
        uv_vl{j} = mex_LDOF(double(vl_1), double(vl_2));
        waitbar((j-.5)/length(images))

        uvi_vl{j} = uint8(flowToColor(uv_vl{j}));

        %uv_ir{i} = estimate_flow_interface(bc(100:end-100,100:end-100,:), bc2(100:end-100,100:end-100,:),method{jj});
        %uv_ir{j} = mex_LDOF(double(bc(trimX:end-trimX,trimY:end-trimY,:)), double(bc2(trimX:end-trimX,trimY:end-trimY,:)));
        uv_ir{j} = mex_LDOF(double(ir_1), double(ir_2));
        uvi_ir{j} = uint8(flowToColor(uv_ir{j}));

        curFigure = figure;
        plot(uv_vl{j}(horzLine,:,1),'LineWidth',2)
        hold all
        plot(uv_ir{j}(horzLine,:,1),'LineWidth',2)
        title(plotTitle);
        legend('RGB', 'IR')
        ylabel('pixels')
        xlabel('Flow in pixels along center horizontal line')
        %axis([0,450,0,20])
        hold off
        filename = strcat('brox',num2str(images(j)),'.jpg');
        saveas(curFigure,filename)
        filename = strcat('brox',num2str(images(j)),'.fig');
        saveas(curFigure,filename)

        waitbar((j)/length(images))
    end
    
%end

close(h)

save('flow','uv_ir','uv_vl','uvi_ir','uvi_vl')
clear 'uv_ir'
clear 'uv_vl'
clear 'uvi_ir'
clear 'uvi_vl'
close all
