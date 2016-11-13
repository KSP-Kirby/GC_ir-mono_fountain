function [ rmsErr ] = computErrorGC( frontFlow, depthMap, minLabel, labels )
% 24 August 2016 Dual Flea Images Fountain GC
% call: computErrorGC(vl, depthMap, minLabel, labels) after running
% depthFromDualFocalLengthFlowDriverRight.m

left = 250;
right = 320;
top = 215;
bottom = 315;

depthMapXY = depthMap;

imshow(depthMapXY/length(labels));

fl = 1293.74844*4.8/6;


hold on
plot(left:right,top:top,'r')
plot(right:right,top:bottom,'r')
plot(left:right,bottom:bottom,'r')
plot(left:left,top:bottom,'r')
hold off;

stepSize = 1;

depth = [];
flow = [];
Vel = [];
if stepSize == 1
    index = 1;
    for i = left:right
        for j = top:bottom
            depth(index) = ((depthMapXY(j,i) - 1) + minLabel)*10;
            flow(index) = frontFlow(j,i);
            Vel(index) = flow(index).*depth(index)/fl;
            index = index + 1;
        end
    end
else
    index = 1;
    for i = left:right
        for j = top:bottom
            depth(index) = ((depthMapXY(j,i) - 1)*10 + minLabel)*10;
            flow(index) = frontFlow(j,i);
            Vel(index) = flow(index).*depth(index)/fl;
            index = index + 1;
        end
    end
end
    

vMean = 20;
rmsErr = sqrt(mean((Vel-vMean).*conj(Vel-vMean)));
disp(strcat('RMS Error:',num2str((rmsErr/20)*100),'%'))

vMean = mean(Vel);
rmsErr = sqrt(mean((Vel-vMean).*conj(Vel-vMean)));
disp(strcat('RMS Error:',num2str((rmsErr/20)*100),'%'))

end

