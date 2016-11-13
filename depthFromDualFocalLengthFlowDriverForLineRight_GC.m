% version 24 August 2016 IR-RGB images no gator
% for a line using the right image as the sensed image


try
    addpath('C:\Users\Richard\Documents\MATLAB\flow-code-matlab')     % this is flow to color
catch
    addpath('C:\Users\richa\Documents\MATLAB\flow-code-matlab')     % this is flow to color
end

try
    addpath('C:\Users\Richard\Documents\MATLAB\gco-v3.0\matlab')      % Graph Cuts
catch
    addpath('C:\Users\richa\Documents\MATLAB\gco-v3.0\matlab')      % Graph Cuts
end
load('flow.mat')

for smoothingFactor = [1]
    for dataFactor = [90]
        for horzNeighborMaskWeight = [5]
            for verticalNeighborMaskWeight = [2]
                %smoothingFactor = 10;    % higher is smoother, zero = no smoothing
                %dataFactor = 90;
                %horzNeighborMaskWeight = 5;
                %verticalNeighborMaskWeight = 2;

                minLabel = 60;
                maxLabel = 300;
                f_l = 6.21*.0048/.006;        % left camera is ir, calibrated 8/25/2016
                f_r = 6.21*.0048/.006;         % right camera is RGB
                pixelSize = .006;
                Vel = 20.0;
                b =  3.0*25.4;         % stereo baseline
                d =  300.0;         % dual focal length baseline

                %error.rms_Z should be multiplied by 100 to get percent depth error
                %error.rms_h is the pixel error, don't multiple by 100 and don't report as
                %a percent
                %error.dispErrPercent is the number of pixels that have an error greater
                %than 1

                % possible disparities range from 11 - 20, use labels 10 to 21.
                % 440 pixels (sites), 12 disparities (labels) one epipolar line
                % include the directory gco-v3.0/matlab
                % run GCO_UnitTest() to make sure everything works, otherwise you might
                % need to specify compiler per instructions in ??
                % run GCO_Delete(h) to get rid of object


                imageSet = 1;

                usePreviousDataCost = 0;

                startRow = 120;
                endRow = 120;
                vl = -uv_vl{imageSet}(startRow:endRow,:, 1);
                ir = -uv_ir{imageSet}(startRow:endRow,:, 1);

                [numRows,numCols] = size(vl);
                optimalLabellingOut = [];

                numSites = numRows*numCols;

                % labels are in cm
                % min label for the coaxial camera setup is about
                % max label is about
                labels = (minLabel:1:maxLabel)*10;  %*10 converts from cm to mm
                numLabels = length(labels);
                h = GCO_Create(numSites,numLabels);

                if usePreviousDataCost == 1
                    load('dataCost.mat')
                else
                    dataCost = zeros(numLabels, numSites);

                    site = 1;

                    h1 = waitbar(0, 'Constructing data cost matrix');
                    for k = 1:numRows
                        for i = 1:numCols
                            for j = 1: numLabels
                                Z = labels(j);
                                %if i + disparity <= numCols
                                    m = (f_r/f_l)*(Z+d)/(Z);
                                    rightPixel_1 = i;
                                    x_r = (rightPixel_1-320)*pixelSize;
                                    x_l = (x_r*Z*f_l+b*f_r*f_l)/(f_r*Z+f_r*d);
                                    leftPixel_1 = 320+x_l/pixelSize;        % this is a fractional pixel
                                    %rp(j) = rightPixel_1;
                                    if (leftPixel_1 < 1)
                                        dataCost(j, site) = 2;
                                    elseif (leftPixel_1 > 640)
                                        dataCost(j, site) = 2;
                                    else
                                        %wb_intrp = wb(k, floor(scaledPixelLocation)) + ((wb(k,ceil(scaledPixelLocation))-wb(k,floor(scaledPixelLocation))) * (scaledPixelLocation-floor(scaledPixelLocation)));                      % this interpolates between the two interger wb values
                                        %wb_intrp = wb(k,round(scaledPixelLocation));
                                        dataCost(j, site) = abs(m*ir(round(leftPixel_1)) - vl(rightPixel_1));

                                    end
                                    
                                %else
                                %    dataCost(j,site) = NaN;
                                %end   
                            end
                            site = site + 1;
                        end
                        waitbar(k/numRows)
                    end
                    close(h1)
                    
%                     figure
%                     plot(labels, energy)
%                     figure         
%                     plot(vl)
%                     hold all
%                     plot(ir)
%                     [val,index] = min(energy);
%                     plot(leftPixel_1,ir(leftPixel_1), '*r')
%                     plot(rp(index),vl(round(rp(index))), '*r')
%                     hold off

                    minDataCost = min(min(dataCost));
                    maxDataCost = max(max(dataCost));

                    %scale data cost between 1 and 1000

                    scaleFactor =50/maxDataCost;
                    dataCost = dataFactor*scaleFactor*dataCost+1;
                end


                GCO_SetDataCost(h,cast(dataCost,'int32'));

                % Smooth cost is number of labels by number of labels
                smoothCost = zeros(numLabels, numLabels);
                for i = 1:numLabels
                    smoothCost(i,:) = abs((1:1:numLabels) - i);
                end

                smoothCost = smoothCost*smoothingFactor;

                GCO_SetSmoothCost(h,cast(smoothCost, 'int32'));

                neighbors = sparse(numSites,numSites);

                colCount = 1;
                for i = 1:numSites-1
                    if i ~=numCols
                        neighbors(i,i+1) = horzNeighborMaskWeight;
                        colCount = colCount + 1;
                    else
                        colCount = 1;
                    end
                end

                % Radial line to radial line neighbors
%                 for i = 1:numSites-numCols - 1
%                     neighbors(i,i+numCols) = verticalNeighborMaskWeight;
%                 end

                % Tie together first and last ray
%                 for i = 1:numCols
%                     neighbors(i,numSites-numCols + i) = verticalNeighborMaskWeight;
%                 end

                GCO_SetNeighbors(h,neighbors);

                GCO_Expansion(h);

                [E, D, S] = GCO_ComputeEnergy(h); 

                optimalLabelling = cast(GCO_GetLabeling(h),'double');

                %figure
                %plot(optimalLabelling);

                GCO_Delete(h);
                
                plot(optimalLabelling)
                depthMap = optimalLabelling*10 + (minLabel-1)*10;
                m = (f_r/f_l)*(depthMap+d)./(depthMap);
                rightPixel_1 = 1:1:640;
                x_r = (rightPixel_1-320)*pixelSize;
                x_l = (x_r.*depthMap'*f_l+b*f_r*f_l)./(f_r*depthMap'+f_r*d);
                leftPixel_1 = 320+x_l/.006;
                
                for i = 1:640
                    if (leftPixel_1(i) >= 1) && (leftPixel_1(i) <=640)
                        ir_adj(i) = (m(i))*ir(round(leftPixel_1(i)));
                    end
                end
                
                plot(vl)
                hold all
                plot(ir)
                plot(ir_adj)
                hold off
                
                title('Flow along an epipolar line adjusted by estimated depth')
                legend('Visible Light Flow (VLF)','IR Flow (IRF)', '(1/c)*VLF(x+h)')
                xlabel('Pixel')
                ylabel('Flow (pixels/frame)')
                
                figure
                plot(vl)
                hold all
                plot(ir)    

                % evaluate at a pixel
                rightPixel_1 = 427;
                plot(rightPixel_1,vl(rightPixel_1),'g*')
                x_r = (rightPixel_1-320)*.006;
                for i = 1:length(labels)
                    z = labels(i);
                    x_l = (x_r*z*f_l+b*f_r*f_l)/(f_r*z+f_r*d);
                    leftPixel_1 = 320+x_l/pixelSize;
                    m = (f_r/f_l)*(z+d)/(z);
                    ir_flow(i) = (1/m)*vl(rightPixel_1);
                    pix(i) = leftPixel_1;
                    %plot(rightPixel_1,v1_flow(i), '*r')
                end
                plot(pix,ir_flow)

                %axis([180,405,5,35])
                title('Energy minimization search path')
                xlabel('Pixels')
                ylabel('Flow (pixels/frame)')
                legend('Flow 1', 'Flow 2',strcat('Eval Point:',num2str(rightPixel_1)), 'Eval Energy')
                
                Zest = (optimalLabelling + minLabel - 1)*10;
                ZfromFlow = Vel./vl/pixelSize*f_r;
                ZfromFlowIR = Vel./ir/pixelSize*f_l + d;
                figure
                plot(Zest)
                hold all
                plot(ZfromFlow)
                plot(ZfromFlowIR)
                
                %dispError 
                
                dispError = vl-ir_adj;
                nonOccludedDispError = [dispError(1:210),dispError(245:322),dispError(392:633)];
                rmsError = sqrt(mean((nonOccludedDispError).*conj(nonOccludedDispError)));
                
                disp(strcat('Non occluded velocity reconstruction error center row:',num2str(rmsError)))
                
                Zerror = ZfromFlow - Zest';
                rmsError = sqrt(mean((Zerror).*conj(Zerror)))/mean(Zest);
                
                      
            end
        end
    end
end
