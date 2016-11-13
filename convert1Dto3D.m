function [ imgOut ] = convert1Dto3D( imgIn )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

        [m,n,p] = size(imgIn);
        if p == 1
            imgOut = zeros(m,n,3);
            imgOut = cast(imgOut,'uint8');
            imgOut(:,:,1) = imgIn;
            imgOut(:,:,2) = imgIn;
            imgOut(:,:,3) = imgIn;
        end
        
end

