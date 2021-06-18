%% generate_REG_chips.m
%Written by JP Dinh on 2021 June 17

% Randomly generates circles of any desired regularity as calculated by Gluckman & Cardosa 2009. 
% From the paper: "A measure of regularity (henceforth REG) is calculated from the 
% aligned image as follows. Each horizontal row of pixels is scored as twice the absolute 
% value of the difference between its average pixel value (black pixels scored as 1 and 
% white as 0) and 0.5; i.e., a row of pixels is scored as 1 when all pixels are black or 
% all are white, and is scored progressively less, until zero, when half the pixels are 
% black and half are white; REG is the average of this score across all rows. Many rows of 
% very regular patterns are nearly monochromatic, and receive high REG scores, while 
% irregular patterns have many rows with a mixture of black and white pixels, and 
% consequently receive lower REG scores (e.g., Fig. 2). 

%How to use: Set the desired regularity (REG), variance (var), and diameter
%of the circle in pixels (diameter). This is in lines 20, 21, and 22. Click
%Run. The program will save a .tif on your desktop of the generated circle.

%% Code to execute the function. Fill in the first three lines to set desired conditions for the image. 
REG = 0.1; %This is the regularity we want to generate
var = 0.09; %variance of regularities between rows
%chips in previous papers were 2.5 cm = 0.984 inches. 600 dpi x 0.984
%inches = 590 pixel image
diameter = 300; %diameter of the chip in pixels
chip_num = 10; %number of unique chips you want to generate

[image, check_regularity_vector, check_REG] = generate_REG_chip(REG, var, diameter, chip_num); %Execute the function below

 %% function used to create chips of a desired regularity. 
function [image, row_regularities, check_REG] = generate_REG_chip(REG, var, diameter, chip_num)
    %make directory if it doesn't already exist
    if ~exist("~/Desktop/REG_chips/" + num2str(REG) + "/")
        mkdir("~/Desktop/REG_chips/" + num2str(REG) + "/")
    end
    
    for n = 1:chip_num %repeat the code until the desired number of chips is created
        check_REG = 2; %keep this some value that will never be used - just to initiate the variable
        while check_REG ~= REG %because of rounding, the regularity isn't always perfect. Retry until desired regularity is reached
            num_pixels_per_row = [];
            for k = 1:diameter
                x = sqrt((diameter/2)^2 - (diameter/2 - k)^2);
                num_pixels_per_row(k) = 2*x;
            end

            %randomly generate rows regularities
            rounded_regularity = 0; %initiate variable for while loop
            while rounded_regularity ~= REG %end the loop when a regularity that is desired is reach (not perfect because randomly generated)
                row_regularities = var.*randn(diameter, 1) + REG; %randomly generate regularities for each row
                for l = 1:length(row_regularities)
                    if row_regularities(l) > 1 % Regularity cannot be greater than 1
                       row_regularities(l) = 1;        	
                    end
                    if row_regularities(l) < 0 %Regularity cannot be less than 0
                       row_regularities(l) = 0;        	
                    end
                end
                mean_regularity = mean(row_regularities); %Calculate regularity
                rounded_regularity = round(mean_regularity, 2); %Round to the nearest hundredth
            end

            image = [];
            %Generate image based on each row's regularity

            if REG == 0 %for some reason, this variable dissapears if REG = 0. Manually set the regularities to 0 if desired. 
                row_regularities = repelem(0, diameter);
            end

            for j = 1:length(row_regularities) 
                black_or_white = randi([0, 1], length(row_regularities)); % Regularity has an absolute value - this line randomly determines if the row is mostly white or mostly black
                if black_or_white(j) == 1 %one of two absolute value conditions
                    avg_pixel_value = -1*row_regularities(j)/2 + 0.5; %average pixel value
                    white_pixels = ceil(num_pixels_per_row(j)*avg_pixel_value); %number of white pixels the avg pixel value*number of pixels in that row (rounded)
                    black_pixels = ceil(num_pixels_per_row(j) - white_pixels); %number of black pixels are the ones remaining
                else %one of two absolute value conditions
                    avg_pixel_value = row_regularities(j)/2 + 0.5;
                    white_pixels = ceil(num_pixels_per_row(j)*avg_pixel_value); %number of white pixels the avg pixel value*number of pixels in that row (rounded)
                    black_pixels = ceil(num_pixels_per_row(j) - white_pixels); %number of black pixels are the ones remaining
                end

                chip_pixels = [repelem(0, white_pixels), repelem(1, black_pixels)]; %create a vector with the correct number of black and white pixels
                chip_pixels = chip_pixels(randperm(length(chip_pixels))); %randomize pixel order

                buffer_size = diameter/2 - length(chip_pixels)/2; %create row buffers
                image(j,:) = [repelem(0.75, ceil(buffer_size)), chip_pixels, repelem(0.75, floor(buffer_size))]; %for each row, create the pixel values
            end
            
            %Check Regularity
            check_regularity_vector = [];
                for m = 1:size(image, 1) %go through each row and calculate regularity
                    row_pixels = image(m,:); %select the row of pixels
                    avg_pixel_value == sum(row_pixels == 1) / (sum(row_pixels == 1)+sum(row_pixels == 0)); %avg pixel values is occurences of white / total pixels
                    check_row_regularity = 2*abs(0.5 - avg_pixel_value); %regularity is absolute value of total pixel value*2
                    check_regularity_vector(m) = check_row_regularity; %store in a vector
                end
                check_REG = mean(check_regularity_vector); %overall regularity is the mean of all row regularities. This value gets returned. 
                check_REG = round(check_REG, 2); %round to the nearest hundredth
                display(check_REG) %check if this is equal to desired regularity. If not, try again. 
        end
        
        average_pixel_value = sum(image(:) == 1) / (sum(image(:) == 1) + sum(image(:) == 0));%calculate the average pixel value just in case it affects perception
        
        imwrite(image, "~/Desktop/REG_chips/" + num2str(REG) + "/REG_" + num2str(REG) + "_Chip" + num2str(n) + ".tif");
        save("~/Desktop/REG_chips/" + num2str(REG) + "/REG_" + num2str(REG) + "_Chip" + num2str(n) + ".mat");

        if n ~= chip_num %clear all the relevant parameters if not the last cycle. If the last cycle, return image, check_REG, and row_regularities
            clear image check_REG row_regularities
        end
    end
end