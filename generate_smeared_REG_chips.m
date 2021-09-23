image = [];
diameter = 100;
shift_mag = 1;
num_shifts = 1;
num_zeroes = 0;
change_percent = 1;
check_REG = NaN;
REG = 0;
p_revert_color = 0;
chip_num = 15; 

%make directory if it doesn't already exist
if ~exist("~/Desktop/REG_chips/" + num2str(REG) + "/")
    mkdir("~/Desktop/REG_chips/" + num2str(REG) + "/")
end

for l = 1:chip_num %repeat the code until the desired number of chips is created
    check_REG = 2;
    while check_REG ~= REG
        chip_pixels = [repelem(0, diameter/2), repelem(1, diameter/2)]; %create a vector with the correct number of black and white pixels
        chip_pixels = chip_pixels(randperm(length(chip_pixels))); %randomize pixel order

        for k = 1:diameter
            x = sqrt((diameter/2)^2 - (diameter/2 - k)^2);
            num_pixels_per_row(k) = 2*x;
        end

        num_pixels_per_row = ceil(num_pixels_per_row);
        image = [];

        for i = 1:diameter

            shift = randsample([repelem(0, num_zeroes), repelem(shift_mag, num_shifts), repelem(-shift_mag, num_shifts)], 1);

            if shift > 0
                move_to_top = chip_pixels(diameter - shift + 1:diameter);
                keep = chip_pixels(1:diameter - shift);

                chip_pixels(1:length(move_to_top)) = move_to_top;
                chip_pixels(length(move_to_top) + 1:end) = keep;
            end

            if shift < 0 
                move_to_bottom = chip_pixels(1:-shift);
                keep = chip_pixels(-shift + 1:diameter);

                chip_pixels(1:diameter + shift) = keep;
                chip_pixels(diameter + shift + 1:end) = move_to_bottom;
            end

            change_pixel_value = [repelem(0, floor(change_percent*(diameter/100))) repelem(1, ceil(diameter - change_percent*(diameter/100)))];
            change_pixel_value = change_pixel_value(randperm(length(change_pixel_value)));
            changed_values = find(change_pixel_value == 0);

            for j = 1:length(changed_values)
                if chip_pixels(changed_values(j)) == 1
                    chip_pixels(changed_values(j)) = 0;
                else
                    chip_pixels(changed_values(j)) = 1;
                end
            end

            image(:, i) = chip_pixels;

            revert_vec = [repelem(1, ceil(length(changed_values*p_revert_color))), repelem(0, floor(length(changed_values)- length(changed_values*p_revert_color)))];
            revert_vec = revert_vec(randperm(length(revert_vec)));

            for k = 1:length(changed_values)
                if revert_vec(k) == 1
                    if chip_pixels(changed_values(j)) == 1
                        chip_pixels(changed_values(j)) = 0;
                    else
                        chip_pixels(changed_values(j)) = 1;
                    end
                end
            end

        end

        %build in buffer
        buffer_sizes = diameter - num_pixels_per_row; %create row buffers
        for n = 1:length(buffer_sizes)
            each_buffer = buffer_sizes(n) / 2;
            buffer_1 = ceil(each_buffer);
            buffer_2 = floor(each_buffer);
            image(n, 1:buffer_1) = repelem(0.75, buffer_1);
            image(n, end - buffer_2 + 1:end) = repelem(0.75, buffer_2);
        end

        check_regularity_vector = [];
        for m = 1:size(image, 1) %go through each row and calculate regularity
            row_pixels = image(m,:); %select the row of pixels
            avg_pixel_value = sum(row_pixels == 1) / (sum(row_pixels == 1)+sum(row_pixels == 0)); %avg pixel values is occurences of white / total pixels
            check_row_regularity = 2*abs(0.5 - avg_pixel_value); %regularity is absolute value of total pixel value*2
            check_regularity_vector(m) = check_row_regularity; %store in a vector
        end

        check_REG = mean(check_regularity_vector, 'omitnan'); %overall regularity is the mean of all row regularities. This value gets returned. 
        check_REG = round(check_REG, 2); %round to the nearest hundredth
        display(check_REG) %check if this is equal to desired regularity. If not, try again. 
    end

    imwrite(image, "~/Desktop/REG_chips/" + num2str(REG) + "/" + num2str(check_REG) + "_" +  num2str(l) + ".tif")
    save("~/Desktop/REG_chips/" + num2str(REG) + "/" + num2str(check_REG) + "_" + num2str(l) + ".mat");

    if l ~= chip_num %clear all the relevant parameters if not the last cycle. If the last cycle, return image, check_REG, and row_regularities
        clear image check_REG row_regularities
    end
end