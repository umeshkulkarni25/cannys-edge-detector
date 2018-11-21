% author: Umesh Kulkarni
% netId: umk214
% This code to accept image-path as commandline argument 
% to compile and execute this code use command:
% matlab -nodesktop -nosplash -r "cannys_edge_detector <insert path to the image file>"
% make sure that matlab binary is added to the PATH in your
% Windows/Linux/Mac machine 
% Also the code stores all the generated images int he same directory as
% that of sourcefile
function cannys_edge_detector(imagePath)
    close all;
    image = imread(imagePath); 

    % display original image%
    fig = figure('Name', 'Original Image','NumberTitle', 'off');  
    imshow(image); 
    saveas(fig,'original_image.bmp');
    
    % convert image to a matrix of double for further matematical operations %
    image = double(image); 
    
    % smoothen the image by applying gaussian filtering %
    smoothened_image = apply_gaussian_filter(image);
    fig1 = figure('Name','gaussian filtred image','NumberTitle', 'off'); 
    imshow(uint8(smoothened_image));
    saveas(fig1,'gaussian_filtred_image.bmp');
    % apply edge detion operator to smoothened image %
    [Gm, Ga] = apply_prewitt_operator(smoothened_image);
    
    % apply non maxima suppression gradient magnitude image %
    normalized_thinned_image = apply_non_maxima_suppression(Gm,Ga);
    fig2 = figure('Name','non-maxima suppressed image','NumberTitle', 'off');
    imshow(uint8(normalized_thinned_image));
    saveas(fig2,'non_maxima_suppressed_image.bmp');
    
    % apply different thrshold and display the image %
    disp("Thresholding Results:");
    disp("----------------------------------------------------");
    threshold50_image = apply_thresholds(normalized_thinned_image, 50);
    fig3 = figure('Name','50% thresholded image','NumberTitle', 'off'); 
    imshow(threshold50_image); 
    saveas(fig3, 'ptile_50_thresholded_image.bmp');
    
    threshold30_image = apply_thresholds(normalized_thinned_image, 30);
    fig4 = figure('Name', '30% thresholded image','NumberTitle', 'off'); 
    imshow(threshold30_image); 
    saveas(fig4, 'ptile_30_thresholded_image.bmp');
    
    threshold10_image = apply_thresholds(normalized_thinned_image, 10);
    fig5 = figure('Name', '10% thresholded image','NumberTitle', 'off');
    imshow(threshold10_image); 
    saveas(fig5, 'ptile_10_thresholded_image.bmp');
    quit;
end

% function to take image as input and returns a gaussian filtered version
% of that image %
function [gaussian_filtered_image] = apply_gaussian_filter(image)
    gaussian_filter = [
        1, 1, 2, 2, 2, 1, 1;
        1, 2 ,2, 4, 2, 2, 1;
        2, 2, 4, 8, 4, 2, 2;
        2, 4, 8, 16, 8, 4, 2;
        2, 2, 4, 8, 4, 2, 2;
        1, 2 ,2, 4, 2, 2, 1;
        1, 1, 2, 2, 2, 1, 1;
        ];
    % apply convolution to the image %
    gaussian_filtered_image = convolution(image, gaussian_filter);

    % normalize the image by dividing it with sum of enties in the gaussian filter mask %
    mask_sum = sum(gaussian_filter, 'all');
    for row = 1: size(gaussian_filtered_image,1)
        for column = 1:size(gaussian_filtered_image,2)
            gaussian_filtered_image(row,column) = gaussian_filtered_image(row,column)/mask_sum;
        end
    end
end

% function to apply prewitt_operator to provided image and return gradient
% magnitude and gradient angle matrices %
function [Gm, Ga] = apply_prewitt_operator(image)
    prewitt_operator_Gx = ([
        -1,  0,  1;
        -1,  0,  1;
        -1,  0,  1
        ]);
    prewitt_operator_Gy = ([
        1,   1,   1;
        0,   0,   0;
       -1,  -1,  -1;
        ]);

    % calculate horizontal gradient %
    Gx = convolution(image, prewitt_operator_Gx);
    
    % normalize and display the horizontal gradient image %
    normalized_Gx = Gx(:,:);
    for row = 1:size(normalized_Gx,1)
        for column = 1:size(normalized_Gx,2)
            pixel_value = normalized_Gx(row, column);
            if pixel_value < 0 
                pixel_value = -1 * pixel_value;
            end
            normalized_Gx(row, column) = pixel_value / 3;
        end
    end
    fig6 = figure('Name', 'horizontal-gradient image','NumberTitle', 'off');
    imshow(uint8(normalized_Gx)); 
    saveas(fig6, 'horizontal_gradient.bmp');
    
    % calculate vertical gradient %
    Gy = convolution(image, prewitt_operator_Gy);
    
    % normalize and display the horizontal gradient image %
    normalized_Gy = Gy(:,:);
    for row = 1:size(normalized_Gy,1)
        for column = 1:size(normalized_Gy,2)
            pixel_value = normalized_Gy(row, column);
            if pixel_value < 0 
                pixel_value = -1 * pixel_value;
            end
            normalized_Gy(row, column) = pixel_value / 3;
        end
    end
    fig7 = figure('Name','vertical-gradient image', 'NumberTitle', 'off');
    imshow(uint8(normalized_Gy)); 
    saveas(fig7, 'vertical_gradient.bmp');

    % calculate gradient maginutde image %
    Gm = zeros(size(Gx));
    for row = 1:size(Gm ,1)
        for column = 1:size(Gm,2)
           Gm(row, column) = sqrt(Gx(row, column)^2+Gy(row, column)^2);
        end
    end
    
    % normalize and show gradient magnitude image %
    normalized_Gm = Gm./(sqrt(2)*3);
    fig8 = figure('Name','normalized gradient magnitude image', 'NumberTitle', 'off'); imshow(uint8(normalized_Gm)); 
    saveas(fig8,'gradient_magnitude.bmp');
    % calculate gradient angles %
    Ga = zeros(size(Gm));
    for row = 1:size(Gm,1)
        for column = 1:size(Gm,2)
            Ga(row, column) = atan(Gy(row, column)/Gx(row, column))*180/pi;
        end
    end
end

% function to apply non maxima suppresion to provided image and
% return the normalized thinned image %
function [normalized_thinned_image] = apply_non_maxima_suppression(Gm, Ga)
    thinned_image = Gm(:,:);
    % calulate nms operation indices using nms window size and the
    % size of the provided image %
    nms_window_rows = 3;
    nms_window_columns = 3;
    [image_rows, image_column] = size(Gm);
    nms_row_start = (nms_window_rows + 1)/2;
    nms_row_end = (image_rows - (nms_window_rows - 1)/2);
    nms_col_start = (nms_window_columns + 1)/2;
    nms_col_end = (image_column - (nms_window_columns - 1)/2);
    
    % apply non-maxima suppression %
    for row=nms_row_start:nms_row_end
        for column=nms_col_start:nms_col_end
            gradient_angle=Ga(row,column);
            if (gradient_angle >67.5 || gradient_angle<-67.5)
                if (Gm(row-1,column)>Gm(row,column)) ||(Gm(row+1,column)>Gm(row,column))
                    thinned_image(row,column)=0;
                end
            elseif(gradient_angle >=22.5 && gradient_angle<=67.5)
                if (Gm(row-1,column+1)>Gm(row,column)) || (Gm(row+1,column-1)>Gm(row,column))
                    thinned_image(row,column)=0;
                end
            elseif (gradient_angle >-22.5 && gradient_angle<22.5)
                if (Gm(row,column+1)>Gm(row,column)) || (Gm(row,column-1)>Gm(row,column))
                    thinned_image(row,column)=0;
                end
            elseif(gradient_angle >=-67.5 && gradient_angle<=-22.5)
                if (Gm(row+1,column+1)>Gm(row,column)) || (Gm(row-1,column-1)>Gm(row,column))
                    thinned_image(row,column)=0;
                end
            end

        end
    end
    
    % remove boundry pixels %
    thinned_image = thinned_image(nms_row_start:nms_row_end, nms_col_start:nms_col_end);
    % this normalization takes both prewitt-normalization and
    % magnitude calulation-normalization into consideration %
    normalized_thinned_image = thinned_image ./ (sqrt(2)*3);
end

function [thresholded_image] = apply_thresholds(image, ptile)
 thresholded_image = image(:,:);
 threshold = 0;
 % Q is number of non-zero pixel in the non-maxima suppressed image 
 % histogram_vector vector is needed to apply ptile method %
 [Q, histogram_vector] = get_histogram(image);
 
 % calculate threshold value for given value %
 number_of_pixels = round(Q * ptile /100);
 pixel_counter = number_of_pixels;
 for gray_value = 255:-1:1
     pixel_counter = pixel_counter - histogram_vector(1, gray_value);
     if pixel_counter <= 0 
         threshold = gray_value;
         break;
     end
 end

 disp(sprintf(" ptile= %d \n threshold= %d \n number_of_edge_pixels= %d", ptile, threshold, number_of_pixels));
 disp("----------------------------------------------------");
 % apply calculated threshold to the provided image %
  for row =1 : size(image,1)
     for column= 1 : size(image, 2)
         if image(row, column) < threshold
            thresholded_image(row, column) = 0;
         else
            thresholded_image(row, column) = 255;
         end
     end
  end
end

% calculate histogram of provided image 
% also return number of non-zero pixel s%
function [Q ,histogram_vector] = get_histogram(image)
 histogram_vector = zeros(1, 255);
 % we ceil the image to nearest integer as histogram will need whole numbers
 % we use ceil because we receive the non maxima suppressed image values 
 % such as 0.3, 0.16 are local mximas and they should not be rounded of to 
 % zero but to 1%
 image = ceil(image);
 
 % Q is number of non-zero pixels %
 Q = 0;
 % iterate through image to calculate histogram %
 for row = 1 : size(image,1)
     for column= 1 : size(image, 2)
         gray_value = image(row, column);
         if(gray_value > 0)
              Q = Q + 1;
              histogram_vector(1,gray_value) = histogram_vector(1,gray_value) + 1;
         end
     end
 end
end

% function to convolve any provided image with any provided filter %
function [filtered_image] = convolution (image, filter)
    filtered_image = image(:,:);
    % calulate from and to indices for the convolution operation based on 
    % provided filter size and image size %
    [filterRow, filterCol] = size(filter);
    [imageRow, imageCol] = size(filtered_image);
    convolutionRowStart = (filterRow + 1)/2;
    convolutionRowEnd = (imageRow - (filterRow - 1)/2);
    convolutionColStart = (filterCol + 1)/2;
    convolutionColEnd = (imageCol - (filterCol - 1)/2);
    for row = convolutionRowStart:convolutionRowEnd
        for col = convolutionColStart: convolutionColEnd
            filtered_image(row, col) = getWeightedSum(image, row, col, filter); 
        end
    end
    filtered_image = filtered_image(convolutionRowStart:convolutionRowEnd,convolutionColStart: convolutionColEnd);
end

% function to multiply the filter window with image at 
% a provided point of reference %
function [weightedSum] = getWeightedSum(image, refRow, refCol, filter)
    weightedSum = 0;
    [filterRow, filterCol] = size(filter);
    refRow = refRow - (filterRow-1)/2;
    refCol = refCol - (filterCol-1)/2;
    for row = 1:filterRow
        for col = 1: filterCol
            imRow = refRow+row-1;
            imCol = refCol+col-1;
            weightedSum = weightedSum + filter(row, col)*image(imRow,imCol);
        end
    end
end