close all
clear

arr_3d = preProcess("NHMS.jpeg");

bit_stream = bitStreamConvert(arr_3d, 4);


function arr_3d = preProcess(imageName)
    image = imread(imageName);
    image_arr = im2double(im2gray(image));
    
    fun = @(block_struct) dct2(block_struct.data);
    image_dct = blockproc(image_arr, [8 8], fun);
    
    image_max = max(image_dct, [], 'all');
    image_min = min(image_dct, [], 'all');
    
    scaled_arr = (image_dct - image_min) / (image_max - image_min);

    [m, n] = size(scaled_arr);
    arr_3d = reshape(scaled_arr, [8 8 (m*n)/64]);
end

function bit_stream = bitStreamConvert(arr_3d, N)
    col_group = reshape(arr_3d(:, :, 1:N), [64*N 1]);
    col_int = round(col_group *  255);
    bit_stream = reshape(int2bit(col_int, 8), [1 64*8*N]);

end