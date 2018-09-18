---
title: "Homework Assignment 1"
date: 2018-09-18T12:34:57+09:00
draft: false
tags: ["Computational Photograpy"]
categories: []
---


# Initials

```Matlab
% Load the image into Matlab
im = imread('data/banana_slug.tiff');
figure, imshow(im); title('Loaded image');
```

This gives a dark image, shown as below.

![Loaded image](https://i.imgur.com/KDWwJI4.png)



```Matlab
% Check and report how many bits per integer the image has, and what its width and height is
disp(['Data type: ' class(im)])
disp(['Image height: ' num2str(size(im,1)) ', width: ' num2str(size(im,2))])
```
Now we can say its bit-per-integer, and image width/height:

`Data type: uint16`

`Image height: 2856, width: 4290`

And finally,

```Matlab
% convert the image into a double-precision array
im_double = double(im);
```

# Linearization
We convert the integer range [2047, 15000] to [0, 1], and cut out over the threshold.

```Matlab
% Convert the image into a linear array within the range [0, 1]
lim = [2047 15000];
im_lin = (im_double-lim(1)) / (lim(2)-lim(1));

% clip negative values to 0, and values greater than 1 to 1
im_lin = max(im_lin, 0);
im_lin = min(im_lin, 1);
```


# Identifying the correct Bayer pattern
```Matlab
% create three sub-images of im as shown in figure below
im1 = im_lin(1:2:end, 1:2:end);
im2 = im_lin(1:2:end, 2:2:end);
im3 = im_lin(2:2:end, 1:2:end);
im4 = im_lin(2:2:end, 2:2:end);

disp('Is the Bayes pattern XggX, or gXXg ?')
disp(['  MSE between UL pixels and DR pixels: ' num2str(immse(im1,im4))])
disp(['  MSE between UR pixels and DL pixels: ' num2str(immse(im2,im3)) newline])
```

`  MSE between UL pixels and DR pixels: 0.00013317`

`  MSE between UR pixels and DL pixels: 5.4164e-05`

MSE between UR-DL is smaller, and we can conclude the Bayes pattern is rather `XggX` pattern.
(Although, green is the strongest signal out of R,G and B. So it is not a bad strategy to select the brightest channel as green.)

Let's visualize each of `rggb` and `bggr`.

```Matlab
% combine the above images into an RGB image, such that im1 is the red, im2 is the green, and im3 is the blue channel
im_rggb = cat(3, im1, im2, im4);
im_bggr = cat(3, im4, im2, im1);

% (check each image, and) report which version you identified
figure, imshow(im_rggb * 5), title('Bayer pattern estimated: rggb');
figure, imshow(im_bggr * 5), title('Bayer pattern estimated: bggr');
```
For no reason I used `im*5` instead of `min(1, im*5)`. Each image is shown below.

* rggb

![Bayer pattern estimated: rggb](https://i.imgur.com/rC0HsX9.jpg)

* bggr

![Bayer pattern estimated: bggr](https://i.imgur.com/OPTu3nC.jpg)

Whatever, I thought `rggb` seems more natural (especially on the shirt), so I chose that one.

# White balancing

```Matlab
% Grey world assumption: force average color of scene to be grey.
avg_color = squeeze( mean(mean(im_rggb)) );
grey_balance(:,:,1) = avg_color(2)/avg_color(1);
grey_balance(:,:,2) = 1;
grey_balance(:,:,3) = avg_color(2)/avg_color(3);
im_gbal = grey_balance .* im_rggb;

% White world assumption: force brightest object in scene to be white.
max_color = squeeze( max(max(im_rggb)) );
white_balance(:,:,1) = max_color(2)/max_color(1);
white_balance(:,:,2) = 1;
white_balance(:,:,3) = max_color(2)/max_color(3);
im_wbal = white_balance .* im_rggb;

% check what the image looks like under both white balancing algorithms, and decide which one you like best
figure, imshow(im_gbal * 5); title('Grey-world assumption');
figure, imshow(im_wbal * 5); title('White-world assumption');
```

These are the results:

* Grey-world assumption

![Grey-world assumption](https://i.imgur.com/1nUWKdc.jpg)

* White-world assumption

![White-world assumption](https://i.imgur.com/b3elEQo.jpg)

Grey-balanced result looks good. I'd go with it.



# Demosicing

```Matlab
im_interp = zeros(size(im_lin));

% Interpolating red and blue
im_interp(:,:,1) = padarray( interp2(im_gbal(:,:,1)), [1 1], 'replicate', 'pre');
im_interp(:,:,3) = padarray( interp2(im_gbal(:,:,3)), [1 1], 'replicate', 'post');

% In case of green, to interpolate checkerboard pattern, I used inpaint_nans function implemented as here:
% https://www.mathworks.com/matlabcentral/fileexchange/4551-inpaint_nans
im_green = nan(size(im_lin));
im_green(1:2:end, 2:2:end) = im2;
im_green(2:2:end, 1:2:end) = im3;
im_interp(:,:,2) = inpaint_nans(im_green, 2);
```

Rather interpolating green with function `interp2`, I downloaded `inpaint_nans` function from matlab central(since I was lazy). Option 2 of `inpaint_nans` is used since it was the fastest option.

# Brightness adjustment and gamma correction
```Matlab
% Brighten the image by linearly scaling it by some percentage of the pre-brightening maximum grayscale value
max_gray_val = max(max(rgb2gray(im_interp)));
some_percentage = 4.5;
C_lin = im_interp * some_percentage/max_gray_val;
figure, imshow(C_lin); title(['C_{linear}, brightness adjusted image (by magnitude of ' num2str(some_percentage) ')']);
```
I chose 450% for, well, some reason. 100% (over `max_gray_val`) would be a intuitive choice, but a camera is in high dynamic range, thereby a white can be bright to infinity (opposed to black recieving no less than zero photon). Why the exact 450% though? No reason, just use 416% or 522% then. I won't care. Here's the result:

![C_{linear}, brightness adjusted image (by magnitude of 4.5)](https://i.imgur.com/ksLRnrL.jpg)


```Matlab
% tone reproduction (gamma correction)
mask_C = imbinarize(C_lin, .0031308);
C_nonlin = (1-mask_C).*(12.92 * C_lin) + mask_C.*( (1+.055)*(C_lin.^(1/2.4))-.055 );
figure, imshow(C_nonlin); title('C_{non-linear}, tone reproducted image');
```

And here's the final result:
![C_{non-linear}, tone reproducted image](https://i.imgur.com/IMTGimZ.jpg)

# Compression
```Matlab
% Use the imwrite command to store the image in .PNG format (no compression), and also in .JPEG format with quality setting 95
imwrite(C_nonlin, 'output.png');
imwrite(C_nonlin, 'output.jpg', 'Quality', 95);

% What is the compression ratio?
out_png = dir('output.png');
out_jpg = dir('output.jpg');
disp(['Compression ratio ( sizeOfJpg/sizeOfPng ) : ' num2str( out_jpg.bytes/out_png.bytes )])
```

Thanks to the high-tech JPEG compression (and high quality factor), We can barely notice any differnece between both image.

The size was 3.9MB versus 18.3MB, and

`Compression ratio ( sizeOfJpg/sizeOfPng ) : 0.21273`

