---
title: "Homework Assignment 4"
date: 2018-11-14T07:48:06+09:00
draft: false
tags: ["Computational Photograpy"]
---

Before start: I was able to do only the first two section (unsuccessfully). 


# Initials

```matlab```

addpath('src');
addpath('dcraw_farhi'); % https://www.mathworks.com/matlabcentral/fileexchange/66927-read-raw-camera-images

stack_path = 'exposure_stack';
weight_scheme = 'uniform'; % 'uniform' or 'hat'

%% Read LDR image stack
clearvars s_in
S_nef = dir(fullfile(stack_path, '*.nef'));
for k = numel(S_nef):-1:1
    F_nef = fullfile(stack_path, ['exposure' num2str(k) '.nef']);
    F_jpg = fullfile(stack_path, ['exposure' num2str(k) '.jpg']); % to read in (reverse) order.
    dc = readraw;
    s_nef(:,:,:,k) = imread(dc, F_nef, '-w T -6 -q 3');
    s_jpg(:,:,:,k) = imread(dc, F_jpg);
end

```


We read both raw(.nef) and rendered(.jpg) image stack. A Matlab implementation of `dcraw` by farhi ( https://www.mathworks.com/matlabcentral/fileexchange/66927-read-raw-camera-images ) is used.

Weight scheme can be chosen, as shown in the code.


# Linearize rendered images


```matlab```

function [g,lE]=gsolve(Z,B,l,w)

n = 256;

A = sparse(size(Z,1)*size(Z,2)+n+1,n+size(Z,1));
b = zeros(size(A,1),1);

%% Include the data?fitting equations

k = 1;
for i=1:size(Z,1)
	for j=1:size(Z,2)
        wij = w(Z(i,j)+1);
        
        A(k,Z(i,j)+1) = wij;
        A(k,n+i) = -wij;
        b(k,1) = wij * B(j);
        
        k=k+1;
	end
end

%% Fix the curve by setting its middle value to 0

A(k,129) = 1;
k=k+1;

%% Include the smoothness equations
for i=1:n-2
	A(k,i)=l*w(i+1);
    A(k,i+1)=-2*l*w(i+1);
    A(k,i+2)=l*w(i+1);
    
	k=k+1;
end

%% Solve the system using SVD

x = A\b;

g = x(1:n);
lE = x(n+1:size(x,1));


```

The code above is `gsolve.m`. Basically copied from Debevec's paper, with slight modificaions and fixed typos.

```matlab```
% resize image stack. We will only use g anyway.
r=64;
for j = size(s_jpg,4):-1:1
    s_in(:,:,:,j) = imresize(s_jpg(:,:,:,j),1/r);
end

ZL = reshape(s_in, size(s_in, 1)*size(s_in, 2), 3, size(s_in,4));
Z = squeeze(ZL(:,1,:));

B = -11 : size(s_in,4)-11-1;
B_exp = 2.^B;
l = 20;

if strcmp(weight_scheme, 'uniform')
	w = [zeros(3,1); ones(250,1); zeros(3,1)];
else
    w = [1:128 128:-1:1]';
end

for i = 3:-1:1
    Z = squeeze(ZL(:,i,:));
    [g, IE] = gsolve(Z, B, l, w); % eq. (2-3)
    
    % use g to convert the non-linear images into linear ones
     s_jl = exp( g(s_jpg(:,:,i,:)+1)/255 ); % eq. (4)
     s_jpglin(:,:,i,:) = uint8((s_jl-min(min(min(s_jl))))/(max(max(max(s_jl)))-min(min(min(s_jl))))*255);
    
    g_all(:,i) = g;
end

```



Although I resize image by 1/64 (to 63x94), size of A is like 95009x6178. Still dumb, but we can barely handle that using sparse matrix. Each of the channel takes a dozen seconds to be approximated by solving one big least-square problem.

I use the smoothness parameter `l = 20` which is big, since we're using downsized image which can cause more error.

Finally we can plot the function g, by each channel:

```matlab```

figure, h = plot(g_all);
set(h,{'color'}, {[1 0 0]; [0 1 0]; [0 0 1]});

```

The result below is got with uniform weight scheme.

![image](https://i.imgur.com/EZgBQeZ.png)


This one is with hat scheme.

![image](https://i.imgur.com/KSyiQxS.png)


Weight fluctuates a lot. It doesn't even monotonically increases, which is unexpected (and pretty annoying). 


# Merge exposure stack into HDR image

```matlab```
im_h = size(s_nef, 1);
im_w = size(s_nef, 2);

I_raw_lin = zeros(im_h, im_w, 3);
I_ren_lin = zeros(im_h, im_w, 3);
I_raw_log = zeros(im_h, im_w, 3);
I_ren_log = zeros(im_h, im_w, 3);

B_ = permute(repmat(B',[1 im_h im_w 3]),[2 3 4 1]);
B_exp_ = permute(repmat(B_exp',[1 im_h im_w 3]),[2 3 4 1]);

I_ = s_nef;
w_I = w(I_+1);
I_raw_lin = sum(w_I.*double(I_)./B_exp_, 4) ./ sum(w_I, 4);
I_raw_log = exp( sum(w_I.*(double(I_)/255-B_), 4) ./ sum(w_I, 4) );

I_ = s_jpglin;
w_I = w(I_+1);
I_ren_lin = sum(w_I.*double(I_)./B_exp_, 4) ./ sum(w_I, 4);
I_ren_log = exp( sum(w_I.*(double(I_)/255-B_), 4) ./ sum(w_I, 4) );

```

This code generates a total of 4 HDR images, 2 from each LDR stack. (We already fixed the weighting scheme on the first part of the code.)

We're dealing with 4016*6016*3*16 = 1160M elements. A quad-nested for loop will take ~1hr per stack. But we can cut it to a couple of minutes, using vectorization. Yet, a simple `w_I = w(I_+1);` takes most of the time.

&nbsp;


# Evaluation

Left column is uniform weighting, and right column is hat weighting.


### RAW - linear merging

<img src="https://i.imgur.com/VQmSN9C.png" width="350"><img src="https://i.imgur.com/cPSqY7C.png" width="350">


### RAW - logarithmic merging

<img src="https://i.imgur.com/hlZb3rG.png" width="350"><img src="https://i.imgur.com/CYAxUwx.png" width="350">


### Rendered - linear merging

<img src="https://i.imgur.com/KiJ1HcO.png" width="350"><img src="https://i.imgur.com/ptSVJ3J.png" width="350">


### Rendered - logarithmic merging

<img src="https://i.imgur.com/zBWnf5W.png" width="350"><img src="https://i.imgur.com/TKymrS4.png" width="350">

&nbsp;


You can see original images at: https://drive.google.com/open?id=1EynNuuWaDgnnvbSEq2C2kUA0LRzX7Gok


Seriously, no image was better than `lE` (output of `gsolve.m`). It should be a code bug, and here are my reasonings of the bug:

* the function g should be monotonically increasing. Somehow it didn't, and the low part fluctuation made a lotta noise.
* Debevec's code works on `uint8` range, but equations on assignment take places on `double` range. I think I messed up conversion, though I tried almost every possibility.

&nbsp;

And I wanted to spare my free late day. Goodbye, my points.
