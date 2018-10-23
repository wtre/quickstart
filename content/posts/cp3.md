---
title: "Homework Assignment 3"
date: 2018-10-23T17:48:58+09:00
draft: false
tags: ["Computational Photograpy"]
---

Today's topic is image blending: connecting two different images smoothly.

&nbsp;

### 1D example of poisson blending

![image](https://i.imgur.com/wpymKJG.png)

&nbsp;

### Poisson blending on example image

![image](https://i.imgur.com/sAlx7bS.png)

&nbsp;

### Formulation 


![image](https://i.imgur.com/wfodCvp.png =50x)


* v : the (blended) result we want.
* s : source image
* t : target image(background)

On left half of the equation, we want t's gradients to be same as that of s.

The right half only counts on boundary pixels of t. (boundary can be direct outside of t, or may have some overlap. set thickness is as you want) In here we want the difference between neighboring pixels of t and v is equal to the s's gradient. In other word, we want the boundary transitions to be smooth.



We'll Manage every objective into a single least squre problem, of minimizing `(Av-b)^2`. Matlab will sove the equation for v, when we execute `v=A\b;`.

&nbsp;

&nbsp;

# Toy problem

Reconstructing an image from its gradient. (plus one pixel value). In other words, We will code left half of the equation.

Instead of reconstructing from the gradient (=1st derivative), I chose laplacian (=2nd derivative). The reason is simple: Instead of dealing with x and y gradient filter seperately (as `[1 -1]` and `[1; -1]`), I can use one symetric odd-sized filter to deal with every laplacian components.


<img src="https://i.imgur.com/nff3PyP.png" width="500">

&nbsp;

```matlab```
im = im2double(imread('data/toy_problem.png'));
f = [0 -1 0; -1 4 -1; 0 -1 0];
```

Since we should apply laplacian kernel on every pixel of image, we need a huge, sparse `A`. `A`'s #column = #pixel, and #row will be usually same or slightly larger than that. Each row will have a laplacian kernel, whose element is redistributed to match the image (which will be resized as a long single vector). I made a function named `filter2matrix.m` which returns sparse matrix `A` for given filter and mask. For example:

```matlab```
>> f = [0 -1 0;-1 4 -1;0 -1 0];
>> mask = [1 0 0; 1 0 1];
>> whole_image = [1 1 1; 1 1 1];
>> out = filter2matrix(f, mask, whole_image);
>> full(out)

ans =

     4    -1     0    -1     0     0
    -1     0    -1     4    -1     0
     0     0    -1     0    -1     4
```
4s are on the 1st, 4th and 6th column, since object `mask` is 1 on its {1, 4, 6}-th pixel. You can ignore the `whole_image` now- it's for maintaining `A`'s column size while combining different masks later.

The rest is easy.

```matlab```
% construct a sparse matrix A
mask = padarray(ones(size(im)), [1 1], 0, 'both');

A_1 = filter2matrix(f, mask);
A_3 = sparse( [1 zeros(1, numel(im)-1)] );

A = [A_1; A_3];
```
`A_3` is a single additional row for top left corner intensity.

```matlab```
% construct a known vector b
im_grad = conv2(im, f, 'same');

b_1 = reshape(im_grad', [], 1);
b_3 = im(1,1);

b = [b_1; b_3];

```

`b_1` is laplacian of the given image, and `b_3` is, you know.

To review the least squre problem, our matrix is as following:

![image](https://i.imgur.com/LdPi4ZQ.gif)

where n = #pixel.

```matlab```
% solve the least squre problem using Matlab's solver
v = A\b;

% compare the result to the original image
im_recon = reshape(v,size(im,2),[])';
figure, imshow([im im_recon],[]), title('input vs recon');
disp(['max pixel difference: '  num2str( max(max(abs(im-im_recon))) )]);
```


![image](https://i.imgur.com/8FmmxBn.png)

`max pixel difference: 5.8953e-14`


&nbsp;

&nbsp;



# Poisson blending

Let's get down to the Poisson blending.


```matlab```
im_background = imresize(im2double(imread('data/hiking.jpg')), 0.5, 'bilinear');
im_object = imresize(im2double(imread('data/penguin-chick.jpeg')), 0.5, 'bilinear');

% get source region mask from the user
objmask = getMask(im_object);
% align im_s and mask_s with im_background
[im_s, mask_s] = alignSource(im_object, objmask, im_background);
```
The interface for getting source region from user is already provided. Then when the desired position on target image is clicked, we can handle `im_s` (aligned s) and `mask_s` (aligned mask) as additional variables. 



```matlab```
f = [0 -1 0; -1 4 -1; 0 -1 0];
f_border = [0 0 0; 0 4 0; 0 0 0];
se = strel('disk', 2);
se_border = strel('disk', 2);

% construct a sparse matrix A
A_1 = filter2matrix(f, imerode(objmask, se), objmask);

bordermask = objmask - imerode(objmask, se_border);
A_2 = filter2matrix(f_border, bordermask, objmask);

A = [A_1; A_2];
A = blkdiag(A, A, A);
```
For given mask, I eroded the outermost pixels in advance. I also set the outermost pixels as `bordermask`. By adjusting the radius of structure element(`strel`), I can make overlap between `bordermask` and the eroded `objmask`.

Computing `A` takes most of the time(10s to minutes depending on the size of `im_object`). I copied it twice and constructed a block diagonal matrix, to compute {R, G, B} components in one least square problem.


```matlab```
% construct a known vector b
bordermask_s = mask_s - imerode(mask_s, se_border);
b = zeros(0,1);
for c = 1:3
	% inner part
	im_bg1 = im_background(:,:,c);
	im_s1 = im_s(:,:,c);

	im_grad = conv2(im_s1, f, 'same');
	im_grad_reshaped = reshape(im_grad', [], 1);
	mask_reshaped = reshape(imerode(mask_s, se)', [], 1);

	b_1 = im_grad_reshaped;
	b_1(boolean(~mask_reshaped)) = [];

	% boundary part
	im_b_grad = circshift(im_bg1, [1 0]) + circshift(im_bg1, [-1 0]) ...
		+ circshift(im_bg1, [0 1]) + circshift(im_bg1, [0 -1]);

	im_b_grad_reshaped = reshape(im_b_grad', [], 1);
	bordermask_reshaped = reshape(bordermask_s', [], 1);

	b_2 = im_grad_reshaped + im_b_grad_reshaped;
	b_2(boolean(~bordermask_reshaped)) = [];

	b = [b; b_1; b_2];
end
```

I then constructed a known vector b. For single channel the matrix is as follows:

![image](https://i.imgur.com/1tm5sdT.gif)

where `b_i` denotes i-th boundary pixel and the diamond symbol stands for summation on the 4 neighboring pixels of `b_i`. The rest is straightforward.


```matlab```
% solve the least squre problem using Matlab's solver
v = A\b;


% fill the mask image with desired value
for c = 1:3
	vs = size(v,1);

	single_channel = double(mask_s)';
	single_channel(single_channel>=1) = v(vs*(c-1)/3+1:vs*c/3);
	im_out(:,:,c) = mask_s.*single_channel' + (1-mask_s).*im_background(:,:,c);
end

figure, imshow(im_out,[]);
```

Note that `A` is exactly a squre matrix, when there are no overlap between object mask and boundary mask. If you think this is 'unstable', you can adjust `strel` size (which kinda helps in case of mixed blending).


# Mixed blending

Everything is same as the Poisson blending, except a small modification:

![image](https://i.imgur.com/7skWsS7.png)

where

![image](https://i.imgur.com/kTSfQBf.png)

In codewise a couple of lines are added, to compute background gradient, and compare its absolute values with those of the target image. 

```matlab```
for c = 1:3
	% inner part
	im_bg1 = im_background(:,:,c);
	im_s1 = im_s(:,:,c);

	im_grad = conv2(im_s1, f, 'same');
	im_grad_reshaped = reshape(im_grad', [], 1);

	bg_grad = conv2(im_bg1, f, 'same');
	bg_grad_reshaped = reshape(bg_grad', [], 1);

	gc_mask = (abs(im_grad_reshaped) > abs(bg_grad_reshaped));
	larger_grad_reshaped = gc_mask .* im_grad_reshaped + (1-gc_mask) .* bg_grad_reshaped;

	mask_reshaped = reshape(imerode(mask_s, se)', [], 1);

	b_1 = larger_grad_reshaped;
	b_1(boolean(~mask_reshaped)) = [];

	% boundary part
	im_b_grad = circshift(im_bg1, [1 0]) + circshift(im_bg1, [-1 0]) ...
		+ circshift(im_bg1, [0 1]) + circshift(im_bg1, [0 -1]);

	im_b_grad_reshaped = reshape(im_b_grad', [], 1);
	bordermask_reshaped = reshape(bordermask_s', [], 1);

	b_2 = larger_grad_reshaped + im_b_grad_reshaped;
	b_2(boolean(~bordermask_reshaped)) = [];

	b = [b; b_1; b_2];
end
```


And now's the fun part.

&nbsp;

&nbsp;

# Results on provided image


### Poisson blending

<img src="https://i.imgur.com/XCbBaKd.jpg" width="340"> <img src="https://i.imgur.com/kFVYDwY.jpg" width="200"> <img src="https://i.imgur.com/C7U4Am6.png" width="200">

From left to right. Background image (target), source image and masked image.

&nbsp;

![image](https://i.imgur.com/dDFQe2H.jpg)

* lower left: Optimal case.
* upper right: Target border is blue, source border is white. By connecting them while maintaining source gradient, the whole image is sent to the blue-ish side. 
* upper left: Will explain later.

&nbsp;

### Mixed gradient

![image](https://i.imgur.com/TtbB2Xb.jpg)

* upper left: The tree in background has stronger gradient (and laplacian) then the belly of penguin. Thereby the tree is selected.
* upper right: When the background has near-0 gradient, mixed gradient blending works as Poisson blending.
* lower left: The mediocre gradient of snowy surface brought some unstable results.

&nbsp;

&nbsp;

# Own examples

### case 1

<img src="https://i.imgur.com/jg0304a.jpg" width="500">

Target image.

<img src="https://i.imgur.com/bqHbXMF.png" width="250"> <img src="https://i.imgur.com/3DdHHz2.png" width="500">

Sorce image, and aligned. Both are AA-based squishy AP carry, so basically the same champ, right?

&nbsp;

![image](https://i.imgur.com/IZeqQqW.png)

Poisson proves it.

&nbsp;

![image](https://i.imgur.com/UTzEwr9.png)

Mixed gradient is not expected to work here, since Teemo's face should be priortized, as always.

&nbsp;


### case 2

<img src="https://i.imgur.com/id2DA0g.jpg" width="500">

my desk.

&nbsp;

<img src="https://i.imgur.com/1h5SDml.png" width="380"><img src="https://i.imgur.com/1BcPUrx.jpg" width="380">

Poisson blending some sketches on the notebook. Although the quality is good, parallel lines on the notebook are covered. 

&nbsp;

<img src="https://i.imgur.com/Fu0mPBr.png" width="380"><img src="https://i.imgur.com/7TFTm4Q.jpg" width="380">

Mixed gradient blending reveals the line of notebook, and expected to be more realistic. However the results are kind of unstable. I assume this is because I compared the absolute of laplacian (instead of gradient), however this can be a code bug.

&nbsp;

### case 3

<img src="https://i.imgur.com/4o7GgH1.jpg" width="480"> <img src="https://i.imgur.com/1MSvqlS.jpg" width="270">

Target and source image.

&nbsp;

![image](https://i.imgur.com/9TCk7bp.png)

![image](https://i.imgur.com/xXtNlSG.jpg)

Poisson blended result. If you pay attention, you can notice the floor mat texture around the cat. 

&nbsp;
