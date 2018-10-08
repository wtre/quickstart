---
title: "Homework Assignment 2"
date: 2018-10-08T19:48:31+09:00
draft: false
tags: ["Computational Photograpy"]
---


한글날을 맞아 한글로 씀.

```Matlab
clear all, addpath('src');
video_name = 'baby2';
save_memory = true;
```

먼저 `video_name`는 `face`와 `baby2` 중에서 선택하여 입력하면 된다.
영상이 좀 큰데, 특히 `baby2`는 대강 600×300×3에 900프레임이다. 엄청 큼. `save_memory`를 boolean 변수로 두어서, 필요 시 이전에 사용한 변수들을 바로바로 지우면서 메모리를 절약할 수 있게 설정했다. 


---
# Initials and color transformation

```Matlab
% Load the video file into Matlab
v = VideoReader(['data/' video_name '.mp4']);

% extract its frames
frames = zeros(v.Height, v.Width, 3, floor(v.Duration/v.Framerate), 'uint8');
framecount = 0;
while hasFrame(v)
	framecount = framecount+1;
    frames(:,:,:,framecount) = readFrame(v); %%%% frames(h,w,channel,frame)
end
```

선택한 비디오를 읽어 `frames` matrix에 저장한다. 좀 크긴 하지만 4-dim matrix를 써서 추후 벡터 연산이 용이해진다.

```Matlab
% convert them to double-precision in the range [0, 1]
frames = double(frames)/255.;

% convert each of the frames to the YIQ color space
frames_yiq = zeros(size(frames), 'double');
for k = 1:framecount
	frames_yiq(:,:,:,k) = rgb2ntsc(frames(:,:,:,k));
end
```

for문을 열어 모든 프레임을 RGB->YIQ로 변환해 준다.

---
# Laplacian pyramid

```Matlab
% call parameters from parameter_table.m
[depth, alpha, a_weight, yiq_dim, w_l, w_h, filter_order] = parameter_table(video_name);
```

피라미드 만들기에 앞서, 골랐던 비디오에 적합한 파라미터를 미리 불러온다. 우선은 depth만 쓸 거다.
파라미터의 값은 별도의 `parameter_table.m`함수에 작성해 놓았다.

```Matlab
function [depth, alpha, a_weight, yiq_dim, w_l, w_h, filter_order] = parameter_table(video_name)

switch video_name
	case 'face'
		depth = 7;
		alpha = 100;
		a_weight = [0 0 0 0 .5 1]; % except two finest(=high freq) spatio frequencies
		yiq_dim = [2 3];
		w_l = 0.83;
		w_h = 1;
		filter_order = 256;
		
	case 'baby2'
		depth = 7;
		...
```

이런 식이다.
이제 매 프레임마다 Laplacian pyramid를 만들어 보자.


```Matlab
% construct a Laplacian pyramid for every single frame in the video sequence
pyramid = cell(depth, 1); %%%% pyramid{depth}(h,w,channel,frame)
for k = framecount:-1:1 % in reverse order, for dynamical pre-allocate
	pyramid_frame = laplacian_pyramid(frames_yiq(:,:,:,k), depth);
	for i = 1:depth
		pyramid{i}(:,:,:,k) = pyramid_frame{i};
	end
end
```

각 cell의 pre-allocation을 위하여 k를 역순으로 돌렸다.
i도 역순으로 돌려서 cell 선언문을 지우고 싶었는데, Matlab에서 {1×5}와 {5×1}을 혼용하게 되어 에러가 났으므로 얌전히 초기선언함.


h×w×3 double 이미지를 받아 Laplacian pyramid를 만드는 함수는 아래와 같다.

```Matlab
function p = laplacian_pyramid(im, depth)
% Input an image, outputs a laplacian pyramid.
% im = h*w*3 (double) image array.
% depth = max depth as integer >1.
% out = {depth} cell array. each contains (h'*w'*3) images.
% by default, uses binomial filter of size 5, as in paper.

% ensure image size is devidable by multiple of 2
maxdiv = 2^(depth-1);
pad_amount = [maxdiv-1-rem(size(im,1)-1, maxdiv) maxdiv-1-rem(size(im,2)-1, maxdiv)];
if nnz(pad_amount) ~= 0
	im = padarray(im, pad_amount, 0, 'post');
end

% init a binomial filter of size 5
b = [1  4  6  4 1;
	 4 16 24 16 4;
	 6 24 36 24 6;
	 4 16 24 16 4;
	 1  4  6  4 1]/256;
 
 % init image pyramid
 p = cell(depth, 1);
 p{1} = im;
 
 % construct laplacian pyramid
 % MAKE IT RECONSTRUCTABLE!! http://www.eng.tau.ac.il/~ipapps/Slides/lecture05.pdf ->19page.
 for k = 1:depth-1
	im_filt = imfilter(p{k}, b);				% blur
	p{k+1} = im_filt(1:2:end, 1:2:end, :);		% subsample
	p{k} = p{k} - imresize(p{k+1}, 2);          % L = orig - upsample(G)
 end
```

복잡해 보이지만 내용은 맨 아래 for문이 전부다.
다만 나중에 피라미드를 다시 합쳤을 때 원본 이미지가 복원 가능하도록, 알고리즘을 수업내용과 약간 달리했다. 자세한 내용은 주석의 링크를 참조.

---
# Temporal filtering

```Matlab
% generate a Butterworth band-pass filter of a particular order
pyramid_amp = cell(depth, 1);
hd_l = butterworthBandpassFilter(v.Framerate, filter_order, w_l/2, w_h/2);
hd_r = butterworthBandpassFilter(v.Framerate, filter_order,(v.Framerate-w_h)/2, (v.Framerate-w_l)/2);
fft_hd = reshape( freqz(hd_l, framecount) + freqz(hd_r, framecount) , 1, 1, 1, [] ); % unsqueezed for vectorized .*
```

굳이 다음 챕터와 내용을 나눌 필요가 있나싶다. 제공된 `src/butterworthBandpassFilter.m`을 이용하여 필터를 생성하는 코드이다. Real 신호를 FFT했을 때 좌우대칭인 값이 나오는 점을 감안하여, filter의 cutoff frequency를 절반으로 설정하여 좌우대칭시켜 주었다. 마지막으로 벡터 곱셈에 대비하여, `reshape` 함수를 사용하여 filter를 1×1×1×f 차원의 변태적인 벡터로 변환하였다.

---
# Extracting the Frequency band of interest

```Matlab
% convert time series to the frequency domain, and apply a band pass filter to this signal
for i = 1:depth 
	pyramid_amp{i} = 0*pyramid{i}; %%%% dirty initialization
	
	iw = min(i, size(a_weight,2));
	if a_weight(iw)>0
		fft_x = fft(pyramid{i}(:,:,yiq_dim,:), [], 4);
		pyramid_amp{i}(:,:,yiq_dim,:) = real(ifft(fft_x.*fft_hd , [], 4));
	end
end
```

논문의 lambda_c값에 따른 알고리즘의 작동여부를 구현하는 대신에, 피라미드의 각 layer별 alpha의 weight를 `a_weight`라는 vector에 미리 저장해 두고 사용하였다. 논문의 구현에서, color magnification시에는 이미지 사이즈가 큰 high spatial frequncy layer에서의 alpha가 0인 경우가 많으므로, 계산량을 아끼기 위해 if문을 추가하였다. `[2 3]`값을 가지는 `yiq_dim`변수를 이용하여 I와 Q채널의 값만을 처리하는 것도 포인트.

```Matlab
% amplify it and add the result back to the original signal
pyramid_out = cell(depth, 1);
for i = 1:depth
	pyramid_out{i} = pyramid{i} + alpha*pyramid_amp{i};
end
```

Bandpass filter를 통과한 신호를 alpha배하여 원본 신호에 더해준다.

---
# Image reconstruction

```Matlab
% collapse the Laplacian pyramids into a single image per frame
pyramid_frame_out = cell(depth, 1);
for k = framecount:-1:1
	for i = 1:depth
		pyramid_frame_out{i} = pyramid_out{i}(:,:,:,k);
	end
	frames_out_yiq(:,:,:,k) = sum_pyramid(pyramid_frame_out);
end
```

피라미드를 도로 합쳤다. `sum_pyramid.m` 함수는 아래에. 사실상 세 줄짜리 함수.

```Matlab
function out = sum_pyramid(p)
% Input a pyramid generated by laplacian_pyramid.m , outpus original image.
% im = {depth}(h*w*3) (double) cell array. each contains (h'*w'*3) images.
% out = (h*w*3) image.

% reconstruct original image
for k = size(p,1)-1:-1:1
	im_upsample = imresize(p{k+1}, 2);          % upsample
	p{k} = im_upsample + p{k};					% add L
end

out = p{1};
```

이후로는 간단한 후처리. 

```Matlab
% convert each of the frames to the RGB color space
for k = framecount:-1:1
	frames_out(:,:,:,k) = ntsc2rgb(frames_out_yiq(:,:,:,k));
end

% fix range to [0 1], and trim video size if needed
frames_out = max(frames_out,0);
frames_out = min(frames_out,1);
frames_out = frames_out(1:v.Height, 1:v.Width, :, :);

% play video
implay(cat(2, frames, frames_out), v.Framerate);

% save video
v = VideoWriter(['output_' video_name '.avi']);
open(v), writeVideo(v,cat(2, frames, frames_out)), close(v);
```

이후 결과물을 보면서, 논문에서 제안한 parameter를 유지하는 선에서 spatial cutoff frequency 등등을 계속 조절했다. `face` 비디오의 경우 논문에서는 two finest layer에서 alpha=0을 사용했다고 하는데, 그대로 하면 결과물이 멍드는 듯이 noisy해서 결국은 상당히 coarse한 layer까지도 alpha=0을 두었다.  


---

결과물 보러 갑시다. 아래는 `face`. [원본 화질](https://giant.gfycat.com/AcidicSkinnyElkhound.webm) 

[![Demo CountPages alpha](https://giant.gfycat.com/AcidicSkinnyElkhound.gif)](https://giant.gfycat.com/AcidicSkinnyElkhound.webm)

Pyramid depth는 7, layer별로 alpah는 (낮은 층부터) `[0 0 0 0 .5 1 1]`로 설정하였다. 그 외에는 논문의 값과 같다. 


아래는 `baby2`. [원본 화질](https://gfycat.com/FirsthandValuableDachshund) <- 이걸로 보시는 걸 추천

[![Demo CountPages alpha](https://giant.gfycat.com/FirsthandValuableDachshund.gif)](https://gfycat.com/FirsthandValuableDachshund)

Pyramid depth는 7, layer별로 alpah는 (낮은 층부터) `[0 0 0 0 0 .5 1]`로 설정하였다. `face`와 마찬가지로 I 및 Q채널만 amplify했으며 그 외에는 논문과 같다.


