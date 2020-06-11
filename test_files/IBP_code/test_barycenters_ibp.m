% Code modified from Solomon et al.'s other code.

function test_barycenters_ibp(n,k,eps_string,test_type,filterSize,img_or_timing)
    % 0.eps_string is the epsilon.
    % img_or_timing = 1 --> save image at each timestep
    % img_or_timing = 2 --> do not save image at each timestep -- save timing info instead
    % times instead
    path(path,'blur_functions/');
    path(path,'convolutional_wasserstein/');
    path(path,'image_blur/');
    path(path,'toolbox/');

    %% Read data

    directory = ['../experiment_results/ibp_results/', test_type, '_n', num2str(n), 'k', num2str(k), 'eps', eps_string, 'f', num2str(filterSize)]
    mkdir(directory)
    file_prefix = directory;

    pims = cell(1,k);
    ps = cell(1,k);
    for i = 1:k
        pims{i} = ['../experiment_data/', test_type, '_n', num2str(n), 'k', num2str(k), 'eps', eps_string, '_', num2str(i),'.png'];
        ps{i} = im2double(imread(pims{i}));
        ps{i} = ps{i} / sum(sum(ps{i}));
    end

    imSize = size(ps{1})
    for i = 1:k
        ps{i} = ps{i}(:);
        ps{i} = ps{i} + 1e-8;
        ps{i} = ps{i} / sum(ps{i});
    end

    n = length(ps{1});
    p = cat(2, ps{:})*n;
    areaWeights = ones(n,1)/n;

    %% Set up blur

    h = fspecial('gaussian',[1 max(imSize)],filterSize);% hsize sigma
    h = h / sum(h);
    imBlur = @(x) imfilter(imfilter(x,h,'replicate'),h','replicate');
    blurColumn = @(x) reshape(imBlur(reshape(x,imSize)),[],1);
    blurAll = @(x) cell2mat(cellfun(blurColumn, num2cell(x,1), 'UniformOutput', false));

    %% Compute barycenter


    entropies = -sum(p.*log(p).*repmat(areaWeights, [1,k]));
    minEntropy = min(entropies);
    targetEntropy = minEntropy/100;

    close all;
    alpha = repmat(1/k, [1,k])



    options.verb = 1;
    options.img_or_timing = img_or_timing
    options.unit_area_projection = 0;
    options.imSize = imSize;
    options.file_prefix = file_prefix;

    [~,~] = convolutionalBarycenter(p,alpha,areaWeights,blurAll,blurAll,targetEntropy, options);

end

