function db = loaddata( numOfSamples, s_modalities, d_modalities, filename)
%
fprintf('Loading data ... ');
tic;
fp=fopen(filename);

for i=1:s_modalities %默认是生成到db{2}
    db{i}.stride = [];
    db{i}.w = [];
    db{i}.supp = [];
end

count = 0;
while ~feof(fp) % 按行读入，直到文件末尾
  count = count +1;
  for i=1:s_modalities      
    fscanf(fp, '%d', 1);
    [d check] = fscanf(fp, '%d', 1); %这里的check是读取成功的字符数
    if check == 0 break; end

    db{i}.stride(end+1) = d;
    we = fscanf(fp, '%f', [1, d]);%读取到[1,d]的数组中
    db{i}.w(1,(end+1):(end+d)) = we/sum(we);
    db{i}.supp(:, (end+1):(end+d)) = fscanf(fp, '%f', [d_modalities(i), d]);      
  end
  if (count == numOfSamples) break; end
end

fclose(fp);
toc;
fprintf('[done]');
end

