function x_out = fun_padarray_2D(x, pad_dim, val)

% pad_dim(1,2,3,4) = (top, bottom, left, right)
% Example:
% x =
% 
%      1     2
%      3     4
% x_out = fun_padarray_2D(x, [2 3 1 4], [100 101 0 50])
% 
% x_out =
% 
%      0   100   100    50    50    50    50
%      0   100   100    50    50    50    50
%      0     1     2    50    50    50    50
%      0     3     4    50    50    50    50
%      0   101   101    50    50    50    50
%      0   101   101    50    50    50    50
%      0   101   101    50    50    50    50

x_12 = [ones([pad_dim(1), size(x,2)])*val(1); x; ones([pad_dim(2), size(x,2)])*val(2);];

x_temp = x_12';
x_1234 = [ones([pad_dim(3), size(x_temp,2)])*val(3);  x_temp; ones([pad_dim(4), size(x_temp,2)])*val(4);];

x_out = x_1234';


% x_out = [x_out; ones([pad_dim(1), size(x,2)])*val(1)];
% x_out = [x_out; ones([pad_dim(1), size(x,2)])*val(1)];
% 
