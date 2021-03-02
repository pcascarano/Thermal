function upU=my_upsampling_bicubic(g_image,rows, columns, up_factor)
   upU=imresize(g_image,up_factor,'Bicubic');
end