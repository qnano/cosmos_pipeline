function img = glimps_readall(folder,vid)


    for i=1:size(vid.field,2)
        a=glimpse_image(folder,vid,i);
        if ~isempty(a)
            img(:,:,i) = a;
        end
    end


end