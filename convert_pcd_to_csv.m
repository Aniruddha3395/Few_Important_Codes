% Author    : Aniruddha Shembekar, Research Engineer, University of
% Southern California

function pts = convert_pcd_to_csv(file_name)

pcdata = pcread(file_name);
pcwrite(pcdata,'temp_file_no_duplicate_caption.ply');
pcdata = pcread('temp_file_no_duplicate_caption.ply');

delete('temp_file_no_duplicate_caption.ply');
pts = zeros(pcdata.Count,3);

for i=1:pcdata.Count
    pts(i,:) = pcdata.Location(i,1:3);
end

end

