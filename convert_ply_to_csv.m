%****************************************************************************************
%
% Author : Aniruddha Shembekar, University of Southern California
%
%****************************************************************************************

function pts = convert_ply_to_csv(file_name)

pcdata = pcread(file_name);
pts = zeros(pcdata.Count,3);

for i=1:pcdata.Count
    pts(i,:) = pcdata.Location(i,1:3);
end

end