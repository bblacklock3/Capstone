data_dir = [pwd '\DATA\'];
files = dir(data_dir);
for i = 1:length(files)
    bag_name = files(i).name;
    if length(bag_name) > 19
        ExtractTorque(bag_name);
    end
end