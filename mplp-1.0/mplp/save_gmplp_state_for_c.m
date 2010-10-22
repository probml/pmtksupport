function save_gmplp_state_for_c(gmplp_state,base_dir)

regions_fn = [base_dir '/regions.txt'];
lambdas_fn = [base_dir '/lambdas.txt'];
intersects_fn = [base_dir '/intersects.txt'];
var_sizes_fn = [base_dir '/var_sizes.txt'];
region_intersects_fn = [base_dir '/region_intersects.txt'];

% Write regions file
fid = fopen(regions_fn,'wt');
for i=1:length(gmplp_state.regions)
    fprintf(fid,'%d ',gmplp_state.regions{i});
    fprintf(fid,'\n');
end
fclose(fid);

fid = fopen(lambdas_fn,'wt');
for i=1:length(gmplp_state.regions)
    tmp = gmplp_state.lambda{i};
    if isempty(tmp)
        fprintf(fid,'\n');
    else    
        fprintf(fid,'%g ',tmp);
        fprintf(fid,'\n');
    end
end
fclose(fid);

fid = fopen(intersects_fn,'wt');
for i=1:length(gmplp_state.intersects)
    fprintf(fid,'%d ',gmplp_state.intersects{i});
    fprintf(fid,'\n');
end
fclose(fid);

fid = fopen(var_sizes_fn,'wt');
fprintf(fid,'%d\n',gmplp_state.var_size);
fclose(fid);


fid = fopen(region_intersects_fn,'wt');
for i=1:length(gmplp_state.region_subsets)
    fprintf(fid,'%d ',gmplp_state.region_subsets{i});
    fprintf(fid,'\n');
end
fclose(fid);

fid = fopen('msgs.txt','wt');
for ri=1:length(gmplp_state.regions)
  for si=1:length(gmplp_state.region_subsets{ri})
    fprintf(fid,'%g ',gmplp_state.curr_msgs{ri,si});
  end
end
fclose(fid);
