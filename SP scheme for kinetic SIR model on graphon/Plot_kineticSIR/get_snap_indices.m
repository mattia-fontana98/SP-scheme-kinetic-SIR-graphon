function snap = get_snap_indices(time, t_snap)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Function finding indices in time closest to requested snapshot times    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


t_snap = t_snap(:).';

snap = zeros(size(t_snap));

for s = 1:numel(t_snap)

    [~, snap(s)] = min(abs(time - t_snap(s)));

end


end