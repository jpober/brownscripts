function oscar_folder_locs, folder_names_in
  
  folder_names = folder_names_in
  for i=0, n_elements(folder_names)-1 do begin
    folder_test = file_test(folder_names[i], /directory)
    if folder_test eq 0 then begin
        start_path = '/users/wl42/scratch/FHD_out/'
        test_name = start_path + folder_names[i]
        folder_test = file_test(test_name, /directory)
        if folder_test eq 1 then folder_names[i] = test_name
    endif
    if folder_test eq 0 then begin
        start_path = '/users/wl42/data/wl42/FHD_out/'
        test_name = start_path + folder_names[i]
        folder_test = file_test(test_name, /directory)
        if folder_test eq 1 then folder_names[i] = test_name
    endif
    if folder_test eq 0 then message, 'folder not found'
  endfor
  return, folder_names
end
