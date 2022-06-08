target_inis = {input_inis};
cd {airlocalize_root}
for i=1:length(target_inis)
    target_ini = convertStringsToChars(target_inis(i)); % you have to provide character vector and NOT string to AIRLOCALIZE!
    AIRLOCALIZE(target_ini)
end