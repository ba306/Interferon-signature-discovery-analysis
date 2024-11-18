# Rename Aybey signature and remove IFNb Aybey

rosetta_new=readRDS(file ="published_signatures/IFN_signature_list_all.RDS")

rosetta_new_valid=rosetta_new
# After validation step
# rename IFNa_Aybey and IFNg_Aybey
names(rosetta_new_valid)[grep("IFNa_Aybey",names(rosetta_new_valid))]=
  "IFN_I_Aybey"
names(rosetta_new_valid)[grep("IFNg_Aybey",names(rosetta_new_valid))]=
  "IFN_II_Aybey"
# rename IFNb
rosetta_new_valid=rosetta_new_valid[!grepl("IFNb_Aybey",names(rosetta_new_valid))]

saveRDS(rosetta_new_valid,file ="published_signatures/IFN_signature_list_all_valid.RDS")

