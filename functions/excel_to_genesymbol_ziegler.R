# Convert excel genes to gene symbol for Ziegler dataset

excel_to_genesymbol=function(data){
  rown=rownames(data)
rown=gsub("1-Dec","DELEC1",rown,fixed = T)
rown=gsub("11-Mar","MARCHF11",rown,fixed = T)
rown=gsub("10-Mar","MARCHF10",rown,fixed = T)
rown=gsub("2-Mar","MARCHF2",rown,fixed = T)
rown=gsub("3-Mar","MARCHF3",rown,fixed = T)
rown=gsub("4-Mar","MARCHF4",rown,fixed = T)
rown=gsub("5-Mar","MARCHF5",rown,fixed = T)
rown=gsub("6-Mar","MARCHF6",rown,fixed = T)
rown=gsub("7-Mar","MARCHF7",rown,fixed = T)
rown=gsub("8-Mar","MARCHF8",rown,fixed = T)
rown=gsub("9-Mar","MARCHF9",rown,fixed = T)
rown=gsub("1-Mar","MARCHF1",rown,fixed = T)

rown=gsub("MARCHF1_1","MTARC1",rown,fixed = T)
rown=gsub("MARCHF2_1","MTARC2",rown,fixed = T)

rown=gsub("12-Sep","SEPTIN12",rown,fixed = T)
rown=gsub("2-Sep","SEPTIN2",rown,fixed = T)
rown=gsub("9-Sep","SEPTIN9",rown,fixed = T)
rown=gsub("11-Sep","SEPTIN11",rown,fixed = T)
rown=gsub("10-Sep","SEPTIN10",rown,fixed = T)
rown=gsub("14-Sep","SEPTIN14",rown,fixed = T)
rown=gsub("15-Sep","SELENOF",rown,fixed = T)

rown=gsub("1-Sep","SEPTIN1",rown,fixed = T)
rown=gsub("3-Sep","SEPTIN3",rown,fixed = T)
rown=gsub("4-Sep","SEPTIN4",rown,fixed = T)
rown=gsub("5-Sep","SEPTIN5",rown,fixed = T)
rown=gsub("6-Sep","SEPTIN6",rown,fixed = T)
rown=gsub("7-Sep","SEPTIN7",rown,fixed = T)
rown=gsub("8-Sep","SEPTIN8",rown,fixed = T)
rown

}
