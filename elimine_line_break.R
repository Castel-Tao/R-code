##this script is used for
##when you copy some sentences from pdf files,
##each row will seperated by line break
##to eliminate the line break
##I build this script
##when you use it
##you first copy sentences from pdf files
##then run this script
##and then paste sentences to your target.

### Note: Windows only
##Writer: Xiangnan Li

ori <- read.table("clipboard",header = FALSE,sep = "\n")
write.table(ori,"tmp_files.txt",quote = FALSE,col.names = FALSE,row.names = FALSE,
            eol = " ")
final <- read.table("tmp_files.txt",sep="\t",header = FALSE,stringsAsFactors = FALSE)
file.remove("tmp_files.txt")
writeClipboard(final[1,1])
