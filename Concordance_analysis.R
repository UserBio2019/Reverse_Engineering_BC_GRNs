Concordance_analysis <- function(fileName, 
                       check.names=FALSE, header = FALSE,
                       stringsAsFactor = FALSE, 
                       sep = "\t",...){

Data <- read.delim(fileName, 
                      header = header,
                     check.names = check.names, 
                     stringsAsFactor = stringsAsFactor, sep = sep, ...)

   tmp2 = read.delim("Metabric_centered.txt", header = T)
   mtch = match(tmp2[,1], Data[,1])
   mtch = as.data.frame(mtch)
      cb = cbind(mtch, tmp2)
      sb = cb[which(cb[,1] != "NA"),]
      sb = sb[,-1]
      c = Data[[1]][1]
         mtch = match(sb[,1], c)
         mtch = as.data.frame(mtch)
         cb = cbind(mtch, sb)
         cb <- cb[order(cb[,1], decreasing = TRUE),  ]
         cb = cb[,-1]
         row.names(cb) = cb[,1]
         cb = cb[,-1]
         cb = t(cb)
            cb[cb>0] <- 1
            cb[cb<0] <- -1
            cb <- cb[order(cb[,1], decreasing = TRUE),  ]
               newdata_a <- cb[which(cb[,1]==1), ]
               newdata_b <- cb[ which(cb[,1]==-1), ]
               newdata_a = as.data.frame(t(newdata_a))
               newdata_b = as.data.frame(t(newdata_b))
               newdata_a_t <- table(stack(newdata_a))
               newdata_b_t <- table(stack(newdata_b))
               newdata_a_t <- as.data.frame(newdata_a_t)
               newdata_b_t <- as.data.frame(newdata_b_t)
               newdata_a_t= reshape(newdata_a_t, timevar = "ind", idvar = "values", direction = "wide")
               newdata_b_t= reshape(newdata_b_t, timevar = "ind", idvar = "values", direction = "wide")
        df_a = newdata_a_t[,-1]
        df_b = newdata_b_t[,-1]
        df_a_freq= df_a[apply(df_a[ncol(df_a)],1,function(z) !any(z==0)),] 
        df_b_freq= df_b[apply(df_b[ncol(df_b)],1,function(z) !any(z==0)),] 
        x = df_a_freq[2,] - df_a_freq[1,]
        row.names(x) = c 
        Hub_name = c
        cb_x = cbind(Hub_name, x)
        write.table(cb_x, paste(c, "_pos_conc_count.txt",sep=""),  col.names = T, row.names = F, sep = "\t")
      y = df_b_freq[1,] - df_b_freq[2,]
      row.names(y) = c
      Hub_name = c
      cb_y = cbind(Hub_name, y)
      write.table(cb_y, paste(c, "_neg_conc_count.txt", sep = ""), col.names = T, row.names = F, sep = "\t")
   cb = cbind(x, y)
   row.names(y) = c
   Hub_name = c
   cb_x_y=cbind(Hub_name, cb)
   write.table(cb_x_y, paste(c, "_cbind_pos_conc_neg_conc_not_norm.txt", sep = ""), col.names = T, row.names = F, sep = "\t")
x = x / nrow(newdata_a)
row.names(x) = c 
Hub_name = c
cb_x_norm = cbind(Hub_name, x)
write.table(cb_x_norm, paste(c,  "_pos_conc_count_norm.txt", sep = ""), col.names = T, row.names = F, sep = "\t")
y = y / nrow(newdata_b)
row.names(y) = c 
Hub_name = c
cb_y_norm = cbind(Hub_name, y)
write.table(cb_y_norm, paste(c,  "_neg_conc_count_norm.txt", sep = ""),  col.names = T, row.names = F, sep = "\t")
cb = cbind(x, y)
row.names(cb) = c 
Hub_name = c
cb_x_y_norm = cbind(Hub_name, cb)
write.table(cb_x_y_norm, paste(c, "_cbind_pos_conc_neg_conc_norm.txt", sep = ""),  col.names = T, row.names = F, sep = "\t")
return(cb_x_y_norm)
#return(tryCatch((z=x+y), error=function(e) NULL))
}
