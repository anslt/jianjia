#Read the data from file(You may change the path)
data_raw<-read.csv("D:/stat/wines.csv")

#Pick the useful part
data<-data_raw[,2:54]

#Make a list store the block division
pick1<-list(seq(1,6),seq(7,12),seq(13,18),seq(19,23),seq(24,29),
            seq(30,34),seq(35,38),seq(39,44),seq(45,49),seq(50,53))

col_name<-names(data)

pick2<-list(col_name[1:6],col_name[7:12],col_name[13:18],col_name[19:23],col_name[24:29],
            col_name[30:34],col_name[35:38],col_name[39:44],col_name[45:49],col_name[50:53])

#Construct the object 
review<-mfa1(data,pick1,scale=1)

#Test the each meathods and functions
#1.summary_eigen
test1<-summary_eigen(review)

#2.contribution
test2<-contribution(review,type="observation")
test3<-contribution(review,type="variable")
test4<-contribution(review,type="table")

#3.RV
num_block<-length(review$weight)
test5<-RV(review$scale_x[,seq(review$divide[1],review$divide[2]-1)],
          review$scale_x[,seq(review$divide[num_block-1],review$divide[num_block]-1)])

#4.LG
test6<-LG(review$scale_x[,seq(review$divide[1],review$divide[2]-1)],
          review$scale_x[,seq(review$divide[num_block-1],review$divide[num_block]-1)])