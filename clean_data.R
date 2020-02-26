
process_schedule <- function(x){
  if (ncol(x)==10) {
    x <- x[,c("V1","V3","V4","V5","V6")]
    colnames(x) = c("Date","HTeam","HScore","ATeam","AScore")
  }
  else{
    x <- x[,c("V1","V2","V3","V4","V5")]
    colnames(x) = c("Date","HTeam","HScore","ATeam","AScore")
  }
  x$Date <- as.Date(x$Date,format="%a, %b %d, %Y")
  #drop playoffs
  ids <- 1:nrow(x)
  to_keep <- ids < which(is.na(x$Date))
  return(x[to_keep,])

}

setwd("~/Documents/Github/draft-policy/data/NBA")

#37 is 1986
filenames = list.files(pattern="*.csv")
files = lapply(filenames, read.csv,header=F)

#36 is 1985
#38 is 1987
#39 is 1988
#40 us 1989

files <- lapply(files,process_schedule)

tms = c(66,67,68,69,70)
yrs = c(2015,2016,2017,2018,2019)
for (i in 1:length(tms)) {
  tm = tms[i] 
  yr = yrs[i] 
  season <- files[[tm]]
  season$hID <- as.integer(as.factor(files[[tm]]$HTeam))-1
  season$aID <- as.integer(as.factor(files[[tm]]$ATeam))-1
  season$win <- season$hID
  season$win[season$AScore>season$HScore]=season$aID[season$AScore>season$HScore]

  teams <- data.frame(team=unique(season$HTeam),id=unique(season$hID))
  teams <- teams[order(teams$id),]

  write.csv(teams,file=paste("../cleaned/teams_",yr,".csv",sep=""),row.names=F)
  write.csv(season,file=paste("../cleaned/season_",yr,".csv",sep=""),row.names=F)
} 
