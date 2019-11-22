
setwd("~/Documents/Github/draft-policy/data/NBA")

filenames = list.files(pattern="*.csv")
files = lapply(filenames, read.csv,header=F)

files <- lapply(files,process_schedule)
season <- files[[37]]
season$hID <- as.integer(as.factor(files[[37]]$HTeam))-1
season$aID <- as.integer(as.factor(files[[37]]$ATeam))-1
season$win <- season$hID
season$win[season$AScore>season$HScore]=season$aID[season$AScore>season$HScore]

teams <- data.frame(team=unique(season$HTeam),id=unique(season$hID))
teams <- teams[order(teams$id),]

write.csv(teams,file="teams_1986.csv",row.names=F)
write.csv(season,file="season_1986.csv",row.names=F)


teams <- as.factor(files[[37]]$HTeam)
records_86 <- transform_schedule(files[[37]])
results <- calculate_bound(records_19) 
loss_pct <- results[[2]]
eliminated <- results[[3]] 
wins <- results[[4]]

team_scores <- rep(0,length(unique(records_19$Team))) 
for (i in 1:length(unique(records_19$Team))) { 
  print(records_19$Team[i]) 
  team_gone <- which(eliminated[,i])[1] 
  print(team_gone)
  team_scores[i] <- loss_pct[team_gone,i]
}

draft_results <- data.frame(season_wins=wins[nrow(wins),],scores=team_scores)
rownames(draft_results) <- unique(records_19$Team)
draft_results <- draft_results[!is.na(draft_results$scores),]
draft_results$probs <- draft_results$scores/sum(draft_results$scores)

p_change <- results[[1]]
p_change_knicks <- p_change[,"New York Knicks"]
p_change_raps <- p_change[,"Toronto Raptors"]
p_change_nets <- p_change[,"Brooklyn Nets"]
matplot(1:length(p_change_nets),data.frame(knicks=p_change_knicks,raps=p_change_nets),type="l") 

team_wins<- function(team,date,season){
  records = season$Team==team & season$Date < date 
  return(sum(season$Win[records])) 
}

games_played <- function(team,date,season){
  records = season$Team==team & season$Date < date 
  return(sum(records)) 
}
# calculate_bound <- function(data,n.games=82) { 
#   dates <- unique(data$Date) 
#   teams <- unique(data$Team)
#   
#   wins <- matrix(0,nrow=length(dates),ncol=length(teams))
#   gp <- matrix(0,nrow=length(dates),ncol=length(teams)) 
#   loss_pct <- matrix(0,nrow=length(dates),ncol=length(teams))
#   can_lose <- matrix(n.games,nrow=length(dates),ncol=length(teams)) 
#   can_lose_w <- matrix(n.games,nrow=length(dates),ncol=length(teams)) 
#   for (d in 2:nrow(wins)) {
#       for (t in 1:ncol(wins)){
#         wins[d,t] <- team_wins(teams[t],dates[d],data) 
#         gp[d,t] <- games_played(teams[t],dates[d],data)
#         loss_pct[d,t] <- (gp[d,t] - wins[d,t])/gp[d,t] 
#       }
#     for (t in 1:ncol(wins)){
#       can_lose[d,t] <- max(n.games - gp[d,t]+wins[d,t] - max(wins[d,]),0) 
#       can_lose_w[d,t] <- max(n.games - gp[d,t]+wins[d,t] +1  - max(wins[d,]),0) 
#     }
#     
#     #for (t in 1:ncol(wins)) { 
#      # can_lose[d,t] <- max(n.games - gp[d,t]+wins[d,t] - max(wins[d,]),0)/max(can_lose[d,])
#      # can_lose_w[d,t] <- max(n.games - gp[d,t]+wins[d,t] +1  - max(wins[d,]),0)/max(can_lose[d,]) 
#     #}
#   }
#   
#   eliminated = (can_lose==0)
#   p_win_l <- can_lose/rowSums(can_lose)
#   p_win_w <- (can_lose_w)/(rowSums(can_lose)) 
#   
#   result <- p_win_w - p_win_l 
#   #result <- can_lose_w - can_lose 
#   rownames(result) = dates
#   colnames(result) = teams 
#   return(list(result,loss_pct,eliminated,wins))  
# }

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

transform_schedule <- function(x) {
  df <- rbind(t(apply(x,MARGIN=1,FUN=home_entry)), t(apply(x,MARGIN=1,FUN=away_entry))) 
  colnames(df) = c("Date","Team","Win") 
  df <- as.data.frame(df)
  df <- df[order(df$Date),]
  df$Date <- as.Date(df$Date)
  df$Win <- as.numeric(df$Win)-1 
  return(df)  
}

away_entry <- function(game) { 
  row2 <- c(game["Date"],game["ATeam"],as.numeric(game["AScore"]>game["HScore"])) 
  return(row2) 
}
home_entry <- function(game) { 
  row1 <- c(game["Date"],game["HTeam"],as.numeric(game["HScore"]>game["AScore"]))
  return(row1)
}