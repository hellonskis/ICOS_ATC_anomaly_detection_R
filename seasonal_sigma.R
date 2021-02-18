#############################################################################
## Calculate sigma curves using a 90-day seasonal window
#############################################################################

seasonal_sigma <- function(data, ndays) {

  sd <- setNames(data.frame(matrix(ncol = 1, nrow = nrow(data))), c("sd")) 
  r <- setNames(data.frame(matrix(ncol = 1, nrow = 0)), c("diff")) 
  r.rs <- setNames(data.frame(matrix(ncol = 1, nrow = 1)), c("diff"))   
  assign.sd <- FALSE
  
  for (i in (nrow(data)-ndays):nrow(data)) {
    print(i)
      if (i == (nrow(data)-ndays)) {
        month <- month(data$sampling_datetime[i])
        day <- day(data$sampling_datetime[i])
      }
      if ((month(data$sampling_datetime[i]) == month) & (day(data$sampling_datetime[i]) == day) & (i != (nrow(data)-ndays))) {
        assign.sd <- TRUE
        print ("1 year complete, assigning sigma values")
        for (j in i:nrow(data)) {
          sd$sd[j] <- sd$sd[j-365]
        }
      }
      if (isTRUE(assign.sd)) {
        print("done")
        break
      }
      if ((yday(data$sampling_datetime[i]) >= 46) & (yday(data$sampling_datetime[i]) <= 320)) { 
        for (j in 1:nrow(data)) {
          if ((yday(data$sampling_datetime[j]) >= yday(data$sampling_datetime[i]) - 45) &
              (yday(data$sampling_datetime[j]) <= yday(data$sampling_datetime[i]) + 45)) {
            r.rs$diff[1] <- data$dC[j]
            r <- rbind(r, r.rs)
          }
        }
        sd$sd[i] <- sd(r$diff)
        r <- setNames(data.frame(matrix(ncol = 1, nrow = 0)), c("diff")) 
        r.rs <- setNames(data.frame(matrix(ncol = 1, nrow = 1)), c("diff"))
      } else if (yday(data$sampling_datetime[i]) < 46) {
        for (j in 1:nrow(data)) {
          if ((yday(data$sampling_datetime[j]) <= yday(data$sampling_datetime[i])) | 
              (yday(data$sampling_datetime[j]) >= ((yday(data$sampling_datetime[i]) - 45) + yearDays(data$sampling_datetime[i]-(86400*46)))) |
              (yday(data$sampling_datetime[j]) <= yday(data$sampling_datetime[i]) + 45)) {
            r.rs$diff[1] <- data$dC[j]
            r <- rbind(r, r.rs)
          }
        }
        sd$sd[i] <- sd(r$diff)
        r <- setNames(data.frame(matrix(ncol = 1, nrow = 0)), c("diff")) 
        r.rs <- setNames(data.frame(matrix(ncol = 1, nrow = 1)), c("diff"))
      } else if (yday(data$sampling_datetime[i]) > 320) {
        for (j in 1:nrow(data)) {
          if ((yday(data$sampling_datetime[j]) >= yday(data$sampling_datetime[i])) | 
              (yday(data$sampling_datetime[j]) <= ((yday(data$sampling_datetime[i]) + 45) - 365)) |
              (yday(data$sampling_datetime[j]) >= yday(data$sampling_datetime[i]) - 45)) {
            r.rs$diff[1] <- data$dC[j]
            r <- rbind(r, r.rs)
          }
        }
        sd$sd[i] <- sd(r$diff)
        r <- setNames(data.frame(matrix(ncol = 1, nrow = 0)), c("diff")) 
        r.rs <- setNames(data.frame(matrix(ncol = 1, nrow = 1)), c("diff"))
      }
  }
  
  return(sd)
  
}
  

###############################################################################################


seasonal_sigma.nr <- function(data, ndays) {
  
  sd <- setNames(data.frame(matrix(ncol = 1, nrow = nrow(data))), c("sd")) 
  r <- setNames(data.frame(matrix(ncol = 1, nrow = 0)), c("diff")) 
  r.rs <- setNames(data.frame(matrix(ncol = 1, nrow = 1)), c("diff"))   
  
  for (i in (nrow(data)-ndays):nrow(data)) {
    print(i)
    if ((yday(data$sampling_datetime[i]) >= 46) & (yday(data$sampling_datetime[i]) <= 320)) { 
      for (j in 1:nrow(data)) {
        if ((yday(data$sampling_datetime[j]) >= yday(data$sampling_datetime[i]) - 45) &
            (yday(data$sampling_datetime[j]) <= yday(data$sampling_datetime[i]) + 45)) {
          r.rs$diff[1] <- data$dC[j]
          r <- rbind(r, r.rs)
        }
      }
      nr <- setNames(data.frame(matrix(ncol = 2, nrow = sum(r$diff < 0))), c("n", "p"))
      nr$n <- subset(r$diff, (r$diff < 0))
      for (k in 1:nrow(nr)) {
        nr$p[k] <- nr$n[k] * (-1)
      }
      normalized <- data.frame(res = c(nr[,"n"], nr[,"p"]))
      sd$sd[i] <- sd(normalized$res)
      # sd$sd[i] <- sd(normalized$res)/2
      r <- setNames(data.frame(matrix(ncol = 1, nrow = 0)), c("diff")) 
      r.rs <- setNames(data.frame(matrix(ncol = 1, nrow = 1)), c("diff"))
    } else if (yday(data$sampling_datetime[i]) < 46) {
      for (j in 1:nrow(data)) {
        if ((yday(data$sampling_datetime[j]) <= yday(data$sampling_datetime[i])) | 
            (yday(data$sampling_datetime[j]) >= ((yday(data$sampling_datetime[i]) - 45) + 365)) |
            (yday(data$sampling_datetime[j]) <= yday(data$sampling_datetime[i]) + 45)) {
          r.rs$diff[1] <- data$dC[j]
          r <- rbind(r, r.rs)
        }
      }
      nr <- setNames(data.frame(matrix(ncol = 2, nrow = sum(r$diff < 0))), c("n", "p"))
      nr$n <- subset(r$diff, (r$diff < 0))
      for (k in 1:nrow(nr)) {
        nr$p[k] <- nr$n[k] * (-1)
      }
      normalized <- data.frame(res = c(nr[,"n"], nr[,"p"]))
      sd$sd[i] <- sd(normalized$res)
      # sd$sd[i] <- sd(normalized$res)/2
      r <- setNames(data.frame(matrix(ncol = 1, nrow = 0)), c("diff")) 
      r.rs <- setNames(data.frame(matrix(ncol = 1, nrow = 1)), c("diff"))
      
    } else if (yday(data$sampling_datetime[i]) > 320) {
      for (j in 1:nrow(data)) {
        if ((yday(data$sampling_datetime[j]) >= yday(data$sampling_datetime[i])) | 
            (yday(data$sampling_datetime[j]) <= ((yday(data$sampling_datetime[i]) + 45) - 365)) |
            (yday(data$sampling_datetime[j]) >= yday(data$sampling_datetime[i]) - 45)) {
          r.rs$diff[1] <- data$dC[j]
          r <- rbind(r, r.rs)
        }
      }
      nr <- setNames(data.frame(matrix(ncol = 2, nrow = sum(r$diff < 0))), c("n", "p"))
      nr$n <- subset(r$diff, (r$diff < 0))
      for (k in 1:nrow(nr)) {
        nr$p[k] <- nr$n[k] * (-1)
      }
      normalized <- data.frame(res = c(nr[,"n"], nr[,"p"]))
      sd$sd[i] <- sd(normalized$res)
      # sd$sd[i] <- sd(normalized$res)/2
      r <- setNames(data.frame(matrix(ncol = 1, nrow = 0)), c("diff")) 
      r.rs <- setNames(data.frame(matrix(ncol = 1, nrow = 1)), c("diff"))
    }
  }
  
  return(sd$sd)
  
}
