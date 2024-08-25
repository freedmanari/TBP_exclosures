# run after lyme_model.R

# code to create Figure 3 in the main text


bin_search <- function(xs, test) {
  b <- bin_search_index(xs,
                        function(x, history) return(ifelse(x %in% history$x,
                                                           history$stat[which(history$x == x)],
                                                           test(x))),
                        c())
  
  if(any(is.na(b))) {
    return(NA)
  }
  
  index <- b["index"]
  frac <- b["frac"]
  return(unname((1-frac)*xs[index-1] + frac*xs[index]))
}

bin_search_index <- function(xs, test, history) {
  if (length(xs) <= 1 || any(is.na(xs))) {
    return(c(index=NA, frac=NA))
  }
  
  if (length(xs) %% 2 == 0) {
    stat1 <- test(xs[length(xs)/2], history)
    if (stat1 > 0) {
      if (length(xs)==2) {
        return(c(index=NA, frac=NA))
      }
      
      stat2 <- test(xs[length(xs)/2-1], history)
      if (stat2 > 0) {
        if (length(xs)==4) {
          return(c(index=NA, frac=NA))
        }
        
        return(bin_search_index(xs[1:(length(xs)/2)], test,
                                rbind(history,
                                      data.frame(x=xs[c(length(xs)/2, length(xs)/2-1)],
                                                 stat=c(stat1, stat2)))))
      } else {
        return(c(index=length(xs)/2,
                 frac=stat2/(stat2-stat1)))
      }
    } else {
      stat2 <- test(xs[length(xs)/2+1], history)
      if (stat2 > 0) {
        return(c(index=length(xs)/2+1,
                 frac=stat1/(stat1-stat2)))
      } else {
        if (length(xs)==2) {
          return(c(index=NA, frac=NA))
        }
        
        return(bin_search_index(xs[(length(xs)/2+1):length(xs)], test,
                                rbind(history,
                                      data.frame(x=xs[c(length(xs)/2, length(xs)/2+1)],
                                                 stat=c(stat1, stat2)))) +
               c(length(xs)/2,0))
      }
    }
  } else {
    stat1 <- test(xs[(length(xs)+1)/2], history)
    if (stat1 > 0) {
      stat2 <- test(xs[(length(xs)-1)/2], history)
      if (stat2 > 0) {
        if (length(xs)==3) {
          return(c(index=NA, frac=NA))
        }
        
        return(bin_search_index(xs[1:((length(xs)-1)/2)], test,
                                rbind(history,
                                      data.frame(x=xs[c((length(xs)+1)/2, (length(xs)-1)/2)],
                                                 stat=c(stat1, stat2)))))
      } else {
        return(c(index=(length(xs)+1)/2,
                 frac=stat2/(stat2-stat1)))
      }
    } else {
      stat2 <- test(xs[(length(xs)+3)/2], history)
      if (stat2 > 0) {
        return(c(index=(length(xs)+3)/2,
                 frac=stat1/(stat1-stat2)))
      } else {
        if (length(xs)==3) {
          return(c(index=NA, frac=NA))
        }
        
        return(bin_search_index(xs[((length(xs)+3)/2):length(xs)], test,
                                rbind(history,
                                      data.frame(x=xs[c((length(xs)+1)/2, (length(xs)+3)/2)],
                                                 stat=c(stat1, stat2)))) +
               c((length(xs)+1)/2,0))
      }
    }
  }
}


test_outside_greater <- function(r_in) {
  out_pdes <- run_pdes(r_ins=r_in, v1s=v1, vTs=vT, beta_mults=beta_mult)
  
  avg_densities <- 
    out_pdes %>% 
    mutate(area = 2*pi*r*dr,
           freq = area*y) %>%
    group_by(outside=r>r_in) %>%
    summarise(avg_dens=sum(freq)/sum(area)) %>%
    ungroup()
  return(avg_densities$avg_dens[avg_densities$outside]-avg_densities$avg_dens[!avg_densities$outside])
}



r_ins <- rs[rs > 0 & rs < r_out * 3/4]

switch_r_ins <- data.frame()

for (beta_mult in c(beta_mult_low,beta_mult_mid,beta_mult_high)) {
for (vT in c(0,5)) {
  print(noquote(paste0("beta multiplier = ",beta_mult,", with",ifelse(vT==0,"out","")," tick movement")))
  print(noquote("v1 ="))
for (v1 in seq(v1_low, v1_high, length=21)) {
  print(v1)
  
  if (test_outside_greater(r_ins[1])>0 || test_outside_greater(r_ins[length(r_ins)])<0) {
    switch_r_ins <- rbind(switch_r_ins,
                          data.frame(beta_mult=beta_mult,
                                     tick_movement=(vT!=0),
                                     v1=v1,
                                     switch_r_in=NA))
  } else {
    switch_r_ins <- rbind(switch_r_ins,
                          data.frame(beta_mult=beta_mult,
                                     tick_movement=(vT!=0),
                                     v1=v1,
                                     switch_r_in=bin_search(r_ins, test_outside_greater)))
  }
}
}
}


switch_r_ins %>% 
  ggplot() +
  geom_line(aes(x=v1, y=switch_r_in,
                linetype=tick_movement, color=as.factor(beta_mult),
                group=interaction(tick_movement,as.factor(beta_mult))), linewidth=.7) +
  coord_cartesian(ylim=c(0,300)) +
  scale_x_continuous(expand=expansion(c(0,0))) +
  scale_y_continuous(expand=expansion(c(0,0)),breaks=seq(0,300,100)) +
  theme_classic() +
  xlab("average daily rodent movement (m/day)") +
  ylab("radius of exclosure (m)\nat which average NQi is equal\ninside & outside of exclosure") +
  scale_color_discrete(name="tick-host\ncontact rates",labels=c("low","medium","high")) +
  scale_linetype_discrete(name="questing tick\nmovement",labels=c("0 m/day","5 m/day"))





