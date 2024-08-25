require(ReacTran)
require(tidyverse)
require(lemon)

odes <- function(t, y, parms) {
  names(y) <- var_names
  with(as.list(c(y,parms)), {
    H1 <- H1_s+H1_i+H1_r
    
    AF1 <- AF1_s+AF1_i
    
    TF1 <-
      LF1_s+LF1_e+
      NF1_s+NF1_e+NF1_i+
      AF1_s+AF1_i
    TF2 <-
      LF2_s+
      NF2_s+NF2_i+
      AF2
    
    mu1 <- d1 + H1/K1 * (b1-d1)
    mu2 <- d2 + H2/K2 * (b2-d2)
    
    
    dH1_s <- b1*H1 - (beta1N*NQ_i+beta1A*AQ_i)*H1_s - mu1*H1_s
    dH1_i <- (beta1N*NQ_i+beta1A*AQ_i)*H1_s - (gamma+mu1+alpha)*H1_i
    dH1_r <- gamma*H1_i - mu1*H1_r
    dH2 <- b2*H2 - mu2*H2
    
    dLQ <- bT/(1+c*TF1)*sigmaA*AF1 + bT/(1+c*TF2)*sigmaA*AF2 -
      (beta1L*H1+beta2L*H2)*LQ - dL*LQ
    dLF1_s <- exp(-(lambdaN*NF1_i+lambdaA*AF1_i)/H1)*beta1L*(H1_s+H1_r)*LQ - sigmaL*LF1_s
    dLF1_e <- beta1L*H1_i*LQ + (1-exp(-(lambdaN*NF1_i+lambdaA*AF1_i)/H1))*beta1L*(H1_s+H1_r)*LQ -
      sigmaL*LF1_e
    dLF2_s <- beta2L*H2*LQ - sigmaL*LF2_s
    
    dNQ_s <- mL*sigmaL*(LF1_s+LF2_s) - (beta1N*H1+beta2N*H2)*NQ_s - dN*NQ_s
    dNQ_i <- mL*sigmaL*LF1_e - (beta1N*H1+beta2N*H2)*NQ_i - dN*NQ_i
    dNF1_s <- exp(-(lambdaN*NF1_i+lambdaA*AF1_i)/H1)*beta1N*(H1_s+H1_r)*NQ_s - sigmaN*NF1_s
    dNF1_e <- beta1N*H1_i*NQ_s + (1-exp(-(lambdaN*NF1_i+lambdaA*AF1_i)/H1))*beta1N*(H1_s+H1_r)*NQ_s -
      sigmaN*NF1_e
    dNF1_i <- beta1N*H1*NQ_i - sigmaN*NF1_i
    dNF2_s <- beta2N*H2*NQ_s - sigmaN*NF2_s
    dNF2_i <- beta2N*H2*NQ_i - sigmaN*NF2_i
    
    dAQ_s <- mN*sigmaN*(NF1_s+NF2_s) - (beta1A*H1+beta2A*H2)*AQ_s - dA*AQ_s
    dAQ_i <- mN*sigmaN*(NF1_e+NF1_i+NF2_i) - (beta1A*H1+beta2A*H2)*AQ_i - dA*AQ_i
    dAF1_s <- beta1A*H1*AQ_s - sigmaA*AF1_s
    dAF1_i <- beta1A*H1*AQ_i - sigmaA*AF1_i
    dAF2 <- beta2A*H2*(AQ_s+AQ_i) - sigmaA*AF2
    
    dy <- lapply(1:N_var, function(i) get(paste0("d", var_names[i])))
    return(list(dy))
  })
}


pdes <- function(t, y, parms) {
  with(as.list(parms), {
    rs_out <- rs[rs >= r_in]
    N_r_in <- sum(r_mids < r_in)
    # indices_in <- rep(0:(N_theta*N_var-1), each=N_r_in)*N_r + 1:N_r_in
    
    for (i in 1:N_var) {
      var_name <- var_names[i]
      var_name_out <- paste0(var_name, "_out")
      var_value <- matrix(data=y[((i-1)*N_grid+1):(i*N_grid)], nrow=N_r, ncol=N_theta)
      var_value_out <- var_value[r_mids > r_in,]
      
      assign(var_name, var_value)
      assign(var_name_out, var_value_out)
      
      if (!grepl("Q", var_name)) {
        if (grepl("2", var_name)) {
          assign(paste0("tran_", var_name),
                 rbind(matrix(nrow=N_r_in, ncol=N_theta, data=0),
                       tran.polar(var_value_out, D.r = rho2, D.theta = rho2, r = rs_out, theta = thetas,
                                  C.r.up = first(var_value_out), C.r.down = last(var_value_out),
                                  cyclicBnd = 2)$dC))
          # assign(paste0("d", var_name), rbind(matrix(nrow=N_r_in, ncol=N_theta, data=0), get(paste0("tran_", var_name))))
        } else {
          assign(paste0("tran_", var_name),
                 tran.polar(var_value, D.r = rho1, D.theta = rho1, r = rs, theta = thetas,
                            C.r.up = first(var_value), C.r.down = last(var_value),
                            cyclicBnd = 2)$dC)
          # assign(paste0("d", var_name), get(paste0("tran_", var_name)))
        }
      } else {
        assign(paste0("tran_", var_name),
               tran.polar(var_value, D.r = rhoT, D.theta = rhoT, r = rs, theta = thetas,
                          C.r.up = first(var_value), C.r.down = last(var_value),
                          cyclicBnd = 2)$dC)
      }
    }
    
    
    
    H1 <- H1_s+H1_i+H1_r
    
    AF1 <- AF1_s+AF1_i
    
    TF1 <-
      LF1_s+LF1_e+
      NF1_s+NF1_e+NF1_i+
      AF1_s+AF1_i
    TF2 <-
      LF2_s+
      NF2_s+NF2_i+
      AF2
    
    mu1 <- d1 + H1/K1 * (b1-d1)
    mu2 <- d2 + H2/K2 * (b2-d2)
    
    
    dH1_s <- b1*H1 - (beta1N*NQ_i+beta1A*AQ_i)*H1_s - mu1*H1_s + tran_H1_s
    dH1_i <- (beta1N*NQ_i+beta1A*AQ_i)*H1_s - (gamma+mu1+alpha)*H1_i + tran_H1_i
    dH1_r <- gamma*H1_i - mu1*H1_r + tran_H1_r
    dH2 <- b2*H2 - mu2*H2 + tran_H2
    
    dLQ <- bT/(1+c*TF1)*sigmaA*AF1 + bT/(1+c*TF2)*sigmaA*AF2 -
      (beta1L*H1+beta2L*H2)*LQ - dL*LQ# + tran_LQ
    dLF1_s <- exp(-(lambdaN*NF1_i+lambdaA*AF1_i)/H1)*beta1L*(H1_s+H1_r)*LQ - sigmaL*LF1_s + tran_LF1_s
    dLF1_e <- beta1L*H1_i*LQ + (1-exp(-(lambdaN*NF1_i+lambdaA*AF1_i)/H1))*beta1L*(H1_s+H1_r)*LQ -
      sigmaL*LF1_e + tran_LF1_e
    dLF2_s <- beta2L*H2*LQ - sigmaL*LF2_s + tran_LF2_s
    
    dNQ_s <- mL*sigmaL*(LF1_s+LF2_s) - (beta1N*H1+beta2N*H2)*NQ_s - dN*NQ_s + tran_NQ_s
    dNQ_i <- mL*sigmaL*LF1_e - (beta1N*H1+beta2N*H2)*NQ_i - dN*NQ_i + tran_NQ_i
    dNF1_s <- exp(-(lambdaN*NF1_i+lambdaA*AF1_i)/H1)*beta1N*(H1_s+H1_r)*NQ_s - sigmaN*NF1_s + tran_NF1_s
    dNF1_e <- beta1N*H1_i*NQ_s + (1-exp(-(lambdaN*NF1_i+lambdaA*AF1_i)/H1))*beta1N*(H1_s+H1_r)*NQ_s -
      sigmaN*NF1_e + tran_NF1_e
    dNF1_i <- beta1N*H1*NQ_i - sigmaN*NF1_i + tran_NF1_i
    dNF2_s <- beta2N*H2*NQ_s - sigmaN*NF2_s + tran_NF2_s
    dNF2_i <- beta2N*H2*NQ_i - sigmaN*NF2_i + tran_NF2_i
    
    dAQ_s <- mN*sigmaN*(NF1_s+NF2_s) - (beta1A*H1+beta2A*H2)*AQ_s - dA*AQ_s + tran_AQ_s
    dAQ_i <- mN*sigmaN*(NF1_e+NF1_i+NF2_i) - (beta1A*H1+beta2A*H2)*AQ_i - dA*AQ_i + tran_AQ_i
    dAF1_s <- beta1A*H1*AQ_s - sigmaA*AF1_s + tran_AF1_s
    dAF1_i <- beta1A*H1*AQ_i - sigmaA*AF1_i + tran_AF1_i
    dAF2 <- beta2A*H2*(AQ_s+AQ_i) - sigmaA*AF2 + tran_AF2
    
    
    
    dy <- do.call(cbind, lapply(1:N_var, function(i) get(paste0("d", var_names[i]))))
    return(list(dy))
    
  })
}



N_r <- 800 # number of slices of r
r_out <- 1000 #radius of total area

rs <- seq(0, r_out, len = N_r+1)
dr <- rs[2] - rs[1]
r_mids <- (rs[-1] + rs[-(N_r+1)]) / 2
N_theta <- 3 # minimum allowable number of slices of theta
thetas <- seq(0, 2*pi, len = N_theta+1)
N_grid <- N_r*N_theta
y0_dens <- # initial conditions for the non-spatial ODEs
  c(H1_s = .0015, H1_i = 0, H1_r = 0,
    LQ = .001, LF1_s = .001, LF1_e = .001,
    NQ_s = .001, NQ_i = .001, NF1_s = .001, NF1_e = .001, NF1_i = .001,
    AQ_s = .001, AQ_i = .001, AF1_s = .001, AF1_i = .001,
    H2 = .00001,
    LF2_s = .001,
    NF2_s = .001, NF2_i = .001,
    AF2 = .001)
var_names <- names(y0_dens)
N_var <- length(y0_dens)
vars_out <- grep("2", var_names) #indices of variables present only outside exclosure
cross_sections_out <- rep(vars_out, each=N_theta)*N_theta + rep(-N_theta:-1,times=length(vars_out))


v1_low <- 10
v1_mid <- 25
v1_high <- 50

beta_mult_low <- 2
beta_mult_mid <- 3
beta_mult_high <- 5

run_pdes <- function(r_ins, v1s=v1_mid, v2=250, vTs=0, beta_mults=beta_mult_mid, beta2N_reductions=1, track_vars=c("NQ_i"),
                     track_time=FALSE, max_time=ifelse(track_time,100,3000), print_progress=TRUE) {
  out <- data.frame()
  
  for (beta_mult in beta_mults) {
    if (print_progress && length(beta_mults) > 1) {
      print(noquote(paste("beta_mult =",beta_mult)))
    }
  for (beta2N_reduction in beta2N_reductions) {
    if (print_progress && length(beta2N_reductions) > 1) {
      print(noquote(paste("   beta2N_reduction =",beta2N_reduction)))
    }
    
    parms <-
      c(b1 = .00821, b2 = 0, d1 = .0037, d2 = 0, K1 = .0015, K2 = .00001,
        bT = 2000, c = 2500, dL = .0365, dN = .015, dA = .00625, sigmaL = .28, sigmaN = .22, sigmaA = .12, mL = .15, mN = .15,
        beta1L = beta_mult*150, beta1N = beta_mult*5, beta1A = beta_mult*.2, beta2L = beta_mult*800, beta2N = beta2N_reduction*beta_mult*2500, beta2A = beta_mult*2500,
        gamma = .3, alpha = .33,
        lambdaN = .55, lambdaA = .55)
      
    ts_odes <- seq(0,3000,length=2)
    out_odes <- ode(y=y0_dens, times=ts_odes, parms=parms, func=odes) # run non-spatial ODEs to equilibrium
    
  for (v1 in v1s) {
    if (print_progress && length(v1s) > 1) {
      print(noquote(paste(" v1 =",v1)))
    }
  for (vT in vTs) {
    if (print_progress && length(vTs) > 1) {
      print(noquote(paste("  vT =",vT)))
    }
  for (r_in in r_ins) {
    if (print_progress && length(r_ins) > 1) {
      print(noquote(paste("   r_in =",r_in)))
    }
    
    r_in <- rs[which.min(abs(r_in-rs))] #ensures r_in is in rs
    N_r_in <- sum(r_mids < r_in)
    
    parms <-
      c(b1 = .00821, b2 = 0, d1 = .0037, d2 = 0, K1 = .0015, K2 = .00001,
        bT = 2000, c = 2500, dL = .0365, dN = .015, dA = .00625, sigmaL = .28, sigmaN = .22, sigmaA = .12, mL = .15, mN = .15,
        beta1L = beta_mult*150, beta1N = beta_mult*5, beta1A = beta_mult*.2, beta2L = beta_mult*800, beta2N = beta2N_reduction*beta_mult*2500, beta2A = beta_mult*2500,
        gamma = .3, alpha = .33,
        lambdaN = .55, lambdaA = .55,
        rho1 = v1^2/pi,
        rho2 = v2^2/pi,
        rhoT = vT^2/pi,
        r_in = r_in)
    
    y0 <- unname(rep(out_odes[length(ts_odes),-1], each=N_grid))
    y0[1:N_r_in + rep(cross_sections_out*N_r, each=N_r_in)] <- 0 # initial conditions for spatial PDEs taken from equilibrium output of non-spatial ODEs
    
    ts <- seq(0, max_time, length=ifelse(track_time, 501, 2))
    out_pdes <- ode.2D(y = y0, times = ts, parms = parms,
                       func = pdes, dimens = c(N_r, N_theta),
                       lrw = 1e7, cyclicBnd = 2) #might need to increase lrw
    
    for (i in which(var_names %in% track_vars)) {
      out <-
        rbind(out,
              out_pdes[, c(1, ((i-1)*N_grid+2):((i-1)*N_grid+N_r+1))] %>%
                as.data.frame() %>%
                pivot_longer(cols=-time, values_to="y", names_to="r") %>%
                mutate(r_in = r_in,
                       v1 = v1,
                       v2 = v2,
                       vT = vT,
                       beta_mult=beta_mult,
                       beta2N_reduction=beta2N_reduction,
                       var = var_names[i],
                       r = r_mids[as.numeric(r) - (i-1)*N_grid],
                       y = y * 1e4)) # converting from densities per m^2 to densities per ha
    }
  }
  }
  }
  }
  }
  
  if (track_time) {
    return(out)
  }
  return(out %>% filter(time > 0))
}



#####
# code to make Figure 2 in the main text

eqs <- run_pdes(r_ins=rs[rs>0 & rs<=300],
                vTs=c(0,5),v1s=c(v1_low,v1_mid,v1_high))

facet_labels <- c("low rodent movement\n(10 m/day)","medium rodent movement\n(25 m/day)","high rodent movement\n(50 m/day)")
names(facet_labels) <- c(v1_low, v1_mid,v1_high)


eqs %>%
  mutate(area = 2*pi*r*dr,
         freq = area*y,
         tick_mvmt = (vT > 0)) %>%
  group_by(v1, r_in, outside=r>r_in, tick_mvmt) %>%
  summarise(avg_dens = sum(freq)/sum(area)) %>%
  ungroup() %>% 
  ggplot() +
  geom_line(aes(x=r_in, y=avg_dens, color=outside, lty=tick_mvmt), linewidth=.7) +
  facet_rep_grid(~ v1, labeller=labeller(v1=facet_labels), scales="free", repeat.tick.labels = T) +
  scale_color_discrete(name="", labels=c("inside exclosure","outside exclosure")) +
  scale_x_continuous(name="radius of exclosure (m)",limits=c(0,300),expand=c(0,0)) +
  scale_y_continuous(name="average density of infected\nquesting nymphs (per ha)",
                     breaks=seq(0,18,3),expand=expansion(c(0,.01))) +
  scale_linetype_discrete(name="questing tick\nmovement",labels=c("0 m/day","5 m/day")) +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line.x=element_line(),
        axis.line.y=element_line(),
        strip.background=element_blank(),
        strip.text=element_text(size=10),
        panel.spacing = unit(1, "lines"),
        legend.position = "bottom",
        legend.text = element_text(size=10))




####
# code to make Figure 4 in the main text

eqs_no_tick_movement <- run_pdes(r_ins = c(50,150))
outs_no_tick_movement <- run_pdes(r_ins = c(50,150), track_time=TRUE, max_time=200)

for (r_in_plot in c(50,150)) {
  print(
    eqs_no_tick_movement %>%
    filter(r <= 200, r_in==r_in_plot) %>%
    ggplot() +
    geom_tile(aes(x=0, y=r, fill=y, color=y)) +
    scale_fill_gradientn(colors=rev(rainbow(7))[-1],
                         name="density of infected\nquesting nymphs\n(per ha)") +
    scale_color_gradientn(colors=rev(rainbow(7))[-1],
                          name="density of infected\nquesting nymphs\n(per ha)") +
    coord_polar() +
    geom_hline(aes(yintercept=r_in), linetype="dashed", linewidth=1) +
    scale_x_continuous(expand=expansion(c(0,0))) +
    scale_y_continuous(expand=expansion(c(0,0)), name="distance from center of exclosure (m)") +
    guides(fill = guide_colorbar(barheight = unit(1.5,"in"),
                                 ticks.colour = "black",
                                 ticks.linewidth = .5,
                                 frame.colour = "black",
                                 frame.linewidth = .5,
                                 title.hjust = .5)) +
    theme(panel.background=element_blank(),
          panel.grid=element_blank(),
          axis.line.y=element_line(),
          axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.title=element_text(size=9),
          axis.line.x=element_blank(),
          legend.position="right"))
  
  print(
    outs_no_tick_movement %>%
    filter(r<=200, r_in==r_in_plot) %>%
    ggplot() +
    geom_tile(aes(x=time, y=r, fill=y, color=y)) +
    scale_fill_gradientn(colors=rev(rainbow(7))[-1],
                         name="density of infected\nquesting nymphs\n(per ha)",
                         breaks=(if (r_in_plot==50) seq(4,12,2) else seq(2,8,2))) +
    scale_color_gradientn(colors=rev(rainbow(7))[-1],
                          name="density of infected\nquesting nymphs\n(per ha)",
                          breaks=(if (r_in_plot==50) seq(4,12,2) else seq(2,8,2))) +
    scale_x_continuous(expand=expansion(c(0,0)), name="time since exclosure installed (days)") +
    scale_y_continuous(expand=expansion(c(0,0)), name="distance from center of exclosure (m)") +
    geom_hline(aes(yintercept=r_in), linetype="dashed", linewidth=1) +
    guides(fill = guide_colorbar(barheight = unit(1.5,"in"),
                                 ticks.colour = "black",
                                 ticks.linewidth = .5,
                                 frame.colour = "black",
                                 frame.linewidth = .5,
                                 title.hjust = .5)) +
    theme(axis.line=element_line(),
          legend.title=element_text(size=9),
          legend.position="right"))
}


#####
# code to make Figure 5 in the main text

eqs_tick_movement <- run_pdes(r_ins = c(50,150), v1s=v1_high, vTs=5)
outs_tick_movement <- run_pdes(r_ins = c(50,150), v1s=v1_high, vTs=5, track_time=TRUE, max_time=200)

for (r_in_plot in c(50,150)) {
  print(
    eqs_tick_movement %>%
      filter(r <= 200, r_in==r_in_plot) %>%
      ggplot() +
      geom_tile(aes(x=0, y=r, fill=y, color=y)) +
      scale_fill_gradientn(colors=rev(rainbow(7))[-1],
                           name="density of infected\nquesting nymphs\n(per ha)") +
      scale_color_gradientn(colors=rev(rainbow(7))[-1],
                            name="density of infected\nquesting nymphs\n(per ha)") +
      coord_polar() +
      geom_hline(aes(yintercept=r_in), linetype="dashed", linewidth=1) +
      scale_x_continuous(expand=expansion(c(0,0))) +
      scale_y_continuous(expand=expansion(c(0,0)), name="distance from center of exclosure (m)") +
      guides(fill = guide_colorbar(barheight = unit(1.5,"in"),
                                   ticks.colour = "black",
                                   ticks.linewidth = .5,
                                   frame.colour = "black",
                                   frame.linewidth = .5,
                                   title.hjust = .5)) +
      theme(panel.background=element_blank(),
            panel.grid=element_blank(),
            axis.line.y=element_line(),
            axis.text.x=element_blank(),
            axis.title.x=element_blank(),
            axis.ticks.x=element_blank(),
            legend.title=element_text(size=9),
            axis.line.x=element_blank(),
            legend.position="right"))
  
  print(
    outs_tick_movement %>%
      filter(r<=200, r_in==r_in_plot) %>%
      ggplot() +
      geom_tile(aes(x=time, y=r, fill=y, color=y)) +
      scale_fill_gradientn(colors=rev(rainbow(7))[-1],
                           name="density of infected\nquesting nymphs\n(per ha)") +
      scale_color_gradientn(colors=rev(rainbow(7))[-1],
                            name="density of infected\nquesting nymphs\n(per ha)") +
      scale_x_continuous(expand=expansion(c(0,0)), name="time since exclosure installed (days)") +
      scale_y_continuous(expand=expansion(c(0,0)), name="distance from center of exclosure (m)") +
      geom_hline(aes(yintercept=r_in), linetype="dashed", linewidth=1) +
      guides(fill = guide_colorbar(barheight = unit(1.5,"in"),
                                   ticks.colour = "black",
                                   ticks.linewidth = .5,
                                   frame.colour = "black",
                                   frame.linewidth = .5,
                                   title.hjust = .5)) +
      theme(axis.line=element_line(),
            legend.title=element_text(size=9),
            legend.position="right"))
}








######
# code to make Figure S1 in the supplement

eqs_beta2N_reduction <- run_pdes(r_ins = 150, v1s=v1_high, vTs=5, beta2N_reductions=c(.1,.25,1))


for (beta2N_reduction_plot in c(.1,.25,1)) {
  print(
    eqs_beta2N_reduction %>%
      filter(r <= 200, beta2N_reduction==beta2N_reduction_plot) %>%
      ggplot() +
      geom_tile(aes(x=0, y=r, fill=y, color=y)) +
      scale_fill_gradientn(colors=rev(rainbow(7))[-1],
                           name="density of infected\nquesting nymphs\n(per ha)") +
      scale_color_gradientn(colors=rev(rainbow(7))[-1],
                            name="density of infected\nquesting nymphs\n(per ha)") +
      coord_polar() +
      geom_hline(aes(yintercept=r_in), linetype="dashed", linewidth=1) +
      scale_x_continuous(expand=expansion(c(0,0))) +
      scale_y_continuous(expand=expansion(c(0,0)), name="distance from center of exclosure (m)") +
      guides(fill = guide_colorbar(barheight = unit(1.5,"in"),
                                   ticks.colour = "black",
                                   ticks.linewidth = .5,
                                   frame.colour = "black",
                                   frame.linewidth = .5,
                                   title.hjust = .5)) +
      theme(panel.background=element_blank(),
            panel.grid=element_blank(),
            axis.line.y=element_line(),
            axis.text.x=element_blank(),
            axis.title.x=element_blank(),
            axis.ticks.x=element_blank(),
            legend.title=element_text(size=9),
            axis.line.x=element_blank(),
            legend.position="right")
  )
}
