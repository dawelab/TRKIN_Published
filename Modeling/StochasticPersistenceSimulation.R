## Test for persistence of Ab10(+) upon invasion by Ab10(-)

for(s in seq(0,0.4,0.01)) {

  # Set selection coefficients and drive strengths
  
  ha <- 0.25  ## Dominance coefficient for Ab10/N10
  hk <- 0.2  ## Dominance coefficient for K10/N10
  a <- 0.6   ## fitness cost of Ab10/Ab10
  k <- 0.225   ## fitness cost of K10/K10
  d1 <- 0.4  ## drive for Ab10 against N10 -> 60-80%
  d2 <- 0.1  ## weak drive for K10 against N10 -> 51-55%
  d3 <- 0.1   ## drive for Ab10 against K10 -> 52-60%
  Ne <- 10000   ## effective population size -> 5000-10000
  
  gen <- 10001  ##Number of generations to run each simulation
  
  
  #Number of 'iterations' to repeat the simulation (for each parameter set)
  iterations <- 10000 
  
  # Initialize an empty dataframe to store genotype counts for each paramater set
  gt_counts_all <- data.frame()
  
  # Repeat the simulation 'iterations' times (for each parameter set)
  for (i in 1:iterations) {
    
    
    # Initialize allele frequencies
    pmplus = 0.06                   # Initial frequency of Ab10_trkin+ in males
    pfplus = 0.06                   # Initial frequency of Ab10_trkin+ in females
    pfminus = 1/Ne                 # Initial frequency of Ab10_trkin- in females
    pmminus = 1/Ne                 # Initial frequency of Ab10_trkin- in males
    qf = 0.06                             # Initial frequency of K10L2 in females
    qm = 0.06                             # Initial frequency of K10L2 in males
    nf = 1-(pfplus+pfminus+qf)            # Initial frequency of N10 in females
    nm = 1-(pmplus+pmminus+qm)            # Initial frequency of N10 in males
    
    # Initialize generation count
    t = 0
    
    # Create a data frame to track counts of alleles in the population (initial state)
    gt_counts = data.frame(t = numeric(),
                           Ab10mplus=numeric(),
                           Ab10fplus=numeric(),
                           Ab10mminus=numeric(),
                           Ab10fminus=numeric(),
                           K10m=numeric(),
                           K10f=numeric(),
                           N10m=numeric(),
                           N10f=numeric(),
                           i=numeric())
    
    # Simulate for up to 'gen' generations
    while (t < gen & pfplus+pmplus>0 & pfminus+pmminus>0) {
      
      # Calculate Population Fitness
      wbar = (nm*nf) +
        (1 - ha*a)*(nm*pfplus + nf*pmplus) +
        (1 - ha*a)*(nm*pfminus + nf*pmminus) +
        (1 - hk*k)*(nm*qf + nf*qm) +
        (1 - a)*(pmplus*pfplus) +
        (1 - a)*(pmplus*pfminus + pfplus*pmminus) +
        (1 - ha*a)*(1 - hk*k)*(pmplus*qf + pfplus*qm) +
        (1 - a)*(pmminus*pfminus) +
        (1 - ha*a)*(1 - hk*k)*(pmminus*qf + pfminus*qm) +
        (1 - k)*(qm*qf)
      
      # Calculate next generation frequencies for all alleles
      pmplusnext = ((1 - ha*a)*(1/2)*(nm*pfplus + nf*pmplus) +
                      (1 - a)*(pmplus*pfplus) +
                      (1 - a)*(1/2)*(pmplus*pfminus + pfplus*pmminus) +
                      (1 - ha*a)*(1 - hk*k)*(1/2)*(pmplus*qf + pfplus*qm))/wbar
      
      pmminusnext = ((1 - ha*a)*(1/2)*(nm*pfminus + nf*pmminus) +
                       (1 - a)*(1/2)*(pmplus*pfminus + pfplus*pmminus) +
                       (1 - a)*(pmminus*pfminus) +
                       (1 - ha*a)*(1 - hk*k)*(1/2)*(pmminus*qf + pfminus*qm))/wbar
      
      qmnext = ((1 - hk*k)*(1/2)*(nm*qf + nf*qm) +
                  (1 - ha*a)*(1 - hk*k)*(1/2)*(pmplus*qf + pfplus*qm) +
                  (1 - ha*a)*(1 - hk*k)*(1/2)*(pmminus*qf + pfminus*qm) +
                  (1 - k)*(qm*qf))/wbar
      
      nmnext = ((nm*nf) +
                  (1 - ha*a)*(1/2)*(nm*pfplus + nf*pmplus) +
                  (1 - ha*a)*(1/2)*(nm*pfminus + nf*pmminus) +
                  (1 - hk*k)*(1/2)*(nm*qf + nf*qm))/wbar
      
      pfplusnext = ((1 - ha*a)*((1 + d1 - s)/2)*(nm*pfplus + nf*pmplus) +
                      (1 - a)*(pmplus*pfplus) +
                      (1 - a)*(1/2)*(pmplus*pfminus + pfplus*pmminus) +
                      (1 - ha*a)*(1 - hk*k)*((1 + d3)/2)*(pmplus*qf + pfplus*qm))/wbar
      
      pfminusnext = ((1 - ha*a)*((1 + d1)/2)*(nm*pfminus + nf*pmminus) +
                       (1 - a)*(1/2)*(pmplus*pfminus + pfplus*pmminus) +
                       (1 - a)*(pmminus*pfminus) +
                       (1 - ha*a)*(1 - hk*k)*((1 + d3)/2)*(pmminus*qf + pfminus*qm))/wbar
      
      qfnext = ((1 - hk*k)*((1 + d2)/2)*(nm*qf + nf*qm) +
                  (1 - ha*a)*(1 - hk*k)*((1 - d3)/2)*(pmplus*qf + pfplus*qm) +
                  (1 - ha*a)*(1 - hk*k)*((1 - d3)/2)*(pmminus*qf + pfminus*qm) +
                  (1 - k)*(qm*qf))/wbar
      
      nfnext = ((nm*nf) +
                  (1 - ha*a)*((1 - d1 + s)/2)*(nm*pfplus + nf*pmplus) +
                  (1 - ha*a)*((1 - d1)/2)*(nm*pfminus + nf*pmminus) +
                  (1 - hk*k)*((1 - d2)/2)*(nm*qf + nf*qm))/wbar
      
      # Use a multinomial distribution to simulate the next generation's allele frequencies
      gt = rmultinom(1, Ne, 
                     c(pmplusnext,
                       pfplusnext,
                       pmminusnext,
                       pfminusnext,
                       qmnext,
                       qfnext,
                       nmnext,
                       nfnext))
      
      # Calculate selective benefit of Ab10_trkin+, Ab10_trkin-, K10 during first generation
      if(t==0){
        selective_benefit_Ab10trkinplus=(pmplusnext+pfplusnext-pmplus-pfplus)/wbar
        selective_benefit_Ab10trkinminus=(pmminusnext+pfminusnext-pmminus-pfminus)/wbar
        selective_benefit_K10=(qmnext+qfnext-qm-qf)/wbar
        selective_benefit_N10=(nmnext+nfnext-nm-nf)/wbar
        
        relative_selective_benefit_Ab10trkinminus=((pfminusnext+pmminusnext)/(pfplusnext+pmplusnext))-1
        relative_selective_benefit_Ab10trkinplus=((pfplusnext+pmplusnext)/(pfminusnext+pmminusnext))-1
      }
      
      # Update the current frequencies for the next iteration
      pmplus = gt[1] / (Ne / 2)
      pfplus = gt[2] / (Ne / 2)
      pmminus = gt[3] / (Ne / 2)
      pfminus = gt[4] / (Ne / 2)
      qm = gt[5] / (Ne / 2)
      qf = gt[6] / (Ne / 2)
      nm = gt[7] / (Ne / 2)
      nf = gt[8] / (Ne / 2)
      
      
      # Add last generation's counts to the tracking data frame
      gt_counts[1,] = c(t, gt, i)
      
      
      # Increment generation count
      t = t + 1
      
    }
    
    
    # Append the genotype counts for each iteration to the dataframe
    gt_counts_all <- rbind(gt_counts_all, gt_counts)
    
    
  }
  
  gt_counts_all$ha = ha
  gt_counts_all$hk = hk
  gt_counts_all$a = a
  gt_counts_all$k = k
  gt_counts_all$d1 = d1
  gt_counts_all$s = s
  gt_counts_all$d2 = d2
  gt_counts_all$d3 = d3
  gt_counts_all$Ne = Ne
  gt_counts_all$selective_benefit_Ab10trkinplus=selective_benefit_Ab10trkinplus
  gt_counts_all$selective_benefit_Ab10trkinminus=selective_benefit_Ab10trkinminus
  gt_counts_all$selective_benefit_K10=selective_benefit_K10
  gt_counts_all$selective_benefit_N10=selective_benefit_N10
  gt_counts_all$relative_selective_benefit_Ab10trkinminus=relative_selective_benefit_Ab10trkinminus
  gt_counts_all$relative_selective_benefit_Ab10trkinplus=relative_selective_benefit_Ab10trkinplus
  
  
  #Write gt_counts file for all iterations of the simulation 
  write.csv(gt_counts_all,
            file=paste("Persistence_gt_counts_",
                       "ha",ha,"hk",hk,
                       "a",a,"k",k,
                       "d1",d1,"s",s,
                       "d2",d2,"d3",d3,
                       "Ne",Ne,
                       "gen",gen,
                       ".csv",
                       sep=""),
            quote=FALSE,
            row.names=FALSE)
  
}

# Access the loaded data frames as before
file_list <- list.files(pattern = "Persistence_gt_counts_ha[0-9.]+hk[0-9.]+a[0-9.]+k[0-9.]+d1[0-9.]+s[0-9.]+d2[0-9.]+d3[0-9.]+Ne[0-9.]+gen[0-9.]+\\.csv$")


combined_data <- data.frame()

# Loop through each file
for (file in file_list) {
  # Read the CSV file
  data <- read.csv(file)
  
  # Append the data to the combined_data data frame
  combined_data <- rbind(combined_data, data)
}

write.csv(combined_data, "Persistence_df_1.csv", row.names = FALSE)

