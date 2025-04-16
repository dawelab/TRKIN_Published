###Maize Ab10trkin+/- model

##Set working directory
setwd("/Users/anjaligupta/Library/CloudStorage/OneDrive-UniversityofKansas/PhD/MaizeMeioticDrive/Persistence")


##1. Deterministic test of frequencies for the parameter subset

#Number of generations to run each simulation
gen <- 5000

# Set selection coefficients and drive strengths

ha <- 0.25  ## Dominance coefficient for Ab10/N10
hk <- 0.2  ## Dominance coefficient for K10/N10
a <- 0.6   ## fitness cost of Ab10/Ab10
k <- 0.225   ## fitness cost of K10/K10
d1 <- 0.4  ## drive for Ab10 against N10 -> 60-80%
s <- 0   ##s here is called delta1 in the model ## amount of Ab10 drive that decreases in presence of trkin 
d2 <- 0.1  ## weak drive for K10 against N10 -> 51-55%
d3 <- 0.1   ## drive for Ab10 against K10 -> 52-60%
Ne <- 10000   ## effective population size -> 5000-10000


# Initialize allele frequencies
pfplus = 1/Ne                        # Initial frequency of Ab10_trkin+ in females
pmplus = 1/Ne                          # Initial frequency of Ab10_trkin+ in males
pfminus = 1/Ne                        # Initial frequency of Ab10_trkin- in females
pmminus = 1/Ne                            # Initial frequency of Ab10_trkin- in males
qf = 1/Ne                             # Initial frequency of K10L2 in females
qm = 1/Ne                             # Initial frequency of K10L2 in males
nf = 1-(pfplus+pfminus+qf)              # Initial frequency of N10 in females
nm = 1-(pmplus+pmminus+qm)              # Initial frequency of N10 in males

# Initialize generation count
t = 0

# Create a data frame to track counts of alleles in the population (initial state)
gt_counts = data.frame(t = 0,
                       Ab10mplus=Ne/2 * pmplus,
                       Ab10fplus=Ne/2 * pfplus,
                       Ab10mminus=Ne/2 * pmminus,
                       Ab10fminus=Ne/2 * pfminus,
                       K10m=Ne/2 * qm,
                       K10f=Ne/2 * qf,
                       N10m=Ne/2 * nm,
                       N10f=Ne/2 * nf)

# Simulate for up to 'gen' generations
while (t < gen) {
  
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
  
  # Store next generation's allele frequencies
  gt = c(pmplusnext*(Ne/2),
         pfplusnext*(Ne/2),
         pmminusnext*(Ne/2),
         pfminusnext*(Ne/2),
         qmnext*(Ne/2),
         qfnext*(Ne/2),
         nmnext*(Ne/2),
         nfnext*(Ne/2))
  
  # Update the current frequencies for the next iteration
  pmplus = pmplusnext
  pfplus = pfplusnext
  pmminus = pmminusnext
  pfminus = pfminusnext
  qm = qmnext
  qf = qfnext
  nm = nmnext
  nf = nfnext
  
  # Increment generation count
  t = t + 1
  
  # Add the current generation's genotype counts to the tracking data frame
  gt_counts[length(gt_counts$t) + 1, ] = c(t, gt)
}


tail(gt_counts)










###2. Look at invasion of trkin- for a range of delta1 deterministically

#initialise an empty dataframe to store last generations genotype count for all s values
TrkinMinusInv <- data.frame()

for(s in seq(0,0.4,0.001)){ ##s here is called delta1 in the model
  
  #Number of generations to run each simulation
  gen <- 5000
  
  # Set selection coefficients and drive strengths
  
  ha <- 0.25  ## Dominance coefficient for Ab10/N10
  hk <- 0.2  ## Dominance coefficient for K10/N10
  a <- 0.6   ## fitness cost of Ab10/Ab10
  k <- 0.225   ## fitness cost of K10/K10
  d1 <- 0.4  ## drive for Ab10 against N10 -> 60-80%
  d2 <- 0.1  ## weak drive for K10 against N10 -> 51-55%
  d3 <- 0.1   ## drive for Ab10 against K10 -> 52-60%
  Ne <- 10000   ## effective population size -> 5000-10000


# Initialize allele frequencies
pfplus = 0.05                         # Initial frequency of Ab10_trkin+ in females
pmplus = 0.05                          # Initial frequency of Ab10_trkin+ in males
pfminus = 1/Ne                       # Initial frequency of Ab10_trkin- in females
pmminus = 1/Ne                            # Initial frequency of Ab10_trkin- in males
qf = 0.05                            # Initial frequency of K10L2 in females
qm = 0.05                             # Initial frequency of K10L2 in males
nf = 1-(pfplus+pfminus+qf)              # Initial frequency of N10 in females
nm = 1-(pmplus+pmminus+qm)              # Initial frequency of N10 in males

# Initialize generation count
t = 0

# Create a data frame to track counts of alleles in the population (initial state)
gt_counts = data.frame(t = 0,
                       Ab10mplus=Ne/2 * pmplus,
                       Ab10fplus=Ne/2 * pfplus,
                       Ab10mminus=Ne/2 * pmminus,
                       Ab10fminus=Ne/2 * pfminus,
                       K10m=Ne/2 * qm,
                       K10f=Ne/2 * qf,
                       N10m=Ne/2 * nm,
                       N10f=Ne/2 * nf,
                       s = s)

# Simulate for up to 'gen' generations
while (t < gen) {
  
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
  
  # Store next generation's allele frequencies
  gt = c(pmplusnext*(Ne/2),
         pfplusnext*(Ne/2),
         pmminusnext*(Ne/2),
         pfminusnext*(Ne/2),
         qmnext*(Ne/2),
         qfnext*(Ne/2),
         nmnext*(Ne/2),
         nfnext*(Ne/2))
  
  # Update the current frequencies for the next iteration
  pmplus = pmplusnext
  pfplus = pfplusnext
  pmminus = pmminusnext
  pfminus = pfminusnext
  qm = qmnext
  qf = qfnext
  nm = nmnext
  nf = nfnext
  
  # Increment generation count
  t = t + 1
  
  # Add the current generation's counts to the tracking data frame
  gt_counts[length(gt_counts$t) + 1, ] = c(t, gt, s)
}

#store last generations genotype count for all s values
TrkinMinusInv <- rbind(TrkinMinusInv, gt_counts[t+1,])


}




###3. Look at invasion of trkin+ for a range of delta1 deterministically

#initialise an empty dataframe to store last generations genotype count for all s values
TrkinPlusInv <- data.frame()

for(s in seq(0,0.4,0.01)){ ##s here is called delta1 in the model
  
  #Number of generations to run each simulation
  gen <- 5000
  
  # Set selection coefficients and drive strengths
  
  ha <- 0.25  ## Dominance coefficient for Ab10/N10
  hk <- 0.2  ## Dominance coefficient for K10/N10
  a <- 0.6   ## fitness cost of Ab10/Ab10
  k <- 0.225   ## fitness cost of K10/K10
  d1 <- 0.4  ## drive for Ab10 against N10 -> 60-80%
  d2 <- 0.1  ## weak drive for K10 against N10 -> 51-55%
  d3 <- 0.1   ## drive for Ab10 against K10 -> 52-60%
  Ne <- 10000   ## effective population size -> 5000-10000
  
  
  # Initialize allele frequencies
  pfplus = 1/Ne                        # Initial frequency of Ab10_trkin+ in females
  pmplus = 1/Ne                          # Initial frequency of Ab10_trkin+ in males
  pfminus = 0.05                        # Initial frequency of Ab10_trkin- in females
  pmminus = 0.05                            # Initial frequency of Ab10_trkin- in males
  qf = 0.05                             # Initial frequency of K10L2 in females
  qm = 0.05                             # Initial frequency of K10L2 in males
  nf = 1-(pfplus+pfminus+qf)              # Initial frequency of N10 in females
  nm = 1-(pmplus+pmminus+qm)              # Initial frequency of N10 in males
  
  # Initialize generation count
  t = 0
  
  # Create a data frame to track counts of alleles in the population (initial state)
  gt_counts = data.frame(t = 0,
                         Ab10mplus=Ne/2 * pmplus,
                         Ab10fplus=Ne/2 * pfplus,
                         Ab10mminus=Ne/2 * pmminus,
                         Ab10fminus=Ne/2 * pfminus,
                         K10m=Ne/2 * qm,
                         K10f=Ne/2 * qf,
                         N10m=Ne/2 * nm,
                         N10f=Ne/2 * nf,
                         s = s)
  
  # Simulate for up to 'gen' generations
  while (t < gen) {
    
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
    
    # Store next generation's allele frequencies
    gt = c(pmplusnext*(Ne/2),
           pfplusnext*(Ne/2),
           pmminusnext*(Ne/2),
           pfminusnext*(Ne/2),
           qmnext*(Ne/2),
           qfnext*(Ne/2),
           nmnext*(Ne/2),
           nfnext*(Ne/2))
    
    # Update the current frequencies for the next iteration
    pmplus = pmplusnext
    pfplus = pfplusnext
    pmminus = pmminusnext
    pfminus = pfminusnext
    qm = qmnext
    qf = qfnext
    nm = nmnext
    nf = nfnext
    
    # Increment generation count
    t = t + 1
    
    # Add the current generation's counts to the tracking data frame
    gt_counts[length(gt_counts$t) + 1, ] = c(t, gt, s)
  }
  
  #store last generations genotype count for all s values
  TrkinPlusInv <- rbind(TrkinPlusInv, gt_counts[t+1,])
  
  
}





####4.	Test the strength of selection for a range of delta1 

#Initialize empty dataframe to store 'se' - selective benefit values
se_df = data.frame()

for(s in seq(0,0.4,0.005)){ ##this s here is delta1 in the model ###amount of Ab10 drive that decreases in presence of trkin
  
  for(Ne in seq(100,10000,100)){ ## loop over different values of Ne
    
    
    #Number of generations to run each simulation
    gen <- 5000
    
    # Set selection coefficients and drive strengths
    
    ha <- 0.25  ## Dominance coefficient for Ab10/N10
    hk <- 0.2  ## Dominance coefficient for K10/N10
    a <- 0.6   ## fitness cost of Ab10/Ab10
    k <- 0.225   ## fitness cost of K10/K10
    d1 <- 0.4  ## drive for Ab10 against N10 -> 60-80%
    d2 <- 0.1  ## weak drive for K10 against N10 -> 51-55%
    d3 <- 0.1   ## drive for Ab10 against K10 -> 52-60%
    
    
    # Initialize allele frequencies
    pfplus = 1/Ne                        # Initial frequency of Ab10_trkin+ in females
    pmplus = 1/Ne                          # Initial frequency of Ab10_trkin+ in males
    pfminus = 0                        # Initial frequency of Ab10_trkin- in females
    pmminus = 0                            # Initial frequency of Ab10_trkin- in males
    qf = 1/Ne                             # Initial frequency of K10L2 in females
    qm = 1/Ne                             # Initial frequency of K10L2 in males
    nf = 1-(pfplus+pfminus+qf)              # Initial frequency of N10 in females
    nm = 1-(pmplus+pmminus+qm)              # Initial frequency of N10 in males
    
    # Initialize generation count
    t = 0
    
    # Simulate for up to 'gen' generations
    while (t < gen) {
      
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
      
      # Store next generation's allele frequencies
      gt = c(pmplusnext*(Ne/2),
             pfplusnext*(Ne/2),
             pmminusnext*(Ne/2),
             pfminusnext*(Ne/2),
             qmnext*(Ne/2),
             qfnext*(Ne/2),
             nmnext*(Ne/2),
             nfnext*(Ne/2))
      
      
      # Update the current frequencies for the next iteration
      pmplus = pmplusnext
      pfplus = pfplusnext
      pmminus = pmminusnext
      pfminus = pfminusnext
      qm = qmnext
      qf = qfnext
      nm = nmnext
      nf = nfnext
      
      # Increment generation count
      t = t + 1
      
    }
    
    
    #Introduce trkin- in females at generation 5000 at 1/Ne
    pfminus=1/Ne
    
    
    ## Run for one generation
    
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
    
    
    #Calculate 'se' - selective benefit of trkin-
    se = (((pfminusnext+pmminusnext)/(pfminus+pmminus))/((pfplusnext+pmplusnext)/(pfplus+pmplus)))-1
    
    #Store s, se, Ne and 2Nes values
    se_temp=data_frame(Ne=Ne,
                       delta1=s,
                       `RelativeSelectiveBenefitTrkin-`=se,
                       `2Nes`=2*Ne*se)
    
    ##Append se to 'se' dataframe
    se_df=rbind(se_df, se_temp)
    
  }
}


library(ggplot2)

ggplot(se_df,
       aes(x=Ne,
           y=delta1,
           z=`2Nes`,
           fill=`RelativeSelectiveBenefitTrkin-`)) +
  geom_tile() +
  scale_fill_gradient(low="pink",
                      high="skyblue") +
  theme_bw() +
  theme(legend.position = "bottom") +
  geom_contour(breaks = 1,
               col="black",
               lwd=1)






###5. Test how long Ab10trkin+ can persist in a population that is being invaded by Ab10trkin-

## Simulation is in a different R file

library(readr)
Data_1 <- read_csv("Persistence_df_1.csv")

library(ggplot2)
library(ggridges)
library(ggpubr)

FigureB <- ggplot(subset(Data_1,
                         !Ab10mminus==0 &
                           s %in% c(seq(0,0.19,0.01),0.29,0.39)), 
                  aes(x = t, 
                      y = as.factor(s), 
                      group = as.factor(s))) +
  geom_density_ridges(fill = "steelblue", alpha = 0.7) + 
  theme_bw() +
  xlim(0,1000) +
  labs(x = "Time to loss of Ab10 Trkin(+) [generations]", 
       y = "Reduction in drive in Ab10 Trkin(+) [0.1 = 5% drive]") +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20))

library(dplyr)


Data_1_1 <- Data_1 %>%
  group_by(s) %>%
  summarise(`Prop_trkin-_wins` = sum(Ab10mminus>0)/max(i))

FigureA <- ggplot(Data_1_1,
                  aes(y=s,
                      x=`Prop_trkin-_wins`)) +
  geom_point(col="steelblue",
             size=2) +
  geom_line(col="blue") +
  theme_bw() +
  labs(y="Reduction in drive in Ab10 Trkin(+) [0.1 = 5% drive]",
       x="Proportion Ab10 trkin(-) outcompetes Ab10 Trkin(+)") +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20))


Figure10 <- ggarrange(FigureA,
                     FigureB,
                     nrow = 1, ncol = 2, labels = "AUTO",
                     font.label = list(size = 20, face="bold"))
Figure10

pdf("Figure10.pdf",
    width = 14,
    height = 7.5)
print(Figure10)
dev.off()


