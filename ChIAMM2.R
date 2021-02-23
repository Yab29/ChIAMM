##g++ version need 4.9 or upper when install rstan package
##recommend use new version of R
##if you meet some problem about libxml-2.0.pc when install xml2 package, you need to find out the folder which include libxml-2.0.pc and add the path of this folder in variable "PKG_CONFIG_PATH"
library(getopt)

cat("==============================ChIAMM2 is running ==============================\n")
spec = matrix(c(
  'input','i', 1, "character","input file",
  'prefix','p', 2, "character","output prefix (default \"out\")",
  'inter','e', 0, "logical","input data is inter-chromosomal interactions,if not, intra-chromosomal interactions",
  'iter','r',2,"integer","iteration number (default 5000)",
  'warmup','w',2,"integer","warmup value (default 750)",
  'help','h', 0, "logical","print help"
), byrow=TRUE, ncol=5)
opt = getopt(spec) 

if (!is.null(opt$help) || is.null(opt$input)) {
   cat(paste(getopt(spec, usage = T), "\n"))
   q()
}
input_file=opt$input
output_prefix="out"
if (!is.null(opt$prefix)) {
   output_prefix=opt$prefix
}
type=4
if (!is.null(opt$inter)) {
   type=3
}
iter_number=5000
if (!is.null(opt$iter) && opt$iter >=0) {
   iter_number=opt$iter
}
warmup_value=750
if (!is.null(opt$warmup)) {
   warmup_value=opt$warmup
}
cat(paste("input file = ",input_file,"\n", sep = ""))
cat(paste("output prefix = ",output_prefix,"\n", sep = ""))
if (type==3) {
   cat("data type = inter chromosome\n")
} else {
   cat("data type = intra chromosome\n")
}
cat(paste("iteration number = ",iter_number,"\n", sep = ""))
cat(paste("warmup value = ",warmup_value,"\n", sep = ""))
cat("==============================================================================\n")

library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(bayesplot)
color_scheme_set("red")

##--------------------------------------------------------
ChIAMM2 <- read.delim(input_file) 

stan_code =paste("
  data {

      int<lower=1> N;               
      int<lower=2> y[N];         
      matrix[N,",type,"] x;               
      vector[N] mci;            
      real<lower=0> Amci;  
      
  }

  transformed data{
      int y_ast[N];     
       for (n in 1:N)
       y_ast[n] = y[n] - 2;
}


  parameters {

      real<lower=0> ci;   
      vector[",type,"] beta;     
      real alpha;         
      vector<lower=0, upper=1>[N] pi;        
      real<lower=0> phi;
      
}

  transformed parameters{
     vector[N] mu = exp(alpha + x*beta);
  }

  model {


  // Priors
    beta ~ normal(0,3);    
    alpha ~ normal(0,3);
    ci ~ normal(0,3); 
    pi ~ beta(mci,Amci);
    phi ~ inv_gamma(0.001,0.001);
    

  // Likelihood

  for(i in 1:N)
     target += log_mix(pi[i],
                    neg_binomial_2_lpmf(y_ast[i]|mu[i],phi),
                    neg_binomial_2_lpmf(y_ast[i]|mu[i]*exp(ci),phi));     

  }


  generated quantities{
    int y_hat[N]; 
    
    for(i in 1:N) {
      int which = bernoulli_rng(pi[i]);
      
      real mu_hat = pi[i] >= 0.500 ? mu[i]*exp(ci) : mu[i];
     
      if (log(mu_hat) > 17){
         y_hat[i] = 9999999;
      }else {
       y_hat[i] = 2+neg_binomial_2_rng(mu_hat,phi);
      }
     }
}
",sep="")

N=nrow(ChIAMM2)
y=ChIAMM2$ipet
mci <- ChIAMM2$tagcou1+ChIAMM2$tagcou2 - ChIAMM2$ipet    
mci[mci <= 0] <- 0.5  #this is the sape parameter should be above zero
Amci= mean(mci)

ChIAMM2[ChIAMM2$mappaAvg == 0, "mappaAvg"] <- 0.1  
ChIAMM2[ChIAMM2$selfAvg == 0, "selfAvg"] <- 0.5  

if (type==4) {
   x=matrix(c(log(ChIAMM2$selfAvg),log(ChIAMM2$gcAv),log(ChIAMM2$mappaAvg),log(ChIAMM2$distance)),
               nrow = nrow(ChIAMM2), byrow = FALSE) 
} else {
   x=matrix(c(log(ChIAMM2$selfAvg),log(ChIAMM2$gcAv),log(ChIAMM2$mappaAvg)),
               nrow = nrow(ChIAMM2), byrow = FALSE) 
}


input_data <- with(ChIAMM2,list(N=N,y=y,x=x,mci=mci,Amci=Amci)) 


resStan <- stan(model_code = stan_code, data = input_data,
                chains = 4,warmup = warmup_value,iter = iter_number, save_warmup = FALSE) 
                # if it is recommended, use control=list(adapt_delta=0.99) 
                # if it is recommended to increase max_treedepth above 10, 
                # use control=list(adapt_delta=0.99, max_treedepth=15)
                # change the iteration number and the warmup based on the trace plot and R-hat value

#to get the trace plot of beta, ci, pi and mu
traceplot(resStan, pars = c("alpha","beta","ci"), inc_warmup = FALSE)  
traceplot(resStan, pars = c("mu[10]","mu[20]","mu[30]","pi[10]","pi[20]","pi[30]"), inc_warmup = FALSE)
                # the numbers may be vary based on your random choice

y_hat <- as.matrix(resStan, pars = "y_hat") 
ppc_dens_overlay(log2(y), log2(y_hat[1:100, ]))  

a <- data.frame(extract(resStan, pars=c("pi"), permuted = TRUE, inc_warmup = FALSE))
pi = c()
for (i in 1:N) {
  pi <- c(pi, mean(a[,i]))
}

pi.round <- round(pi,1)

ChIAMM2_all_interaction <- data.frame(ChIAMM2$chrom1,ChIAMM2$start1,ChIAMM2$end1,
                               ChIAMM2$chrom2,ChIAMM2$start2,ChIAMM2$end2, ChIAMM2$ipet,pi.round)
colnames(ChIAMM2_all_interaction) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "ipet", "W1i")

ChIAMM2_significant = ChIAMM2_all_interaction %>% filter(ChIAMM2_all_interaction$W1i >= 0.5)
write.table(ChIAMM2_significant, paste(output_prefix,"significant_interaction.txt",sep="_"), sep="\t", row.names = F, quote=FALSE)
