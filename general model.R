#' ---
#' title: “Angry animals - linking changes in aggression after habitat loss to demographic changes at the population level using IBM "
#' author: “LBE"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output: github_document
#' ---
#'
#'primer of ecology (in R)- -elli

#' - **1. Initiate habitat - binary** 
#' * outside loop so that habitat is the same for entire model...*
#' For loop 1 - repeat for each replicate
#' -  2. Move
#' -  3. Fight - dependent o energy
#' -  4. Feed - 
#' -  5. Reproduce
#' -  6. Energy gain/loss
#' End FL1 
#' 
#' =======================================================================
#' # 1. Workspace preparation 
#' =======================================================================

rm(list=ls()) 

#' Load packages
packages = c('reshape','gstat','ggplot2',
             'viridis','tictoc', 'progress')
load.pack = lapply(packages,require,char=T)
load.pack

#' =======================================================================
#' # 2. Parametres
#' =======================================================================

#' ## IBM Parameters *should include constraints on each parametre* #double check how many are needed in reduced model
ngenerations  = 10    # No. generations
replicates    = 10    # No. replicates (first half trans. second half intrans.)
dim           = 100   # dimension of square habitat array
hab_dim       = dim^2 # total no. cells
nspecies      = 20    # No. species
nindiv        = 20    # No. individuals per species
tot_indiv     = nspecies*nindiv # Total individuals
eloss         = 30    # Time-step energy loss
fight_eloss   = 15    # Energy loss from aggression
exp_fight     = 20    # Roulette selection exponent
a             = 1
b             = 100
v             = 0.08  # Logistic mortality parameters
a_rep         = 0.4   # Asymptotic reproduction probability
offspring_pen = 2     # factor by which to scale offspring energy
repro_pen     = 20    # energy loss from reproduction
diff_ag       = 0.05  # Differential aggression \in [0,0.5)

#' Potential one time-step movement. Can move to any adjacent cell (inc diagonal) 
step_moves = c(0,-1,dim-1,dim,dim+1,1,-dim+1,-dim,-dim-1)      #added option to stay in the same place

#' Potential aggression matrix - M1 assume all species have equal competitive potential - so outcome will only depend on energy levels
aggression = matrix(0.5, ncol=nspecies, nrow = nspecies)
diag(aggression) = 0.5 # For interaction between individuals of the same species 

#' Store species richness at each generation
rich = data.frame(gen=rep(1:ngenerations,replicates),
                  rep=rep(1:replicates,each=ngenerations),
                  rich=NA,
                  type=NA)

#' Progress bar  *use pb$tick() in for-loop to run*
pb = progress_bar$new(total=100, width=100, clear=F)

#' =======================================================================
#' #3a. Function - habitat
#' ======================================================================= 
# WORKING IN ISOLATION
#' ## 0. Initiate habitat and individuals --------------------------------
init_habitat <- function(){
  
  # Simulate habitat 
  hab_grid = expand.grid(1:dim, 1:dim)                                     # The habitat arena
  names(hab_grid) = c('x','y')
  hab_bin = sample(x=c(1,0), size = hab_dim, replace = T)                  # samples binary habitat values
  hab_bin = cbind(hab_grid, hab_bin)
  
  #' Make matrix
  hab_vals = as.matrix(cast(hab_bin,
                            x~y,
                            value='hab_bin'))
  
  #' ## Plot habitat values
  hab.plot = ggplot(hab_bin) + theme_bw() +
    geom_raster(aes(x,y,fill=hab_bin)) +
    scale_x_continuous(expand=expand_scale(add=0)) +
    scale_y_continuous(expand=expand_scale(add=0)) +
    theme(legend.position = 'top',
          axis.text=element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  hab.plot
  
  #' 
  #' Randomly place individuals of each species in habitat
  #' 
  
  locs = cbind(sample(hab_dim,tot_indiv),rep(1:nspecies,each=nindiv))       # Randomly picks locations (V1-cell) for all individuals (V2) 
  locs = cbind(locs,hab_vals[locs[,1]])                                     # habitat value
  
  #' add energetic value
  locs = cbind(locs,rep(100,nindiv))                                        # starting energetic value =100?
  colnames(locs)=c("loc", "sp", "hab_val", "e_val")
  
  #' Place individuals in habitat
  habitat = matrix(0,ncol=dim,nrow=dim)                                     # make 100 x 100 matrix
  habitat[locs[,1]] = 1                                                     # pulls out uniqe cell number from locs, and places individual (=1) in the cell
  
  #' Sort by location
  locs = locs[order(locs[,1]),] 
  
  return(locs, habitat)
}

#' =======================================================================
#' #3b. Function - move
#' =======================================================================
#WORKS IN ISOLATION + HABITAT
move <- function(locs){      
  edge_corner <- unique(c(1:dim,                                       # top row
                          which(1:hab_dim%%dim%in%0:1),                 # sides
                          (hab_dim-dim+1):hab_dim))                     # bottom
  ec_occ      <- edge_corner[which(habitat[edge_corner]==1)]           # which edge corner cells are occupied
  ec_indiv    <- locs[locs[,1]%in%ec_occ,]                             # which individuals are on edge corner cells (V2&V1)
  offec_indiv <- locs[locs[,1]%in%setdiff(locs[,1],ec_indiv[,1]),]     # individuals  NOT on edge corner cells (V1&v2)
  
  offec_indiv <- cbind(offec_indiv, sample(step_moves, dim(offec_indiv)[1], replace = TRUE)) # add move step (V3)
  
  #-------------------------------------------- edge and corner individuals ---------------------------------------
  ec_indiv_new <- NULL
  if(sum(ec_indiv[,1]==1)>0){                                         # Q1 find if any inds are in top left corner
    top_left_ind = ec_indiv[ec_indiv[,1]==1,]                         # VALUE IF TRUE -pull out that data
    if(is.null(dim(top_left_ind))){                                     # Q2 if top left individual does not exist
      top_left_ind = c(top_left_ind,sample(step_moves[3:5],1))          # VALUE IF TRUE move down, right or diagonal
    }else{                                     
      top_left_ind = cbind(top_left_ind,sample(step_moves[3:5],         # VALUE IF FALSE
                                               dim(top_left_ind)[1],r=T))
    }
    ec_indiv_new = rbind(ec_indiv_new,top_left_ind)
  }
  if(sum(ec_indiv[,1] == hab_dim-dim+1)>0){   #top right
    top_right_ind = ec_indiv[ec_indiv[,1]==hab_dim-dim+1,]
    if(is.null(dim(top_right_ind))){
      top_right_ind = c(top_right_ind,sample(step_moves[5:7],1))
    }else{
      top_right_ind = cbind(top_right_ind,sample(step_moves[5:7],
                                                 dim(top_right_ind)[1],r=T))
    }
    ec_indiv_new = rbind(ec_indiv_new,top_right_ind)
  }
  if(sum(ec_indiv[,1] == dim)>0){              #bottom left
    bottom_left_ind = ec_indiv[ec_indiv[,1]==dim,]
    if(is.null(dim(bottom_left_ind))){
      bottom_left_ind = c(bottom_left_ind,sample(step_moves[1:3],1))
    }else{
      bottom_left_ind = cbind(bottom_left_ind,sample(step_moves[1:3],
                                                     dim(bottom_left_ind)[1],r=T))
    }
    ec_indiv_new = rbind(ec_indiv_new,bottom_left_ind)
  }
  if(sum(ec_indiv[,1] == hab_dim)>0){           #bottom righ
    bottom_right_ind = ec_indiv[ec_indiv[,1]==hab_dim,]
    if(is.null(dim(bottom_right_ind))){
      bottom_right_ind = c(bottom_right_ind,sample(step_moves[c(1,7:8)],1))
    }else{
      bottom_right_ind = cbind(bottom_right_ind,sample(step_moves[c(1,7:8)],
                                                       dim(bottom_right_ind)[1],r=T))
    }
    ec_indiv_new = rbind(ec_indiv_new,bottom_right_ind)
  }
  if(sum(ec_indiv[,1]%in%(2:(dim-1)))>0){        #left edge
    left_edge_ind = ec_indiv[ec_indiv[,1]%in%(2:(dim-1)),]
    if(is.null(dim(left_edge_ind))){
      left_edge_ind = c(left_edge_ind,sample(step_moves[1:5],1))
    }else{
      left_edge_ind = cbind(left_edge_ind,sample(step_moves[1:5],
                                                 dim(left_edge_ind)[1],r=T))
    }
    ec_indiv_new = rbind(ec_indiv_new,left_edge_ind)
  }
  if(sum(ec_indiv[,1]%in%((hab_dim-dim+2):(hab_dim-1)))>0){    #right edge
    right_edge_ind = ec_indiv[ec_indiv[,1]%in%((hab_dim-dim+2):(hab_dim-1)),]
    if(is.null(dim(right_edge_ind))){
      right_edge_ind = c(right_edge_ind,sample(step_moves[c(1,5:8)],1))
    }else{
      right_edge_ind = cbind(right_edge_ind,sample(step_moves[c(1,5:8)],
                                                   dim(right_edge_ind)[1],r=T))
    }
    ec_indiv_new = rbind(ec_indiv_new,right_edge_ind)
  }
  if(sum(ec_indiv[,1]%in%setdiff(which(1:hab_dim%%dim==1),c(1,hab_dim-dim+1)))>0){ #top edge
    top_edge_ind = ec_indiv[ec_indiv[,1]%in%setdiff(which(1:hab_dim%%dim==1),
                                                    c(1,hab_dim-dim+1)),]
    if(is.null(dim(top_edge_ind))){
      top_edge_ind = c(top_edge_ind,sample(step_moves[3:7],1))
    }else{
      top_edge_ind = cbind(top_edge_ind,sample(step_moves[3:7],
                                               dim(top_edge_ind)[1],r=T))
    }
    ec_indiv_new = rbind(ec_indiv_new,top_edge_ind)
  }
  if(sum(ec_indiv[,1]%in%setdiff(which(1:hab_dim%%dim==0),c(dim,hab_dim)))>0){ #bottom edge
    bottom_edge_ind = ec_indiv[ec_indiv[,1]%in%setdiff(which(1:hab_dim%%dim==0),
                                                       c(dim,hab_dim)),]
    if(is.null(dim(bottom_edge_ind))){
      bottom_edge_ind = c(bottom_edge_ind,sample(step_moves[c(1:3,7:8)],1))
    }else{
      bottom_edge_ind = cbind(bottom_edge_ind,sample(step_moves[c(1:3,7:8)],
                                                     dim(bottom_edge_ind)[1],r=T))
    }
    ec_indiv_new = rbind(ec_indiv_new,bottom_edge_ind)
    
  }
  
  # ------------------------------------- move -------------------------------------
  rownames(ec_indiv_new) = NULL
  
  locs_new = rbind(offec_indiv,ec_indiv_new)
  locs_new = locs_new[order(locs_new[,1]),]
  colnames(locs_new)[5]="move"
  
  #Move
  locs_new[,1] = locs_new[,1] + locs_new[,5]                              # add move to cell value
  locs_new[,3] = hab_vals[locs_new[,1]]                                   # get hab vals for new cells
  locs_new = locs_new[order(locs_new[,1]),]  
  
  return(locs_new)                                                        # Should now be merged with `locs`? but then all subsequent locs_new will need to be changed
  #here might be a good place to add a loop recording position for each timestep so we can simulate their movement over time?
}

#' =======================================================================
#' #3b. Function - fight
#' =======================================================================

#WORKS & minimised - not read for params

fight <- function(locs_new){
  
  #Interacting species are those occupying the same cell -----------------
  
  loc_ind <- locs_new[,1]                                              # extracts cell locations
  int_loc <- loc_ind[duplicated(loc_ind)]                              # extracts duplicated values from new locations (ie where individuals meet)
  fighters <- locs_new[locs_new[,1]%in%int_loc,]                       # extracts values for fighters from locs_new
  #should be an option here to avoid fight if eneregy low?
  fighters[,4] = fighters[,4]-fight_eloss                              # Aggression is energetically expensive
  locs_new = locs_new[locs_new[,1]%in%setdiff(locs_new[,1],fighters[,1]),]  # Remove fighters from set of individuals 
  
  #Roulette selection for winners ----------------------------------------
  all_winners = NULL
  
  for(i in unique(fighters[,1])){
    sp_in_cell = fighters[fighters[,1]==i,]                              # extracts info about only fighting species in this loop
    sp_ids = sp_in_cell[,2]                                              # extracts species number
    nosp = length(sp_ids)                                               
    win_vect = rep(0,nosp)                                               # creates vector for the number of fighting species 
    
    #Fights occur in random pairs -------------------------------------------
    #Losers removed until one individual remains 
    for(j in 1:(nosp-1)){
      fightj = sp_in_cell[sample(1:dim(sp_in_cell)[1],2),]               # extracts out rows of two random individs in same cell
      
      aggress_vals = c(aggression[fightj[1,2],fightj[2,2]],              # extracts out aggression level based on species number
                       aggression[fightj[2,2],fightj[1,2]])              # depending on which is x,y location in aggresion matrix? don't understand this matrix...
      
      energy_levels = c(fightj[1,4],fightj[2,4])                         # extracts energy levels of both fighters
      
      status = aggress_vals*energy_levels                                # Status of individual is product of aggression and energy level. 
      # Highly aggressive but lower energy individuals can lose to less aggressive...
      # individuals who have high energy
      
      out_prob = cumsum(status^exp_fight)/sum(status^exp_fight)          # cumulative sum - not sure what's going on here. Because cumsum, last umber is always coing to be equal to 
      # sum, so =1 ??
      
      if(out_prob[1]>runif(1)){                                          # roulette? runif() determines win vs lose?
        loser = fightj[2,]
      }else{
        loser = fightj[1,]
      }
      
      rowcount = 0
      rowsum = 0
      while(rowsum!=6){                                                  # Don't know where 6 comes from?
        rowcount = rowcount + 1
        rowsum = sum(sp_in_cell[rowcount,] == loser)                     # identifies individ randomly allocated as loser
      }
      sp_in_cell = sp_in_cell[-rowcount,]                                # removes loser, then loops back for next fight
    }
    
    initial = fighters[fighters[,1]==i,]                                 # identifies next cell with multiple individs, extracts only their info 
    
    rowcount = 0
    rowsum = 0
    while(rowsum!=6){                                                    # not sure...
      rowcount = rowcount + 1
      rowsum = sum(initial[rowcount,] == sp_in_cell)
    }
    win_vect[rowcount] = 1                                               # changes 0 to 1 in win vs lose column
    all_winners = c(all_winners,win_vect)
  }
  fighters = cbind(fighters,all_winners)                               # merges 
  
  #Lose energy
  
  winners = fighters[fighters[,7]==1,]                                   # Differentiate winners and losers 
  losers = fighters[fighters[,7]==0,]
  
  losers = losers[losers[,4]>0,]                                         # deletes dead individuals (energy =0)

  # Combines winner with those who didn't engage in aggression (i.e. were sole occupants of cell)
  locs_new = cbind(locs_new,all_winners=1)
  locs_new = rbind(locs_new,winners)
  locs_new = locs_new[order(locs_new[,1]),]
  
    return(locs, locs_new)
}

  

#' =======================================================================
#' #3c. Function - Feed, energy loss, die
#' =======================================================================

die <- function(locs){
  locs[,4] = locs[,4]-eloss+locs[,3]*100;                               # -energy from one step, + energy from habitat (feed) #order??
  locs[locs[,4]>100,4] = 100                                            # make sure does not exceed 100
  mort = as.logical(rbinom(dim(locs)[1],1,a/(1+b*exp(-v*locs[,4]))))    # Death = Bernoulli trial w probability logistically dependent on current energy
  locs = locs[mort,]                                                    # removes individuals from `locs` picked in `mort`
  return(locs)
}

#' =======================================================================
#' #3d. Reproduce
#' =======================================================================

reproduce <- function(locs){

#Low-energy individuals (I sound like Trump!) are unlikely to reproduce.
locs_new = cbind(locs_new,rbinom(dim(locs_new)[1],1,a_rep/(1+b*exp(-v*locs_new[,4])))) # Reproduction is Bernouilli trial dependent on current energy levels. 
offspring = locs_new[locs_new[,8]==1,]                                 # extracts offspring details (same as parent)
offspring[,4] = offspring[,4]/offspring_pen                            # Set offspring energy to 50% of parents to account for lower chance of survival

locs_new[,4] <- ifelse(locs_new[,8]==1,                                # if reproducing
                       locs_new[,4]-repro_pen,                         # remove energy
                       locs_new[,4])                                   # non-breeders no change

locs = rbind(locs_new[,1:4],offspring[,1:4],losers[,1:4])              # merges output with new locs
locs = locs[order(locs[,1]),]
return(locs)
}

#' =======================================================================
#' #4 Simulate
#' =======================================================================


loop = 0
start_time = proc.time()
for(reps in 1:replicates){
  if(reps<=replicates/2){                                                  # Aggression relationships (transitive or intransitive)  
    intransitive = T 
  }else{
    intransitive = F
  }
  
  #' Define aggression values
  if(intransitive){
    aggression = matrix(0.5-diff_ag + runif(nspecies*nspecies)*(diff_ag*2),ncol=nspecies)
    diag(aggression) = 0.5 # For interaction between individuals of the same species
  }else{
    aggression = matrix(ncol=nspecies,nrow=nspecies,
                        sort(0.5-diff_ag + runif(nspecies*nspecies)*(diff_ag*2),d=T))
  }
  
  
  
  
  
  
  loop=0  
  for(generation in 1:ngenerations){
    loop = loop + 1
    
    #Progress bar update
    pb$tick()
    
    ## 1. Feed, energy gain, energy loss, die (from energy loss) 
    locs2 <- die(locs=locs)
    
    ## 2. Move
    locs3 <- move(locs=locs2, hab_vals=hab_vals)
    
    ## 3. Fight, feed, energy gain, loss, reproduce die
    locs4 <- fight(locs_new = locs3, hab_vals=hab_vals)   
    
    ## 4. update sp richness
    rich[rich$gen==generation & rich$rep==reps,3]=length(unique(locs[,2])) #updates species richness
    rich[rich$gen==generation & rich$rep==reps,4]=intransitive
    
    close(pb)
    
    ## 5. plot & outputs
    
    labs = c('Intransitive','Transitive')
    count = 0
    for(type in c(T,F)){
      count = count + 1
      subData = rich[rich$type==type,]
      
      frame = data.frame(gen=1:ngenerations,
                         mean=tapply(subData$rich,subData$gen,mean),
                         se=tapply(subData$rich,subData$gen,function(x) sd(x)/sqrt(length(x))),
                         type=labs[count])
      
      assign(labs[count],frame)
      
    }
    
    
    allRich = rbind(Intransitive,Transitive)
    
    rich.plot = ggplot(allRich) + theme_bw() +
      geom_ribbon(aes(x=gen,ymin=mean-se,ymax=mean+se,fill=type),alpha=0.8) +
      scale_y_continuous(limits=c(0,nspecies)) +
      labs(y='Species richness',x='Generation') +
      scale_fill_manual(values=c("#35B779FF","#440154FF"),
                        name='Aggression type',
                        guide=guide_legend(keywidth = 1,
                                           keyheight=0.5)) +
      theme(legend.position='top',
            panel.grid=element_line(size=0.2),
            axis.title = element_text(size=8),
            axis.text = element_text(size=7),
            legend.title = element_text(size=8),
            legend.text = element_text(size=7))
    
    cairo_pdf('aggression.pdf',width=5,height=5.2)
    print(rich.plot)
    dev.off()
    
  }
  
  
  
  
  #' =======================================================================
  #' #?. Not sure where this fits in
  #' =======================================================================
  #
  #'  
  
  
  
  
  
  
  
  
  
  
  close(pb)
  
  labs = c('Intransitive','Transitive')
  count = 0
  for(type in c(T,F)){
    count = count + 1
    subData = rich[rich$type==type,]
    
    frame = data.frame(gen=1:ngenerations,
                       mean=tapply(subData$rich,subData$gen,mean),
                       se=tapply(subData$rich,subData$gen,function(x) sd(x)/sqrt(length(x))),
                       type=labs[count])
    
    assign(labs[count],frame)
    
  }
  
  
  allRich = rbind(Intransitive,Transitive)
  
  rich.plot = ggplot(allRich) + theme_bw() +
    geom_ribbon(aes(x=gen,ymin=mean-se,ymax=mean+se,fill=type),alpha=0.8) +
    scale_y_continuous(limits=c(0,nspecies)) +
    labs(y='Species richness',x='Generation') +
    scale_fill_manual(values=c("#35B779FF","#440154FF"),
                      name='Aggression type',
                      guide=guide_legend(keywidth = 1,
                                         keyheight=0.5)) +
    theme(legend.position='top',
          panel.grid=element_line(size=0.2),
          axis.title = element_text(size=8),
          axis.text = element_text(size=7),
          legend.title = element_text(size=8),
          legend.text = element_text(size=7))
  
  cairo_pdf('aggression.pdf',width=5,height=5.2)
  print(rich.plot)
  dev.off()
  
  
  #'
  #' ## Aggress ------------------------------------------------------------
  
  aggression = matrix(0.5-diff_ag + runif(nspecies*nspecies)*(diff_ag*2),ncol=nspecies)
  diag(aggression) = 0.5 # For interaction between individuals of the same species 
  aggression = matrix(ncol=nspecies,nrow=nspecies,
                      sort(0.5-diff_ag + runif(nspecies*nspecies)*(diff_ag*2),d=T))
}



# Aggression relationships (transitive or intransitive) ##Gareth working on this? ##
#' 
if(reps<=replicates/2){
  intransitive = T 
}else{
  intransitive = F
}

#' Define aggression values
if(intransitive){
  aggression = matrix(0.5-diff_ag + runif(nspecies*nspecies)*(diff_ag*2),ncol=nspecies)
  diag(aggression) = 0.5 # For interaction between individuals of the same species 
  aggression = matrix(ncol=nspecies,nrow=nspecies,
                      sort(0.5-diff_ag + runif(nspecies*nspecies)*(diff_ag*2),d=T))
}


#deleted from aggresion function
mort = as.logical(rbinom(dim(losers)[1],1,a/(1+b*exp(-v*losers[,4])))) # randomly selects some of remianing indivuals to die
losers = losers[mort,]                                                 # deletes those individs



#These individuals lose energy, consume, die and reproduce
locs_new[,4] = locs_new[,4]-eloss+locs_new[,3]*100                     # current energy= before energy - energy loss (~aging?) + energy gain (feeding)
locs_new[locs_new[,4]>100,4] = 100                                     # reset energy above 100 to 100
locs_new = locs_new[locs_new[,4]>0,]                                   # deletes individs with energy 0
mort = as.logical(rbinom(dim(locs_new)[1],1,a/(1+b*exp(-v*locs_new[,4])))) # randomly selects some of remaining indivuals to die
locs_new = locs_new[mort,]

#Low-energy individuals (I sound like Trump!) are unlikely to reproduce.
locs_new = cbind(locs_new,rbinom(dim(locs_new)[1],1,a_rep/(1+b*exp(-v*locs_new[,4])))) # Reproduction is Bernouilli trial dependent on current energy levels. 
offspring = locs_new[locs_new[,8]==1,]                                 # extracts offspring details (same as parent)
offspring[,4] = offspring[,4]/offspring_pen                            # Set offspring energy to 50% of parents to account for lower chance of survival

locs = rbind(locs_new[,1:4],offspring[,1:4],losers[,1:4])              # merges output with new locs
locs = locs[order(locs[,1]),]

habitat = matrix(0,ncol=dim,nrow=dim)                                  # creates new habitat matrix
habitat[locs[,1]] = 1                                                  # with individuals based on locations in `locs`