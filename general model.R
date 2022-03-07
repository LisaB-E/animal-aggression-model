#' ---
#' title: “Angry fish - linking changes in aggression after habitat loss to demographic changes at the population level using IBM "
#' author: “LBE"
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' output: github_document
#' ---

  

#+ echo = F
#---- 1. WORKSPACE PREP --------------------------------------------------  
#' # 1. Workspace preparation 

rm(list=ls()) 


#' Load packages
#+ packages, message=F, results='hide'
packages = c('tidyverse', 'reshape','gstat','ggplot2',
             'viridis','tictoc', 'progress',
             'dplyr', 'gganimate')
load.pack = lapply(packages,require,char=T)
load.pack

#+ echo = F 
#----- 2. PARAMETERS -----------------------------------------------------  
#' # 2. Parameters  

#' IBM Parameters *should include constraints on each parametre*. Double check how many are needed in reduced model

#simulation
ngenerations  = 100    # No. generations
replicates    = 5     # No. replicates (first half trans. second half intrans.)

#initialisation
dim           = 100   # dimension of square habitat array
hab_dim       = dim^2 # total no. cells
nspecies      = 30    # No. species
nindiv        = 40    # No. individuals per species
tot_indiv     = nspecies*nindiv # Total individuals

#energy
#e-value of repro>feed>aggro>somatic>


e_feed        = 100    # energy gain from habitat factor
eloss         = e_feed * 0.98     #[DEB] # Time-step energy loss (ie somatic growth)
fight_eloss   = 20    # Energy loss from aggression
offspring_pen = (4*10^-6)*e_feed  #[DEB]   # factor by which to scale offspring energy 
repro_pen     = 0.02*e_feed   #[DEB]  # energy loss from reproduction


#aggression
exp_fight     = 20    # Roulette selection exponent
diff_ag       = 0.05  # Differential aggression \in [0,0.5)

#death
a             = 1
b             = 100
v             = 0.08  # Logistic mortality parameters

#repro
a_rep         = 0.4   # Asymptotic reproduction probability



#' Potential one time-step movement. Can move to any adjacent cell (inc diagonal) 
step_moves = c(0,-1,dim-1,dim,dim+1,1,-dim+1,-dim,-dim-1)      #added option to stay in the same place

#' Potential aggression matrix - M1 assume all species have equal competitive potential - so outcome will only depend on energy levels
aggression = matrix(0.5, ncol=nspecies, nrow = nspecies)
diag(aggression) = 0.5 # For interaction between individuals of the same species 

#' Store species richness at each generation IS THIS NEEDED, CAN'T GET FROM END LISTS?
rich = data.frame(gen=rep(1:ngenerations,replicates),
                  rep=rep(1:replicates,each=ngenerations),
                  rich=NA,
                  type=NA)

#' Progress bar 
#+ results='hide'
pb <- txtProgressBar(min = 0, max = replicates*ngenerations, style = 3)

#' unique ID's for each individual
IDs <- c(1:(nspecies*nindiv*ngenerations*replicates)) #does it matter if this is unique per replicate, or unique throughout? right now unique throughout

#+echo = F
#----- 3. INITIALISE  ----------------------------------------------------

#' # 3. Initialise - habitat  

#' Simulate habitat 
  hab_grid = expand.grid(1:dim, 1:dim)                                     # The habitat arena
  names(hab_grid) = c('x','y')
  hab_bin = sample(x=c(1,0), size = hab_dim, replace = T)                  # samples binary habitat values
  hab_bin = cbind(hab_grid, hab_bin)
  
#' Make matrix
  hab_vals = as.matrix(cast(hab_bin,
                            x~y,
                            value='hab_bin'))
 
#' Plot habitat values
#+ habitat.matrix
  hab.plot = ggplot(hab_bin) + theme_bw() +
    geom_raster(aes(x,y,fill=hab_bin)) +
    scale_x_continuous(expand=expansion(add=0)) +
    scale_y_continuous(expand=expansion(add=0)) +
    theme(legend.position = 'none',
          axis.text=element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())+
    scale_fill_gradient(low = "gray50", high = "gray20")
 hab.plot

#' Randomly place individuals of each species in habitat
locs = cbind(sample(hab_dim,tot_indiv),rep(1:nspecies,each=nindiv))       # Randomly picks locations (V1-cell) for all individuals (V2) 
locs = cbind(locs,hab_vals[locs[,1]])                                     # habitat value

#' Add energetic value
locs = cbind(locs,rep(100,nindiv))                                        # starting energetic value =100?
locs = locs[order(locs[,1]),] 

locs <- cbind(locs, IDs[1:nrow(locs)])
IDs <- setdiff(IDs,locs[,5])
colnames(locs)=c("loc", "sp", "hab_val", "e_val", "ID")

OG_locs <- locs                                                           # each rep in sim starts with same initialised habitat 



#+ echo=F
#----- 4. FUNCTIONS ---------------------------------------------------------  
#' # 4. Functions
#' ### Function - move  

move <- function(){      
  
# Place individuals in habitat
  habitat = matrix(0,ncol=dim,nrow=dim)                                # make 100 x 100 matrix
  habitat[locs[,1]] = 1                                                # places individual in the matrix
  
# Find edge and corner individuals
  edge_corner <- unique(c(1:dim,                                       # top row
                          which(1:hab_dim%%dim%in%0:1),                # sides
                          (hab_dim-dim+1):hab_dim))                    # bottom
  
  
ec_occ      <- edge_corner[which(habitat[edge_corner]==1)]             # which edge corner cells are occupied
ec_indiv    <- locs[locs[,1]%in%ec_occ, , drop=FALSE]                  # which individuals are on edge corner cells (V2&V1)
  offec_indiv <- locs[locs[,1]%in%setdiff(locs[,1],ec_indiv[,1]),]     # individuals  NOT on edge corner cells (V1&v2)
  
  offec_indiv <- cbind(offec_indiv, sample(step_moves, dim(offec_indiv)[1], replace = TRUE)) # add move step to non edge inds (V3)
  
# Move edge and corner individuals
#goal: figure out how to do this neater (not sure the if else statements are neede for samppling of the step_moves, surely a dim=1 woudl replace the if TRUE part?)*
  
  ec_indiv_new <- NULL
  if(sum(ec_indiv[,1]==1)>0){                                          # Q1 find if any inds are in top left corner
    top_left_ind = ec_indiv[ec_indiv[,1]==1,]                          # VALUE IF TRUE -pull out that data
    if(is.null(dim(top_left_ind))){                                    # Q2 if top left individual does not exist
      top_left_ind = c(top_left_ind,sample(step_moves[c(1,4:6)],1))    # VALUE IF TRUE move down, right or diagonal
    }else{                                     
      top_left_ind = cbind(top_left_ind,sample(step_moves[c(1,4:6)],   # VALUE IF FALSE
                                               dim(top_left_ind)[1],r=T))
    }
    ec_indiv_new = rbind(ec_indiv_new,top_left_ind)
  }
  if(sum(ec_indiv[,1] == hab_dim-dim+1)>0){   #top right               
    top_right_ind = ec_indiv[ec_indiv[,1]==hab_dim-dim+1,]
    if(is.null(dim(top_right_ind))){
      top_right_ind = c(top_right_ind,sample(step_moves[c(1,6:8)],1))
    }else{
      top_right_ind = cbind(top_right_ind,sample(step_moves[c(1,6:8)],
                                                 dim(top_right_ind)[1],r=T))
    }
    ec_indiv_new = rbind(ec_indiv_new,top_right_ind)
  }
  if(sum(ec_indiv[,1] == dim)>0){                                                     # bottom left
    bottom_left_ind = ec_indiv[ec_indiv[,1]==dim,]
    if(is.null(dim(bottom_left_ind))){
      bottom_left_ind = c(bottom_left_ind,sample(step_moves[c(1, 2:4)],1))
    }else{
      bottom_left_ind = cbind(bottom_left_ind,sample(step_moves[c(1,2:4)],
                                                     dim(bottom_left_ind)[1],r=T))
    }
    ec_indiv_new = rbind(ec_indiv_new,bottom_left_ind)
  }

  if(sum(ec_indiv[,1] == hab_dim)>0){                                                  # bottom right
    bottom_right_ind = ec_indiv[ec_indiv[,1]==hab_dim,]
    if(is.null(dim(bottom_right_ind))){
      bottom_right_ind = c(bottom_right_ind,sample(step_moves[c(1:2,8:9)],1))
    }else{
      bottom_right_ind = cbind(bottom_right_ind,sample(step_moves[c(1:2,8:9)],
                                                       dim(bottom_right_ind)[1],r=T))
    }
    ec_indiv_new = rbind(ec_indiv_new,bottom_right_ind)
  }
  if(sum(ec_indiv[,1]%in%(2:(dim-1)))>0){                                              # left edge
    left_edge_ind = ec_indiv[ec_indiv[,1]%in%(2:(dim-1)),]
    if(is.null(dim(left_edge_ind))){
      left_edge_ind = c(left_edge_ind,sample(step_moves[c(1,2:6)],1))
    }else{
      left_edge_ind = cbind(left_edge_ind,sample(step_moves[c(1,2:6)],
                                                 dim(left_edge_ind)[1],r=T))
    }
    ec_indiv_new = rbind(ec_indiv_new,left_edge_ind)
  }

  if(sum(ec_indiv[,1]%in%((hab_dim-dim+2):(hab_dim-1)))>0){                            # right edge

    right_edge_ind = ec_indiv[ec_indiv[,1]%in%((hab_dim-dim+2):(hab_dim-1)),]
    if(is.null(dim(right_edge_ind))){
      right_edge_ind = c(right_edge_ind,sample(step_moves[c(1:2,6:9)],1))
    }else{
      right_edge_ind = cbind(right_edge_ind,sample(step_moves[c(1:2,6:9)],
                                                   dim(right_edge_ind)[1],r=T))
    }
    ec_indiv_new = rbind(ec_indiv_new,right_edge_ind)
  }

  if(sum(ec_indiv[,1]%in%setdiff(which(1:hab_dim%%dim==1),c(1,hab_dim-dim+1)))>0){      # top edge
   top_edge_ind = ec_indiv[ec_indiv[,1]%in%setdiff(which(1:hab_dim%%dim==1),
                                                    c(1,hab_dim-dim+1)),]
    if(is.null(dim(top_edge_ind))){
      top_edge_ind = c(top_edge_ind,sample(step_moves[c(1, 4:8)],1))
    }else{
      top_edge_ind = cbind(top_edge_ind,sample(step_moves[c(4:8)],
                                               dim(top_edge_ind)[1],r=T))
    }
    ec_indiv_new = rbind(ec_indiv_new,top_edge_ind)
  }

  if(sum(ec_indiv[,1]%in%setdiff(which(1:hab_dim%%dim==0),c(dim,hab_dim)))>0){          # bottom edge

    bottom_edge_ind = ec_indiv[ec_indiv[,1]%in%setdiff(which(1:hab_dim%%dim==0),
                                                       c(dim,hab_dim)),]
    if(is.null(dim(bottom_edge_ind))){
      bottom_edge_ind = c(bottom_edge_ind,sample(step_moves[c(1:4,8:9)],1))
    }else{
      bottom_edge_ind = cbind(bottom_edge_ind,sample(step_moves[c(1:4,8:9)],
                                                     dim(bottom_edge_ind)[1],r=T))
    }
    ec_indiv_new = rbind(ec_indiv_new,bottom_edge_ind)
    

  }                                                                                     # end find edge individs section

# move 
rownames(ec_indiv_new) = NULL

locs_new = rbind(offec_indiv,ec_indiv_new)
locs_new = locs_new[order(locs_new[,1]),]
locs_new[6:7] <- as.numeric(locs_new[6:7])

# Move
locs_new[,1] = locs_new[,1] + locs_new[,6]                            # add move to cell value
locs_new[,3] = hab_vals[locs_new[,1]]                                 # get hab vals for new cells
locs_new = locs_new[order(locs_new[,1]),]  

gen <- rep(generation, nrow(locs_new))                                # output cleanup 
repl <- rep(replic, nrow(locs_new)) 
locs_new <- cbind(locs_new, gen, repl)                                # add generation number V6, repl V7

# Getting x & y coords for animation
temp_habitat1 <- matrix(0,ncol=dim,nrow=dim)                          # make 100 x 100 matrix
temp_habitat1[locs_new[,1]] = 1                                       # place inds on habitat, ignores later dupes
temp_coords1 <- which(temp_habitat1>0, arr.ind=TRUE)                  # get x & y coords

triple_threat_list <- list()                                          # get duplicates, including triples etc
temp_dupes <- matrix(locs_new[,1][duplicated(locs_new[,1], fromLast = TRUE)])
  
for (dupe in 1:nrow(temp_dupes)) {
  if(nrow(temp_dupes)==0) {                                           #for cases where there are no duplicates
    triple_threat_list <- list()
  } else {
  temp_habitat <- matrix(0,ncol=dim,nrow=dim)
  temp_habitat[temp_dupes[dupe,1 ]] = 1
  temp_coordsX <- which(temp_habitat>0, arr.ind=TRUE)
  triple_threat_list[[dupe]] <- temp_coordsX
}}



temp_coords2 <- do.call(rbind,triple_threat_list)                      # extract each list
temp_coords <- rbind(temp_coords1, temp_coords2)                       # merge coords inclduing dupes
temp_coords <- temp_coords[order(temp_coords[,2], temp_coords[,1]),]   # reorder

locs_new <- cbind(locs_new, temp_coords)

return(locs_new)                                                        
  
}

  
#' ### Function - fight

fight <- function(){

  locs_new <- locs                                                     # pulls from move endpoint in sim            
  # Interacting species are those occupying the same cell 
  loc_ind <- locs_new[,1]                                              # extracts cell locations
  int_loc <- loc_ind[duplicated(loc_ind)]                              # extracts cell contains multiple individuals ('fighters')
  fighters <- locs_new[locs_new[,1]%in%int_loc,]                       # extracts values for fighters

  locs_new = locs_new[locs_new[,1]%in%setdiff(locs_new[,1],fighters[,1]),]             # Remove fighters from set of individuals            #BEWARE!!! if dplyr is loaded - changes results! hould chekc this forst! NEED to fix
  
  fighters[,4] = fighters[,4]-fight_eloss                              # Aggression is energetically expensive
  # *DOES THIS MAKE SENSE? IT'S ASSUMING THEY MAKE THIS DECISION AHEAD OF TIME, I WILL OR WILL NOT ENGAGE IN  AGGRESISON BASED ON ENERGY LEVELS, OR SHOULD THIS BE AFTER??*
 
  # Roulette selection for winners
  all_winners = NULL

  # fight in pairs 
  for(i in unique(fighters[,1])){
    sp_in_cell <- fighters[fighters[,1]==i,]                             # extracts info about fighters in same cell
    win_vect = rep(0,length(sp_in_cell[,2]))                             # creates vector for the number of fighting indiviuals in cell
    
# Fights occur in random pairs 
    while(nrow(sp_in_cell)>1) {                                          # Losers removed until one individual remains 
      fightj = sp_in_cell[sample(1:dim(sp_in_cell)[1],2),]               # extracts out rows of two random individs in same cell
       if (sum(fightj[,4]>0)!=2){                                        # dead (e<0) inds lose straight away
          low_e_loser <- fightj[fightj[,4] <= 0,,drop = F]
          sp_in_cell <- sp_in_cell[!sp_in_cell[,5] %in% low_e_loser[,5],,drop=F]
        } else {
        
      aggress_vals = c(aggression[fightj[1,2],fightj[2,2]],              # extracts out aggression level based on species number
                       aggression[fightj[2,2],fightj[1,2]])              # depending on which is x,y location in aggresion matrix? don't understand this matrix...
      
      energy_levels = c(fightj[1,4],fightj[2,4])                         # extracts energy levels of both fighters
      
      status = aggress_vals*energy_levels                                # Status of individual is product of aggression and energy level. 
      # Highly aggressive but lower energy individuals can lose to less aggressive...
      # individuals who have high energy
      
      out_prob = cumsum(status^exp_fight)/sum(status^exp_fight)          # cumulative sum - not sure what's going on here. Because cumsum, last umber is always coing to be equal to 
      # sum, so =1 ??
      
      if(out_prob[1]>runif(1)){                                          # roulette? runif() determines win vs lose? Why involve the aggression values and energy levels if then outcome is ultimately determined by randomly generated number?
        loser = fightj[2,,drop=F]
      }else{
        loser = fightj[1,,drop=F]
      }
       
# remove loser from cell
      
      rowcount = 0
      rowsum = 0
                                                          
      while(rowsum!=5){                                                  # 5 = all cells matches with the loser
        rowcount = rowcount+1
        rowsum = sum(sp_in_cell[rowcount,] == loser)                     # checks whether top row = loser
      }
       sp_in_cell = sp_in_cell[-rowcount,,drop=F]                        # removes loser from sp in cell
      }      
      }                                                                  # loops back to next fightJ
    
    initial = fighters[fighters[,1]==i,]                                 # identifies starting individs in cell (.ie OG sp_in_cell)
   
# identify which individ is winner to add to win_vect tally
   if(nrow(sp_in_cell)>0) {                                              # if all were low energy
     rowcount = 0
    rowsum = 0
    while(rowsum!=5){                                 
      rowcount = rowcount + 1
      rowsum = sum(initial[rowcount,] == sp_in_cell)                     # finds winner in sp_in_cell
    }
    win_vect[rowcount] = 1                                               # changes 0 to 1 in win vs lose column
   }
    all_winners = c(all_winners,win_vect)

    }                                                                    # loops back to next cell with multiples
  fighters = cbind(fighters,all_winners)                                 # merges 
  
# Lose energy post fight
  fighters[,4] = fighters[,4]-fight_eloss                                # Aggression is energetically expensive DOES THIS MAKE SENSE? IT'S ASSUMING THEY MAKE THIS DECISION AHEAD OF TIME, I WILL OR WILL NOT ENGAGE  IN AGGRESISON BASED ON ENERGY LEVELS, OR SHOULD THIS BE AFTER?? also, why arwe there two of these, once before and once after?


  winners = fighters[fighters[,6]==1, , drop=FALSE]                      # Differentiate winners and losers 
  losers = fighters[fighters[,6]==0, ,drop=FALSE]
  
  losers = losers[losers[,4]>0, ,drop=FALSE]                             # deletes dead individuals (energy =0). I wonder if there should alos be random mortality in these guys?

  # Combines winner with those who didn't engage in aggression (i.e. were sole occupants of cell)
  locs_new = cbind(locs_new,all_winners=1)
  locs_new = rbind(locs_new,winners)
  
  locs_new = locs_new[order(locs_new[,1]),]
  
  locs_new[,4] = locs_new[,4]-eloss+locs_new[,3]*100                     # winners -e from engaging? +e by feeding. Check whetehr this doesn't already occur up top?
  locs_new[locs_new[,4]>100,4] = 100                                     # e > 100 = 100
  locs_new = locs_new[locs_new[,4]>0,]                                   # e < 0 = dead
  
# Death 
  mort = as.logical(rbinom(dim(locs_new)[1],1,a/
                             (1+b*exp(-v*locs_new[,4]))))                # Death = Bernoulli trial w probability logistically dependent on current energy
  locs_new = locs_new[mort,]
  
# Reproduction 
  locs_new = cbind(locs_new,rbinom(dim(locs_new)[1],1,a_rep/
                              (1+b*exp(-v*locs_new[,4]))))               # Reproduction is Bernouilli trail dependent on current energy levels. READ INTO THIS
  offspring = locs_new[locs_new[,7]==1, , drop=F]
  offspring[,5] = IDs[1:nrow(offspring)]                                 # give unique ID
  offspring[,4] = offspring[,4]*offspring_pen                            # Set offspring energy to 50% of parents [CHANGED TO DEB THEORY, * instead of /]
  
  locs_new = rbind(locs_new[,1:5],offspring[,1:5],losers[,1:5])          # Combines all after fight, feed, death, reproduction
  locs_new = locs_new[order(locs_new[,1]),]

# save data
  gen <- rep(generation, nrow(locs_new))
  # Store gen (WARNING - change 1 to generation, replic)
  repl <- rep(replic, nrow(locs_new))
  locs_new <- cbind(locs_new, gen,repl)  
  
    return(locs_new)
  
  }

#----- 5. SIMULATE -----------------------------------------------------------------  
#' # 5. Simulate
#+ results = 'hide'
loop=0
list_move=list()                                                      
list_fight=list()
list_abund=list()
list_rich = list()
move_hist = tibble()


# start replication loop
for (replic in 1:replicates) { 
locs <- OG_locs #reset to starting point - is this right level of replication?

# start timestep loop
for(generation in 1:ngenerations){
  gen_rep <- paste("G",generation,"R",replic, sep = "_")
  loop = loop + 1
  
# *Move*
  after_move <- move()                        # after_move = new_locs
  locs = after_move[,1:5]                     # updates locs to endpoint (final timestep) 
  list_move[[gen_rep]] <- after_move          # stores moves per generation as list
# print(paste('move_',' gen', generation, ' rep', replic, ' DONE', sep = ""))
  
# *Fight*
  after_fight <- fight()                      # after_fight = locs_new
# print(paste('fight_','gen', generation, ', rep', replic, ' DONE', sep = ""))
  locs <- after_fight[,1:5]                   # updates locs to endpoint (final timestep) #IS THIS THE PROBLEM?!
  IDs <- setdiff(IDs,locs[,5])                # removes IDs from pool of available
  
  rich <- length(unique(locs[,2]))            # No species
  list_rich[[gen_rep]] <- rich
 
  locs.df <-data.frame(locs)                  # Abundance
  locs.df$sp <- as.factor(locs.df$sp)
  abund <- locs.df %>%
    group_by(sp) %>%
    tally() %>%
    mutate(gen = generation,
           rep = replic)#;
  list_abund[[gen_rep]] <- abund           # stores abundance after fight per generation as list (move down)
 
   #SOMETHING IS WRONG HERE -- THE FINAL LIST IS CUT DOWN TO ONLY THOSE SPECIES REMAINING AT GEN 100...
  
  list_fight[[gen_rep]] <- after_fight     # stores outcomes from fights per generation as list
  
  setTxtProgressBar(pb, loop)
   }                                          # end gen

#save move data #WHY IS THIS ONE THE ONLY ONE THAT GETS COLLECTED IN SIM??
move_temp <- data.frame(do.call(rbind,list_move)) %>% #unnest nested list
  mutate(sp = str_c("SP", sp, sep = " "),
         ID = as_factor(ID),
         U_ID = str_c(sp, ID, sep="_"))
move_hist <- bind_rows(move_temp, move_hist)          #adds each gen rep to each other #IS THIS THE PROBLEM HERE - ARE WE MERGING WEIRDLY?

#save fight data 


}                                             # end rep
close(pb)

#exploring model output

#number of individuals per gen & rep
hist_sp <- move_hist %>% 
  group_by(gen, repl) %>% 
  dplyr::summarise(n.inds=length(unique(U_ID)),
                   n.sp = length(unique(sp)))



#----- 6. VISUALISE --------------------------------------------------
#can get rid of U_ID?

# Movement ----
random_select = sample(x=move_hist$U_ID,size=100)                               # reduced set for plotting n individs
sub_movers <- subset(move_hist, U_ID %in% c(random_select))

#+ movement.plot
(movement_fig <- ggplot(data=hab_bin, aes(x=x,y=y)) +
  geom_raster(aes(fill=hab_bin), alpha=0.7) +
  scale_fill_gradient( low="grey10", high="grey5") +
  geom_point(data=sub_movers, aes(x=col, y=row, col=sp),size=2, inherit.aes = FALSE)+
  scale_colour_viridis(option="viridis", discrete = TRUE)+
  geom_path(data=sub_movers, aes(x=col, y=row, col=U_ID),size=1.5, lineend = "round", inherit.aes = FALSE, alpha = 0.8)+
  ggtitle('Individual movement')+
  theme_dark()+
  theme(title = element_text(colour = "white"),
        plot.background = element_rect(fill = "grey10"),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "grey20", size = 0.2),
        panel.grid.minor = element_line(color = "grey20", size = 0.2),
        axis.title = element_text(colour="grey30")
            ) +
  # theme_dark()+
  theme(legend.position="none")+
 # labs(title = "Generation {previous_state} of 20" )+
  coord_fixed())


# Sp. Abundance ----
sp_hist <- as_tibble(data.frame(do.call(rbind,list_abund))) %>% #unlist sim
  mutate(sp = str_c("SP", sp)) %>% 
  group_by(gen, sp) %>% 
  summarise(mean_abund = mean(n),
            SEM_abund = sd(n)/sqrt(length(n)),
            reps = length(n))

sp_hist_n <- as_tibble(data.frame(do.call(rbind,list_abund))) %>% #unlist sim
  mutate(sp = str_c("SP", sp)) 

#+ abundance.plot
sp_fig <- ggplot( data = sp_hist, aes(x=gen, y=mean_abund, color=sp)) +
 geom_line(data = sp_hist_n, aes(x=gen, y=n, group=sp), colour="gray", alpha = 0.4)+
  #species specific abundance
  geom_line() +
  ggtitle('Species abundance')+
  ylab("Abundance (by species)")+
  xlab("Time step")+
  scale_colour_viridis(option="viridis", discrete = TRUE) +
  theme_classic()+
  theme(panel.background = element_rect(fill='black'), 
        plot.background = element_rect(fill='black', color=NA),
        panel.grid.major = element_line(colour=NA),
        panel.grid.minor = element_line(color =NA),
        axis.line = element_line(colour = 'white', size = 0.8),
        axis.text = element_text(colour = 'white', family = "Arial", size = 12),
        axis.ticks = element_line(colour = "white"),
        axis.title = element_text(colour = "white", size = 16, family = "Arial"),
        plot.title = element_text(colour = "white")
  ) +
  scale_y_continuous( expand=c(0,0)) +
  scale_x_continuous( expand=c(0,0))+
  theme(legend.position = "none")

#' *all species persist here because over the replicates, there's always at least 2 runs where each species survives 'til the bitter end, even though it's extinct in most other runs.*
sp_fig
ggsave("test2.png")

# SAD ----------
# (contd'd) from above
# idea: animate this thorugh time using gganimate()

#' mean should probbaly not be mean abundance, but mean species per bin?

library(sads)
SAD_data <- as_tibble(subset(sp_hist, gen==ngenerations)) #extract final species (but as mean??)

#plotting SAD using final gen & final rep
SAD_data <- as_tibble(data.frame(do.call(rbind,list_abund))) %>% #unlist sim
  filter(gen==ngenerations, rep == replicates) %>% #  keeps final gen final rep , but figure out how to deal with a mean instead
  select(sp,n)
SAD_data_rep <- as.numeric(with(SAD_data, rep(sp,n))) #make into format for octav()
SAD_octave <- octav(SAD_data_rep)

plot(SAD_octave)

#' Using mean instead
SAD_mean <- as_tibble(data.frame(do.call(rbind,list_abund))) %>% #unlist sim
 filter(gen==ngenerations) %>% #  keeps final gen 
  group_by(sp, gen) %>% 
  summarise(n = mean(n)) %>% 
  pull(n)

SAD_octave <- octav(SAD_mean)

SAD_mean_plot <- plot(SAD_octave)

#' RAD instead
#' 

(mean_rad <- rad(SAD_mean))
plot(mean_rad)

#' Plotting manually
#' 
#+ SAD_mean
SAD_fig <- ggplot(data = SAD_mean, aes(x = n)) +
  geom_histogram(binwidth = 0.147)+ #logged
  # geom_histogram(binwidth =4.2)+
  geom_density(colour="red")+
  theme_classic()+
  ggtitle('SAD final')+
  ylab("Number of species")+
  xlab("Abundance (log)")+
  theme(panel.background = element_rect(fill='black'), 
        plot.background = element_rect(fill='black', color=NA),
        panel.grid.major = element_line(colour=NA),
        panel.grid.minor = element_line(color =NA),
        axis.line = element_line(colour = 'white', size = 0.8),
        axis.text = element_text(colour = 'white', family = "Arial", size = 12),
        axis.ticks = element_line(colour = "white"),
        axis.title = element_text(colour = "white", size = 16, family = "Arial"),
        plot.title = element_text(colour = "white")
  )+
  scale_y_continuous( expand=c(0,0)) +
  scale_x_continuous( expand=c(0,0))

SAD_fig

#' as a ridgeline plot

library(ggridges)
gen_filter <- c(1, seq(from = 10, to=ngenerations, by = 10))
  
SAD_mean_ridge <- as_tibble(data.frame(do.call(rbind,list_abund))) %>% #unlist sim
  group_by(sp, gen) %>% 
  summarise(n = mean(n)) %>% 
  mutate(gen=as.factor(gen)) %>% 
  filter(gen%in%gen_filter)

#+ SAD.ridge.plot
(SAD_ridge_plot <- ggplot(data = SAD_mean_ridge, aes(x=log(n), y=fct_rev(gen), fill = gen))+
  geom_density_ridges(color="white")+
  theme_ridges()+
  theme(legend.position = "none",
        panel.background = element_rect(fill='black'), 
        plot.background = element_rect(fill='black', color=NA),
        panel.grid.major = element_line(colour=NA),
        panel.grid.minor = element_line(color =NA),
        axis.line = element_line(colour = 'white', size = 0.8),
        axis.text = element_text(colour = 'white', family = "Arial", size = 12),
        axis.ticks = element_line(colour = "white"),
        axis.title = element_text(colour = "white", size = 16, family = "Arial"),
        plot.title = element_text(colour = "white"))+
  labs(x="log(Mean abundance)", y = "Time Step"))
ggsave( "test4.png", bg = 'transparent')
# Sp. Diversity ----
sp_div_hist <- as_tibble(data.frame(do.call(rbind,list_abund))) %>% #unlist sim
  mutate(sp = str_c("SP", sp)) %>% 
  group_by(rep, gen) %>% 
  tally() %>% 
  group_by(gen) %>% 
  summarise(mean_n=mean(n),
            SEM_n=sd(n)/sqrt(length(n)),
            max_n=max(n),
            min_n=min(n))

#+ sp.div.plot
div_fig <- ggplot(data=sp_div_hist, aes(x=gen,y=mean_n)) +
  geom_ribbon(aes(ymin=mean_n-SEM_n, ymax=mean_n+SEM_n), fill = "#35B779",  alpha = 0.6)+
  geom_line(colour = "#35B779") +
  theme_classic()+
  ggtitle('Species diversity')+
  ylab("Number of species")+
  xlab("Time step")+
  theme(panel.background = element_rect(fill='black'), 
        plot.background = element_rect(fill='black', color=NA),
        panel.grid.major = element_line(colour=NA),
        panel.grid.minor = element_line(color =NA),
        axis.line = element_line(colour = 'white', size = 0.8),
        axis.text = element_text(colour = 'white', family = "Arial", size = 12),
        axis.ticks = element_line(colour = "white"),
        axis.title = element_text(colour = "white", size = 16, family = "Arial"),
        plot.title = element_text(colour = "white"))+
  scale_y_continuous(limits = c(0, nspecies), expand=c(0,0)) +
  scale_x_continuous(limits = c(0, ngenerations), expand=c(0,0))

div_fig

# # Energy ----
e_data <- move_hist %>%                                       # calculating mean/se e-level
  group_by(repl, gen, sp) %>%
  summarise(e_avg = mean(e_val),
            e_sem = sd(e_val)/sqrt(length(e_val)))

e_trend <- move_hist %>%
  group_by(gen) %>%
  summarise(e_avg = mean(e_val),
            e_sem = sd(e_val)/sqrt(length(e_val)))

#+ energy.plot
e_fig <- ggplot(e_data, aes(x=gen, y=e_avg, colour=sp))+
  geom_line()+
  scale_colour_viridis(option="viridis", discrete = TRUE) +
  geom_line( data = e_trend,
             aes(x=gen,
                 y=e_avg),
             colour="black",
             inherit.aes = FALSE)+
  ggtitle('E-levels (mean)')+
  ylab("Energy")+
  xlab("Generation")+
  theme(panel.background = element_rect(fill='black'), 
        plot.background = element_rect(fill='black'),
        panel.grid.major = element_line(colour=NA),
        panel.grid.minor = element_line(color =NA),
        axis.line = element_line(colour = 'white', size = 0.8),
        axis.text = element_text(colour = 'white', family = "Arial", size = 12),
        axis.ticks = element_line(colour = "white"),
        axis.title = element_text(colour = "white", size = 16, family = "Arial"),
        plot.title = element_text(colour = "white"))+

  scale_y_continuous(limits=c(0,100), expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  theme(legend.position = "none")
e_fig
ggsave('test3.png')
# 
# # Panel ----
library(patchwork)

#+ panel.plot
history_4  <-  (sp_fig | div_fig)/( SAD_fig| e_fig)
history_panel <- history_4 / (SAD_ridge_plot)+
  plot_layout(heights = c(1,1,2))

ggsave('model_output_7mar.png', bg = 'transparent')


# history_panel  <-   ((sp_fig | SAD_fig) +
#                         plot_layout(c(2,1))
#                      /(div_fig | e_fig) +
#                        plot_layout(widths = c(1,1))
#                      /(SAD_ridge_plot)) +
#   plot_layout(heights =  c(2,1,3) ) +
#   plot_annotation(tag_levels = 'a')
# history_panel

# Animate ----
# movement_animation <- movement_fig +
#   transition_reveal(gen) +
#   labs(title = "Time step (gen) {as.integer(frame_along)} of 100")
# move_gif <- animate(movement_animation, end_pause = 3, width=800, height=800)
# anim_save(filename="move_anim.gif")

# sp_abund_animation <- sp_fig + 
#   transition_reveal(gen) 
# abund_gif <- animate(movement_animation, end_pause = 3, width=400, height=400)
# # anim_save(filename="move_anim.gif")
# 
# # stitch animations together
# # read the first image (frame) of each animation
# a <- image_read(p_dist_gif[[1]])
# b <- image_read(p_ROC_gif[[1]])
# c <- image_read(p_PR_gif[[1]])
# # combine the two images into a single image
# combined <- image_append(c(a, b, c))
# new_gif <- c(combined)
# for (i in 2:100) { # combine images frame by frame
#   a <- image_read(p_dist_gif[[i]])
#   b <- image_read(p_ROC_gif[[i]])
#   c <- image_read(p_PR_gif[[i]])
#   combined <- image_append(c(a, b, c))
#   new_gif <- c(new_gif, combined)
# }
# 
# # make an animation of the combined images
# combined_gif <- image_animate(new_gif)
# # save as gif
# image_write(combined_gif, "../animations/imbalance3.gif")





