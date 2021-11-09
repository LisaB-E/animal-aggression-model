“Angry fish - linking changes in aggression after habitat loss to
demographic changes at the population level using IBM”
================
“LBE”
09 November, 2021

## Current structure

For loop 1 - repeat for each replicate  
- 2. Move  
- 3. Fight - dependent o energy  
- 4. Feed -  
- 5. Reproduce  
- 6. Energy gain/loss  
End FL1

# 1. Workspace preparation

``` r
rm(list=ls()) 
```

Load packages

``` r
packages = c('tidyverse', 'reshape','gstat','ggplot2',
             'viridis','tictoc', 'progress',
             'dplyr', 'gganimate')
load.pack = lapply(packages,require,char=T)
load.pack
```

# 2. Parameters

IBM Parameters *should include constraints on each parametre*. Double
check how many are needed in reduced model

``` r
ngenerations  = 10    # No. generations
replicates    = 10    # No. replicates (first half trans. second half intrans.)
dim           = 100   # dimension of square habitat array
hab_dim       = dim^2 # total no. cells
nspecies      = 100    # No. species
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
```

Potential one time-step movement. Can move to any adjacent cell (inc
diagonal)

``` r
step_moves = c(0,-1,dim-1,dim,dim+1,1,-dim+1,-dim,-dim-1)      #added option to stay in the same place
```

Potential aggression matrix - M1 assume all species have equal
competitive potential - so outcome will only depend on energy levels

``` r
aggression = matrix(0.5, ncol=nspecies, nrow = nspecies)
diag(aggression) = 0.5 # For interaction between individuals of the same species 
```

Store species richness at each generation

``` r
rich = data.frame(gen=rep(1:ngenerations,replicates),
                  rep=rep(1:replicates,each=ngenerations),
                  rich=NA,
                  type=NA)
```

Progress bar *use pb$tick() in for-loop to run*

``` r
pb <- txtProgressBar(min = 0, max = replicates*ngenerations, style = 3)
```

    ##   |                                                                              |                                                                      |   0%

unique ID’s for each individual

``` r
IDs <- c(1:(nspecies*nindiv*ngenerations*replicates)) #does it matter if this is unique per replicate, or unique throughout? right now unique throughout
```

# 3. Initialise - habitat

Simulate habitat

``` r
  hab_grid = expand.grid(1:dim, 1:dim)                                     # The habitat arena
  names(hab_grid) = c('x','y')
  hab_bin = sample(x=c(1,0), size = hab_dim, replace = T)                  # samples binary habitat values
  hab_bin = cbind(hab_grid, hab_bin)
```

Make matrix

``` r
  hab_vals = as.matrix(cast(hab_bin,
                            x~y,
                            value='hab_bin'))
```

## Plot habitat values

``` r
  hab.plot = ggplot(hab_bin) + theme_bw() +
    geom_raster(aes(x,y,fill=hab_bin)) +
    scale_x_continuous(expand=expand_scale(add=0)) +
    scale_y_continuous(expand=expand_scale(add=0)) +
    theme(legend.position = 'top',
          axis.text=element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
```

    ## Warning: `expand_scale()` is deprecated; use `expansion()` instead.

    ## Warning: `expand_scale()` is deprecated; use `expansion()` instead.

``` r
 hab.plot
```

![](general-model_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
#FIGURE OUT WHIHC LOCS TO KEEP ORIGINAL VS UPDATED...
```

Randomly place individuals of each species in habitat

``` r
locs = cbind(sample(hab_dim,tot_indiv),rep(1:nspecies,each=nindiv))       # Randomly picks locations (V1-cell) for all individuals (V2) 
locs = cbind(locs,hab_vals[locs[,1]])                                     # habitat value
```

add energetic value

``` r
locs = cbind(locs,rep(100,nindiv))                                        # starting energetic value =100?
locs = locs[order(locs[,1]),] 

locs <- cbind(locs, IDs[1:nrow(locs)])
IDs <- setdiff(IDs,locs[,5])
colnames(locs)=c("loc", "sp", "hab_val", "e_val", "ID")
OG_locs <- locs
```

# 4. Functions

### Function - move

``` r
move <- function(){      
  
  # Place individuals in habitat
  habitat = matrix(0,ncol=dim,nrow=dim)                                # make 100 x 100 matrix
  habitat[locs[,1]] = 1                                                # places individual in the matrix
  
# Find edge and corner individuals
  edge_corner <- unique(c(1:dim,                                       # top row
                          which(1:hab_dim%%dim%in%0:1),                # sides
                          (hab_dim-dim+1):hab_dim))                    # bottom
  
  
ec_occ      <- edge_corner[which(habitat[edge_corner]==1)]          # which edge corner cells are occupied
ec_indiv    <- locs[locs[,1]%in%ec_occ, , drop=FALSE]              # which individuals are on edge corner cells (V2&V1)
  offec_indiv <- locs[locs[,1]%in%setdiff(locs[,1],ec_indiv[,1]),]     # individuals  NOT on edge corner cells (V1&v2)
  
  offec_indiv <- cbind(offec_indiv, sample(step_moves, dim(offec_indiv)[1], replace = TRUE)) # add move step to non edge inds (V3)
  
# Move edge and corner individuals
#goal: figure out how to do this neater (not sure the if else statements are neede for samppling of the step_moves, surely a dim=1 woudl replace the if TRUE part?)*
  
  ec_indiv_new <- NULL
  if(sum(ec_indiv[,1]==1)>0){                                         # Q1 find if any inds are in top left corner
    top_left_ind = ec_indiv[ec_indiv[,1]==1,]                         # VALUE IF TRUE -pull out that data
    if(is.null(dim(top_left_ind))){                                     # Q2 if top left individual does not exist
      top_left_ind = c(top_left_ind,sample(step_moves[c(1,4:6)],1))          # VALUE IF TRUE move down, right or diagonal
    }else{                                     
      top_left_ind = cbind(top_left_ind,sample(step_moves[c(1,4:6)],         # VALUE IF FALSE
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
  if(sum(ec_indiv[,1] == dim)>0){                                                     #bottom left
    bottom_left_ind = ec_indiv[ec_indiv[,1]==dim,]
    if(is.null(dim(bottom_left_ind))){
      bottom_left_ind = c(bottom_left_ind,sample(step_moves[c(1, 2:4)],1))
    }else{
      bottom_left_ind = cbind(bottom_left_ind,sample(step_moves[c(1,2:4)],
                                                     dim(bottom_left_ind)[1],r=T))
    }
    ec_indiv_new = rbind(ec_indiv_new,bottom_left_ind)
  }
  if(sum(ec_indiv[,1] == hab_dim)>0){                                                  #bottom righ
    bottom_right_ind = ec_indiv[ec_indiv[,1]==hab_dim,]
    if(is.null(dim(bottom_right_ind))){
      bottom_right_ind = c(bottom_right_ind,sample(step_moves[c(1:2,8:9)],1))
    }else{
      bottom_right_ind = cbind(bottom_right_ind,sample(step_moves[c(1:2,8:9)],
                                                       dim(bottom_right_ind)[1],r=T))
    }
    ec_indiv_new = rbind(ec_indiv_new,bottom_right_ind)
  }
  if(sum(ec_indiv[,1]%in%(2:(dim-1)))>0){                                              #left edge
    left_edge_ind = ec_indiv[ec_indiv[,1]%in%(2:(dim-1)),]
    if(is.null(dim(left_edge_ind))){
      left_edge_ind = c(left_edge_ind,sample(step_moves[c(1,2:6)],1))
    }else{
      left_edge_ind = cbind(left_edge_ind,sample(step_moves[c(1,2:6)],
                                                 dim(left_edge_ind)[1],r=T))
    }
    ec_indiv_new = rbind(ec_indiv_new,left_edge_ind)
  }
  if(sum(ec_indiv[,1]%in%((hab_dim-dim+2):(hab_dim-1)))>0){                             #right edge
    right_edge_ind = ec_indiv[ec_indiv[,1]%in%((hab_dim-dim+2):(hab_dim-1)),]
    if(is.null(dim(right_edge_ind))){
      right_edge_ind = c(right_edge_ind,sample(step_moves[c(1:2,6:9)],1))
    }else{
      right_edge_ind = cbind(right_edge_ind,sample(step_moves[c(1:2,6:9)],
                                                   dim(right_edge_ind)[1],r=T))
    }
    ec_indiv_new = rbind(ec_indiv_new,right_edge_ind)
  }
  if(sum(ec_indiv[,1]%in%setdiff(which(1:hab_dim%%dim==1),c(1,hab_dim-dim+1)))>0){        #top edge
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
  if(sum(ec_indiv[,1]%in%setdiff(which(1:hab_dim%%dim==0),c(dim,hab_dim)))>0){                #bottom edge
    bottom_edge_ind = ec_indiv[ec_indiv[,1]%in%setdiff(which(1:hab_dim%%dim==0),
                                                       c(dim,hab_dim)),]
    if(is.null(dim(bottom_edge_ind))){
      bottom_edge_ind = c(bottom_edge_ind,sample(step_moves[c(1:4,8:9)],1))
    }else{
      bottom_edge_ind = cbind(bottom_edge_ind,sample(step_moves[c(1:4,8:9)],
                                                     dim(bottom_edge_ind)[1],r=T))
    }
    ec_indiv_new = rbind(ec_indiv_new,bottom_edge_ind)
    
  }                                                                       #end find edge individs section
  
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
locs_new <- cbind(locs_new, gen, repl)                               # add generation number V6, repl V7

# Getting x & y coords for animation
temp_habitat1 <- matrix(0,ncol=dim,nrow=dim)                          # make 100 x 100 matrix
temp_habitat1[locs_new[,1]] = 1                                       # place inds on habitat, ignores later dupes
temp_coords1 <- which(temp_habitat1>0, arr.ind=TRUE)                  # get x & y coords

triple_threat_list <- list()                                          # get duplicates, including triples etc
triple_threat_list <- NULL
temp_dupes <- matrix(locs_new[,1][duplicated(locs_new[,1], fromLast = TRUE)])
  
for (dupe in 1:nrow(temp_dupes)) {
  temp_habitat <- matrix(0,ncol=dim,nrow=dim)
  temp_habitat[temp_dupes[dupe,1 ]] = 1
  temp_coordsX <- which(temp_habitat>0, arr.ind=TRUE)
  triple_threat_list[[dupe]] <- temp_coordsX
}

temp_coords2 <- do.call(rbind,triple_threat_list)                      # extract each list
temp_coords <- rbind(temp_coords1, temp_coords2)                       # merge coords inclduing dupes
temp_coords <- temp_coords[order(temp_coords[,2], temp_coords[,1]),]     # reorder
locs_new <- cbind(locs_new, temp_coords)

return(locs_new)                                                        
  
}
```

### Function - fight

``` r
fight <- function(){
  locs_new <- locs     #pulls from move endpoint in sim                 
  # Interacting species are those occupying the same cell 
  loc_ind <- locs_new[,1]                                                 # extracts cell locations
  int_loc <- loc_ind[duplicated(loc_ind)]                                 # extracts duplicated values from new locations (ie where individuals meet)
  fighters <- locs_new[locs_new[,1]%in%int_loc,]                          # extracts values for fighters
  locs_new = locs_new[locs_new[,1]%in%setdiff(locs_new[,1],fighters[,1]),]             # Remove fighters from set of individuals            #BEWARE!!! if dplyr is loaded - changes results! hould chekc this forst! NEED to fix
  
  fighters[,4] = fighters[,4]-fight_eloss                              # Aggression is energetically expensive
  # *DOES THIS MAKE SENSE? IT'S ASSUMING THEY MAKE THIS DECISION AHEAD OF TIME, I WILL OR WILL NOT ENGAGE IN  AGGRESISON BASED ON ENERGY LEVELS, OR SHOULD THIS BE AFTER??*
 
  # Roulette selection for winners
  all_winners = NULL
 
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
                                                          
      while(rowsum!=5){                                                 # 5 = all cells matches with the loser
        rowcount = rowcount+1
        rowsum = sum(sp_in_cell[rowcount,] == loser)                    # checks whether top row = loser
      }
       sp_in_cell = sp_in_cell[-rowcount,,drop=F]                       # removes loser from sp in cell
      }      
      }                                                                 # loops back to next fightJ
    
    initial = fighters[fighters[,1]==i,]                                # identifies starting individs in cell (.ie OG sp_in_cell)
   
# identify which individ is winner to add to win_vect tally
   if(nrow(sp_in_cell)>0) {                                              # if all were low energy
     rowcount = 0
    rowsum = 0
    while(rowsum!=5){                                 
      rowcount = rowcount + 1
      rowsum = sum(initial[rowcount,] == sp_in_cell)                    # finds winner in sp_in_cell
    }
    win_vect[rowcount] = 1                                               # changes 0 to 1 in win vs lose column
   }
    all_winners = c(all_winners,win_vect)
    }                                                                   #loops back to next cell with multiples
  fighters = cbind(fighters,all_winners)                               # merges 
  
 # Lose energy
  fighters[,4] = fighters[,4]-fight_eloss                              # Aggression is energetically expensive DOES THIS MAKE SENSE? IT'S ASSUMING THEY MAKE THIS DECISION AHEAD OF TIME, I WILL OR WILL NOT ENGAGE  IN AGGRESISON BASED ON ENERGY LEVELS, OR SHOULD THIS BE AFTER?? also, why arwe there two of these, once before and once after?
  winners = fighters[fighters[,6]==1, , drop=FALSE]                                   # Differentiate winners and losers 
  losers = fighters[fighters[,6]==0, ,drop=FALSE]
  
  losers = losers[losers[,4]>0, ,drop=FALSE]                                       # deletes dead individuals (energy =0). I wonder if there should alos be random mortality in these guys?

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
  offspring[,4] = offspring[,4]/offspring_pen                            # Set offspring energy to 50% of parents
  
  locs_new = rbind(locs_new[,1:5],offspring[,1:5],losers[,1:5])          # Combines all after fight, feed, death, reproduction
  locs_new = locs_new[order(locs_new[,1]),]


  

  gen <- rep(generation, nrow(locs_new))
  # Store gen (WARNING - change 1 to generation, replic)
  repl <- rep(replic, nrow(locs_new))
  locs_new <- cbind(locs_new, gen,repl)  
  
    return(locs_new)
  
  }
```

# 5. Simulate

``` r
loop=0
list_move=list()                                                      
list_fight=list()
list_abund=list()
list_rich = list()
move_hist = tibble()

## for trouible shooting - delet when done

######
#replication loop
for (replic in 1:replicates) { 
locs <- OG_locs #reset to starting point - is this right level of replication?

#generation loop
for(generation in 1:ngenerations){
  gen_rep <- paste(generation,replic, sep = "_")
  loop = loop + 1
  
# Move 
  after_move <- move()                        # after_move = new_locs
  locs = after_move[,1:5]                     # updates locs to endpoint (final timestep) 
  list_move[[gen_rep]] <- after_move       # stores moves per generation as list
# print(paste('move_',' gen', generation, ' rep', replic, ' DONE', sep = ""))
  
# Fight
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

#save move data
move_temp <- data.frame(do.call(rbind,list_move)) %>% #unnest nested list
  mutate(sp = str_c("SP", sp, sep = " "),
         ID = as_factor(ID),
         U_ID = str_c(sp, ID, sep="_"))
move_hist <- bind_rows(move_temp, move_hist)          #adds each gen rep to each other #IS THIS THE PROBLEM HERE - ARE WE MERGING WEIRDLY?

#save fight data 


}                                             # end rep
```

    ##   |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  23%  |                                                                              |=================                                                     |  24%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  63%  |                                                                              |=============================================                         |  64%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  77%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |===========================================================           |  84%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================| 100%

``` r
close(pb)
```

``` r
#to do next, figure out why species are disapreading between gens - locs_new has 1-100, but output is reduced to 7 species, what's goin on there? The dataset is somehow cropped to only incldue whatever species make it through to gen100... must be an issue with the reps vs generations? figure out that looping to figure out what is goin on...
```

``` r
#can get rid of U_ID?

  # Movement ----
random_select = sample(x=move_hist$U_ID,size=100)                               # reduced set for plotting n individs
sub_movers <- subset(move_hist, U_ID %in% c(random_select))

movement_fig <- ggplot(data=hab_bin, aes(x=x,y=y)) +
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
  coord_fixed()

movement_fig
```

![](general-model_files/figure-gfm/-%206.%20VISUALISE-1.png)<!-- -->

``` r
# Sp. Abundance ----
sp_hist <- as_tibble(data.frame(do.call(rbind,list_abund))) %>% #unlist sim
  mutate(sp = str_c("SP", sp)) %>% 
  group_by(gen, sp) %>% 
  summarise(mean_abund = mean(n),
            SEM_abund = sd(n)/sqrt(length(n)))
```

    ## `summarise()` has grouped output by 'gen'. You can override using the `.groups` argument.

``` r
abund_running <- data.frame()                                                   # generate overall mean and se (across spp)
for (genx in 1:ngenerations) {
  sub_abund_data <- subset(sp_hist, gen==genx)
  sub_abund <- data.frame(gen = genx,
                          meanX = mean(sub_abund_data$n),
                          sem = sd(sub_abund_data$n)/sqrt(length(sub_abund_data$n)))
  abund_running <- rbind(abund_running, sub_abund)

}
```

    ## Warning: Unknown or uninitialised column: `n`.

    ## Warning in mean.default(sub_abund_data$n): argument is not numeric or logical:
    ## returning NA

    ## Warning: Unknown or uninitialised column: `n`.

    ## Warning: Unknown or uninitialised column: `n`.

    ## Warning: Unknown or uninitialised column: `n`.

    ## Warning in mean.default(sub_abund_data$n): argument is not numeric or logical:
    ## returning NA

    ## Warning: Unknown or uninitialised column: `n`.

    ## Warning: Unknown or uninitialised column: `n`.

    ## Warning: Unknown or uninitialised column: `n`.

    ## Warning in mean.default(sub_abund_data$n): argument is not numeric or logical:
    ## returning NA

    ## Warning: Unknown or uninitialised column: `n`.

    ## Warning: Unknown or uninitialised column: `n`.

    ## Warning: Unknown or uninitialised column: `n`.

    ## Warning in mean.default(sub_abund_data$n): argument is not numeric or logical:
    ## returning NA

    ## Warning: Unknown or uninitialised column: `n`.

    ## Warning: Unknown or uninitialised column: `n`.

    ## Warning: Unknown or uninitialised column: `n`.

    ## Warning in mean.default(sub_abund_data$n): argument is not numeric or logical:
    ## returning NA

    ## Warning: Unknown or uninitialised column: `n`.

    ## Warning: Unknown or uninitialised column: `n`.

    ## Warning: Unknown or uninitialised column: `n`.

    ## Warning in mean.default(sub_abund_data$n): argument is not numeric or logical:
    ## returning NA

    ## Warning: Unknown or uninitialised column: `n`.

    ## Warning: Unknown or uninitialised column: `n`.

    ## Warning: Unknown or uninitialised column: `n`.

    ## Warning in mean.default(sub_abund_data$n): argument is not numeric or logical:
    ## returning NA

    ## Warning: Unknown or uninitialised column: `n`.

    ## Warning: Unknown or uninitialised column: `n`.

    ## Warning: Unknown or uninitialised column: `n`.

    ## Warning in mean.default(sub_abund_data$n): argument is not numeric or logical:
    ## returning NA

    ## Warning: Unknown or uninitialised column: `n`.

    ## Warning: Unknown or uninitialised column: `n`.

    ## Warning: Unknown or uninitialised column: `n`.

    ## Warning in mean.default(sub_abund_data$n): argument is not numeric or logical:
    ## returning NA

    ## Warning: Unknown or uninitialised column: `n`.

    ## Warning: Unknown or uninitialised column: `n`.

    ## Warning: Unknown or uninitialised column: `n`.

    ## Warning in mean.default(sub_abund_data$n): argument is not numeric or logical:
    ## returning NA

    ## Warning: Unknown or uninitialised column: `n`.

    ## Warning: Unknown or uninitialised column: `n`.

``` r
sp_fig <- ggplot( data = sp_hist, aes(x=gen, y=mean_abund, color=sp)) +
  #pverall mean abundance
  geom_ribbon(data = abund_running,
              aes(x=gen,
                  ymin=meanX-sem,
                  ymax=meanX+sem),
              inherit.aes = FALSE)+
  geom_line(data = abund_running,
            aes(x=gen,
                y=meanX),
            colour="black")+
  #species specific abundance
  geom_ribbon(aes(ymin = mean_abund-SEM_abund,
                  ymax = mean_abund+SEM_abund),
              alpha =.5,
              linetype = 0)+
  geom_line() +
  ggtitle('Species abundance')+
  ylab("Abundance (by species)")+
  xlab("Generation")+
  scale_colour_viridis(option="viridis", discrete = TRUE) +
  theme_dark()+
  theme(title = element_text(colour = "white"),
        plot.background = element_rect(fill = "grey10"),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour=NA),
        panel.grid.minor = element_line(color =NA),
        axis.title = element_text(colour="white", family = "Arial"),
        axis.line = element_line(colour = "grey30"),
        axis.text = element_text(family = "Arial")
  ) +
  theme(legend.position = "none")

sp_fig
```

    ## Warning in max(ids, na.rm = TRUE): no non-missing arguments to max; returning
    ## -Inf

    ## Warning: Removed 10 row(s) containing missing values (geom_path).

![](general-model_files/figure-gfm/-%206.%20VISUALISE-2.png)<!-- -->

``` r
# SAD ----------
# (contd'd) from above
# idea: animate this thorugh time using gganimate()
# SAD_data <- data.frame(subset(sp_hist, gen==ngenerations))
# SAD_missing <- data.frame(sp=setdiff(sp_hist[,1], SAD_data[,1]), # finds extinct species
#                           n=0,
#                           gen=ngenerations)
# SAD_data <- rbind(SAD_data, SAD_missing)
# 
# SAD_fig <- ggplot(data = SAD_data, aes(x = n)) +
#   geom_histogram(aes(y=..density..),binwidth = 5)+
#   geom_density(colour="red")+
#   theme_dark()+
#   ggtitle('Species abundance distribution')+
#   ylab("Number of species")+
#   xlab("Abundance")+
#   theme(title = element_text(colour = "white"),
#         plot.background = element_rect(fill = "grey10"),
#         panel.background = element_blank(),
#         panel.grid.major = element_line(colour=NA),
#         panel.grid.minor = element_line(color =NA),
#         axis.title = element_text(colour="white", family = "Arial"),
#         axis.line = element_line(colour = "grey30"),
#         axis.text = element_text( family = "Arial")
#   )+
#   scale_y_continuous( expand=c(0,0)) +
#   scale_x_continuous( expand=c(0,0))
# SAD_fig

# # Sp. Diversity ----
# sp_div <-  data.frame()                                          # Calculating no.sp in each generation
# for (genX in 1:ngenerations) {
#   gen_div <- data.frame(gen=genX, 
#   div=length(unique(subset(sp_hist, gen==genX))$sp))
#   sp_div=rbind(sp_div,gen_div)
#   }
# 
#   
# div_fig <- ggplot(sp_div, aes(x=gen,y=div)) +
#   geom_line(colour = "#35B779") +
#   theme_dark()+
#   ggtitle('Species diversity')+
#   ylab("Number of species")+
#   xlab("Generation")+
#   theme(title = element_text(colour = "white"),
#         plot.background = element_rect(fill = "grey10"),
#         panel.background = element_blank(),
#         panel.grid.major = element_line(colour=NA),
#         panel.grid.minor = element_line(color =NA),
#         axis.title = element_text(colour="white", family = "Arial"),
#         axis.line = element_line(colour = "grey30"),
#         axis.text = element_text( family = "Arial")
#         )+
#   scale_y_continuous(limits = c(0, nspecies), expand=c(0,0)) +
#   scale_x_continuous(limits = c(0, ngenerations), expand=c(0,0))
#   
# div_fig
# 
# # Energy ----
# e_data <- move_hist %>%                                       # calculating mean/se e-level
#   group_by(gen, sp) %>%
#   summarise(e_avg = mean(e_val),
#             e_sem = sd(e_val)/sqrt(length(e_val)))
# 
# e_trend <- move_hist %>%
#   group_by(gen) %>%
#   summarise(e_avg = mean(e_val),
#             e_sem = sd(e_val)/sqrt(length(e_val)))
# 
# e_fig <- ggplot(e_data, aes(x=gen, y=e_avg, colour=sp))+
#   geom_line()+
#   geom_ribbon(data= e_trend, 
#               aes(x=gen, 
#                   ymin=e_avg-e_sem, 
#                   ymax=e_avg+e_sem), 
#               alpha=0.7, 
#               inherit.aes = FALSE)+
#   geom_line( data = e_trend,
#              aes(x=gen,
#                  y=e_avg),
#              colour="black",
#              inherit.aes = FALSE)+
#   theme_dark()+
#   ggtitle('Energy levels (mean per species)')+
#   ylab("Energy")+
#   xlab("Generation")+
#   theme(title = element_text(colour = "white"),
#         plot.background = element_rect(fill = "grey10"),
#         panel.background = element_blank(),
#         panel.grid.major = element_line(colour=NA),
#         panel.grid.minor = element_line(color =NA),
#         axis.title = element_text(colour="white", family = "Arial"),
#         axis.line = element_line(colour = "grey30"),
#         axis.text = element_text( family = "Arial")
#   )+
#   scale_y_continuous(limits=c(0,100), expand=c(0,0)) +
#   scale_x_continuous(expand=c(0,0)) +
#   theme(legend.position = "none") 
# e_fig
# 
# # Panel ----
# library(patchwork)
# 
# history_panel  <-   ((sp_fig | SAD_fig +
#                         plot_layout(c(2,1))) 
#                      /(div_fig | e_fig) +
#                        plot_layout(widths = c(1,1))) +
#   plot_layout(heights =  c(2,1) ) +
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
```
