stack_subexps <- function(data, timeID, groupID, adoptionTime, kappa_pre, kappa_post){
  # create the sub-experimental data sets
  events = data[is.na(adopt_year) == FALSE, funique(adopt_year)]
  
  # make a list to store the sub experiments in.
  sub_experiments = list()
  
  # Loop over the events and make a data set for each one
  for (j in events) {
    sub_name = paste0("sub_", j) 
    sub_experiments[[sub_name]] = create_sub_exp(
      dataset = data,
      timeID = timeID,
      groupID = groupID, 
      adoptionTime = adoptionTime, 
      focalAdoptionTime = j,
      kappa_pre = kappa_pre,
      kappa_post = kappa_post)
  }
  
  # Vertically concatenate the sub-experiments
  stackfull = rbindlist(sub_experiments)
  
  # Remove the sub-experiments that are not feasible
  stacked_data = stackfull[feasible == 1]
  
  # Treated, control, and total count by stack
  stacked_data[event_time==0, 
               .(N_treated = fsum(treat), 
                 N_control = fsum(1-treat), 
                 N_total = .N
               ), 
               by = sub_exp][order(sub_exp)]  
  
  return(stacked_data)
}


create_sub_exp = function(dataset, timeID, groupID, adoptionTime, focalAdoptionTime, kappa_pre, kappa_post){
  
  # Copy dataset 
  dt_temp = copy(dataset)
  
  # Determine earliest and latest time in the data. 
  # Used for feasibility check later
  minTime = dt_temp[, fmin(get(timeID))]
  maxTime = dt_temp[, fmax(get(timeID))]
  
  # Include only the treated groups and the clean controls that adopt at least kappa_post periods after the focal atime.
  dt_temp = dt_temp[get(adoptionTime) == focalAdoptionTime | get(adoptionTime) > focalAdoptionTime + kappa_post | get(adoptionTime) == TRUE | is.na(get(adoptionTime))]
  
  # Limit to time periods inside the event window defined by the kappas
  dt_temp = dt_temp[get(timeID) %in% (focalAdoptionTime - kappa_pre):(focalAdoptionTime + kappa_post)]
  
  # Make treatment group dummy
  dt_temp[, treat := 0]
  dt_temp[get(adoptionTime) == focalAdoptionTime, treat := 1] 
  
  # Make a post variable
  dt_temp[, post := fifelse(get(timeID) >= focalAdoptionTime, 1, 0)]
  
  # Make event time variable
  dt_temp[, event_time := get(timeID) - focalAdoptionTime]
  
  # Create a feasible variable
  dt_temp[, feasible := fifelse(focalAdoptionTime - kappa_pre >= minTime & focalAdoptionTime + kappa_post <= maxTime, 1, 0)]
  
  # Make a sub experiment ID.
  dt_temp[, sub_exp := focalAdoptionTime]
  
  return(dt_temp)
} 



compute_weights = function(dataset) {
  
  # Create a copy of the underlying dataset
  stack_dt_temp = copy(dataset)
  
  # Step 1: Compute stack - time counts for treated and control
  stack_dt_temp[, `:=` (stack_n = .N,
                        stack_treat_n = sum(treat),
                        stack_control_n = sum(1 - treat)), 
                by = event_time
  ]  
  # Step 2: Compute sub_exp-level counts
  stack_dt_temp[, `:=` (sub_n = .N,
                        sub_treat_n = sum(treat),
                        sub_control_n = sum(1 - treat)
  ), 
  by = list(sub_exp, event_time)
  ]
  
  # Step 3: Compute sub-experiment share of totals
  stack_dt_temp[, sub_share := sub_n / stack_n]
  
  stack_dt_temp[, `:=` (sub_treat_share = sub_treat_n / stack_treat_n,
                        sub_control_share = sub_control_n / stack_control_n
  )
  ]
  
  # Step 4: Compute weights for treated and control groups
  stack_dt_temp[treat == 1, stack_weight := 1]
  stack_dt_temp[treat == 0, stack_weight := sub_treat_share/sub_control_share]
  
  return(stack_dt_temp)
}  


create_stacked_sample = function(data, timeID, groupID, adoptionTime, kappa_pre, kappa_post){
  stacked_data = stack_subexps(data, timeID = timeID, groupID = groupID, adoptionTime = adoptionTime, 
                                kappa_pre = kappa_pre, kappa_post = kappa_post)
  stacked_data = compute_weights(stacked_data)
}


stacked_sample <- create_stacked_sample(
  dtc,
  timeID = "year",
  groupID = "statefips",
  adoptionTime = "adopt_year",
  kappa_pre = 3,
  kappa_post = 2)

stacked_sample


fepois(unins ~ i(event_time, treat, ref = -1) | statefip^sub_exp + event_time^sub_exp, 
      data = stacked_dtc2, 
      cluster = ~statefip,
      weights = ~stack_weight)


