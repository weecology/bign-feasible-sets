library(VGAM)

par_est_zipf = function(ab, N = round(sum(ab))){
 # Function to estimate the parameter for the zipf distribution
 fit = vglm(round(ab) ~ 1, zipf(link = identity, N = N))
 return(coef(fit))
}

ppf_zipf = function(cdf, N, s){
  cdf_eqn = function(x, cdf, N, s){
    return(pzipf(x, N, s) - cdf)
  }
  x = uniroot(cdf_eqn, c(1, N), cdf = cdf, N = N, s = s)
  return(round(x$root))
}

get_zipf_rad = function(ab){
  # Function to obtain the rad predicted by Zipf distribution
  s = par_est_zipf(ab)
  S = length(ab)
  N = round(sum(ab))
  ans = c()
  for (i in 1:S){
    cdf = (S-i+0.5)/S
    if (cdf <= pzipf(1, N, s)){
      ans = c(ans, 1)
    }
    else {
      ans = c(ans, ppf_zipf(cdf, N, s))
    }
  }
  return(ans)
}

get_pred_zipf_single = function(folder_name){
  #Given the directory of dataset, obtain predicted abundances
  #of the Zipf distribution and save to file
  dat_dir = paste('./', folder_name, '/', folder_name, '-data-cleaned.txt', sep = "")
  out_dir = paste('./', folder_name, '/', folder_name, '-obs-pred-zipf.txt', sep = "")
  dat = read.delim(dat_dir, colClasses = c('character', 'numeric'), header = F)
  names(dat) = c('site', 'obs')
  row_no_zero = apply(dat, 1, function(row) all(row != 0))
  dat = dat[row_no_zero, ] # Remove rows with zeros
  out_file = data.frame(site = character(0), obs = numeric(0), pred = numeric(0))
  site_list = unique(dat$site)
  for (i_site in 1:length(site_list)){
    dat_site = dat[dat$site == site_list[i_site], ]
    obs_site = dat_site$obs
    if (length(obs_site) >= 5){
      pred_site = get_zipf_rad(obs_site)
      out_site = cbind(dat$site[dat$site == site_list[i_site]],
        sort(obs_site, decreasing = T), pred_site)
      names(out_site) = c('site', 'obs', 'pred')
      out_file = rbind(out_file, out_site)
      write.table(out_file, file = out_dir, quote = F, row.names = F, col.names = F)
    }
  }
}
