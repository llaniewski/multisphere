
multisphere_inertia = function(tab, centering=TRUE, align=FALSE) {
  scale = max(tab$r)
  tab[,c("x","y","z","r")] = tab[,c("x","y","z","r")]/scale
  ret = multisphere::sphere_integrals(tab)
  center = ret$center
  if (centering) {
    tab$x = tab$x - center[1]
    tab$y = tab$y - center[2]
    tab$z = tab$z - center[3]
    center[] = 0
  }
  second = ret$second
  if (align) {
    se = eigen(second)
    tab[,c("x","y","z")] = as.matrix(tab[,c("x","y","z")]) %*% se$vectors
    second = diag(x=se$values) # second = t(se$vectors) %*% second %*% se$vectors
  }
  tab[,c("x","y","z","r")] = tab[,c("x","y","z","r")]*scale
  volume = ret$volume*scale^3
  second = second*scale^5
  inertia = sum(diag(second))*diag(nrow = 3) - second
  list(tab=tab, volume=volume, center=center, centering=centering, second=second, inertia=inertia, align=align)
}
