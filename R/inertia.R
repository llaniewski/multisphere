
multisphere_inertia = function(tab) {
  scale = max(tab$r)
  tab[,c("x","y","z","r")] = tab[,c("x","y","z","r")]/scale
  ret = multisphere::sphere_integrals(tab)
  tab$x = tab$x - ret$center[1]
  tab$y = tab$y - ret$center[2]
  tab$z = tab$z - ret$center[3]
  se = eigen(ret$second)
  tab[,c("x","y","z")] = as.matrix(tab[,c("x","y","z")]) %*% se$vectors
  tab[,c("x","y","z","r")] = tab[,c("x","y","z","r")]*scale
  volume = ret$volume*scale^3
  second = se$values*scale^5
  inertia = as.vector((1-diag(3)) %*% second)
  list(tab=tab, volume=volume, center=rep(0,3), second=diag(second), inertia=diag(inertia))
}
