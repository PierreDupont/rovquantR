# 
# lowerCoords <- habitat$scaledLowerCoords
# upperCoords <- habitat$scaledUpperCoords

makeGridFromCorners <- function(lowerCoords,
                                upperCoords,
                                crs = NULL){
  if(!all(dim(lowerCoords) == dim(upperCoords))){
    stop("the number of upper and lower coordinates do not match!")
  }
  
  n.cells <- nrow(lowerCoords)
  ret = vector("list", n.cells)
  square = function(lower, upper){
    upper <- as.numeric(upper)
    lower <- as.numeric(lower)
    corners = matrix(c(lower,
                       upper[1],lower[2],
                       upper,
                       lower[1],upper[2],
                       lower),
                     ncol = 2, byrow = TRUE)
    pts = list(corners)
    st_polygon(pts)
  }
  


  for (c in 1:n.cells){
  ret[[c]] = square(lower = lowerCoords[c,], 
                    upper = upperCoords[c,])
  } #c
  sf::st_sfc(ret, crs = crs)
}

# 
# test <- makeGridFromCorners(st_coordinates(habitat$lower.hab.sp),
#                             st_coordinates(habitat$upper.hab.sp),
#                             st_crs(habitat$habitat.sp))
# 
# plot(test)
# plot(test[habitat$buffered.habitat.poly])
# 
