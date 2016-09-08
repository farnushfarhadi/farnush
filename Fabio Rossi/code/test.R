#test 1

v1 <- c (10 , 12 , 4 , 20 , 16 , 18)
v2 <- c(9 , 13 , 6 , 2 , -2 , 1)
cor (v1 , v2)
cor (v1[1:3] , v2[1:3])
cor (v1[4:6] , v2[4:6])
stdVector <- function(x) { return ((x - mean(x))/sd(x) )}

stdVector(v1[1:3] ) -> v1_h1_std
stdVector(v1[4:6] ) -> v1_h2_std
stdVector(v2[1:3] ) -> v2_h1_std
stdVector(v2[4:6] ) -> v2_h2_std

par (mfrow = c(2, 2))
plot ( v1 , x = c (1:6) , col = "blue" , type = "b")
plot ( c(v1_h1_std , v1_h2_std) , x = c (1:6) , col = "blue" , type = "b")
plot ( v2 , x = c (1:6) , col = "blue" , type = "b")
plot ( c(v2_h1_std , v2_h2_std) , x = c (1:6) , col = "blue" , type = "b")

v1_com <-  c(v1_h1_std , v1_h2_std)
v2_com <- c(v2_h1_std , v2_h2_std)
cor (v1_com , v2_com)

#test 2

#v1 <- c (8 , 12 , 10 , 2 , 14 , 11)
v1 <- c (8 , 12 , 10 , 4 , 1 , 6)
#v2 <- c(7 , 13 , 8 , 10 , 16 , 14)
v2 <- c(8 , 12 , 10 , 16 , 9 , 18)
cor (v1 , v2)
cor (v1[1:3] , v2[1:3])
cor (v1[4:6] , v2[4:6])

stdVector(v1[1:3] ) -> v1_h1_std
stdVector(v1[4:6] ) -> v1_h2_std
stdVector(v2[1:3] ) -> v2_h1_std
stdVector(v2[4:6] ) -> v2_h2_std

par (mfrow = c(2, 2))
plot ( v1 , x = c (1:6) , col = "blue" , type = "b")
plot ( c(v1_h1_std , v1_h2_std) , x = c (1:6) , col = "blue" , type = "b")
plot ( v2 , x = c (1:6) , col = "blue" , type = "b")
plot ( c(v2_h1_std , v2_h2_std) , x = c (1:6) , col = "blue" , type = "b")

v1_com <-  c(v1_h1_std , v1_h2_std)
v2_com <- c(v2_h1_std , v2_h2_std)
cor (v1_com , v2_com)





