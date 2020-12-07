
ext <- read.csv("extraction_variations.csv")
head(ext)
par(mfrow=c(2,1), mar=c(4,4,2,4))
plot(concentration ~ treatment, ext, pch = 19, cex = 2, col = sample, xlab="", xaxt="n", ylab = "concentration")
axis(1, at = 1:5, labels = c("10 ul\nBlood", "20ul 30Prok\nBlood", "10ul\nTissue",
                             "10ul\nBlood", "10ul\nTissue"))
mtext("-------------22h-----------                        ----7h-----", side=1, line = 3, cex=1.5)

plot(ratio ~ treatment, ext, pch = 19, cex = 2, col = sample, xaxt="n", xlab="", ylab="260/230")
axis(1, at = 1:5, labels = c("10 ul\nBlood", "20ul 30Prok\nBlood", "10ul\nTissue",
                             "10ul\nBlood", "10ul\nTissue"))
mtext("-------------22h-----------                        ----7h-----", side=1, line = 3, cex=1.5)
