require(hexSticker)
require(dplyr)
require(ggpubr)
require(gridExtra)
require(GenomicRanges)
require(ggbio)
require(RColorBrewer)


seed <- 110
set.seed(seed)

p <- 11
q <- 6


# hotspot_propensity <- rbeta(p, shape1 = 1.5, shape2 = 1.5)
# m <- sweep(matrix(rnorm(p*q, mean = 1, sd = 0.1), nrow = p), 1, hotspot_propensity, "*")
# m[m<0] <- 0
# m[m>1] <- 1
# m[m < 0.8] <- 0 

m <- matrix(rbeta(p*q, shape1 = 1, shape2 = 1), nrow = p)
colnames(m) <- paste("Col", 1:q)
rownames(m) <- paste("Row", 1:p)
m[m < 0.80] <- 0
m[sample(c(T,F), p, prob = c(0.5, 0.5), replace = T),] <- 0

# Transform the matrix in long format
df <- reshape::melt(m)
colnames(df) <- c("x", "y", "value")

logo_top <- ggplot(df, aes(x = x, y = y, fill = value)) +  
  scale_fill_gradientn(colours=c("white", brewer.pal(7, "Dark2"))) +
  geom_tile(color = "grey80",
            lwd = 0.2,
            linetype = 1) +
  coord_fixed()  +theme(axis.title.x=element_blank(),
                        axis.text.x=element_blank(),
                        axis.ticks.x=element_blank(),
                        axis.title.y=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks.y=element_blank(),
                        legend.position="none")


seed <- 321
set.seed(seed)

N <- 5

gr <- GRanges(seqnames = rep("chr1", N),
              IRanges(
                start = sample(1:(20*N), size = N, replace = TRUE),
                width = sample(25:50, size = N, replace = TRUE)),
              strand = sample(c("+", "-", "*"), size = N,
                              replace = TRUE),
              value = rnorm(N, 2, 20), score = rnorm(N, 100, 30),
              sample = sample(c("Normal", "Tumor"),
                              size = N, replace = TRUE),
              pair = sample(letters, size = N,
                            replace = TRUE)) 

logo_bottom <- ggplot() + geom_alignment(gr,fill = "grey80") + theme_void() + scale_y_reverse()

dir.create("man/figures/", showWarnings = FALSE)

sticker(grid.arrange(logo_top, logo_bottom, heights=c(11,2.4)), 
        package="locus", 
        p_size=6.3, 
        p_color = "grey25",
        s_x=0.97, 
        s_y=0.985, 
        s_width=1.5, 
        s_height=1.35,
        p_x = 1.355, 
        p_y = 0.65, 
        h_size = 1,
        h_fill="white", 
        h_color= "grey25", #"grey80",
        filename="man/figures/locus_logo.png",
        dpi = 1200)

