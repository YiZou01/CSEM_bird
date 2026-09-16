#Appendix data analysis and visualization for "Citizen-science and empirical field data yield limited agreement about avian diversity"  

#Last update: 2026-Sep-16
library(readxl)
library(reshape2)
library(ggplot2)
library(patchwork)
library(iNEXT)
library(vegan)

#1. Read data
# Location-month here means month of year (1-12); the same month across years is pooled.
CS <- read.csv("CSdata5km.csv")
# The current workbook calls the month column "目标月份".
if (!"month" %in% names(CS) && "目标月份" %in% names(CS)) {CS$month <- CS[["目标月份"]]}
CS$O_month <- as.integer(substr(CS$month, 6, 7))
CSag <- aggregate(number~Location+O_month+IOC_14.1, data = CS, FUN=sum)
CSag$Type <- "CS"

EM <- read.csv("EMdata.csv")
EM$Date <- as.POSIXct(EM$Date,format = "%Y/%m/%d %H:%M",tz = "Asia/Shanghai")
EM$O_month <- as.integer(format(EM$Date, "%m"))

#EM$O_month <- as.integer(format(EM$Date, "%m"))
#Sum number of individuals by species, location and month
EMag <- aggregate(number~Location+O_month+IOC_14.1, data = EM, FUN=sum)
EMag$Type <- "EM"
#Data combined
df <- rbind(CSag,EMag)
# Analysis 1. Diversity indices between two data sources
shannon_fun <- function(x) {p <- x / sum(x);-sum(p * log(p))}

shannon_long <- aggregate(number ~ Location + O_month + Type,data = df,FUN = shannon_fun)
names(shannon_long)[names(shannon_long) == "number"] <- "Shannon"

# Put the CS and EM values for each location-month in the same row
shannon <- dcast(shannon_long, Location + O_month ~ Type, value.var = "Shannon")
shannon <- shannon[order(shannon$Location, shannon$O_month), ]
# Compare only location-months observed by BOTH sources; missing is not zero.
shannon <- shannon[complete.cases(shannon[, c("EM", "CS")]), ]

#Correlation between two diversity index
# Pearson
r_h <- cor.test(shannon$EM,shannon$CS,method = "pearson")
# Spearman
s_h <- cor.test(shannon$EM,shannon$CS,method = "spearman",exact = FALSE)
# Create the text shown in the figure
cor_text_h <- sprintf(    "r = %.2f, P = %.2f\n\u03c1 = %.2f, P = %.2f",
    unname(r_h$estimate), r_h$p.value, unname(s_h$estimate), s_h$p.value)


mytheme <- theme_classic()+theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background=element_rect(fill = "white"),
                 aspect.ratio = 1)
# Keep the two-line statistics in the upper-left margin inside every panel.
h <- -0.08
v <- 1.1

# Plot
h_xmin <- min(shannon$EM, na.rm = TRUE) * 0.9
h_xmax <- max(shannon$EM, na.rm = TRUE) * 1.1
h_ymax <- max(shannon$CS, na.rm = TRUE) * 1.22
p <- ggplot(shannon, aes(x = EM, y = CS)) +
    geom_abline(slope = 1,intercept = 0,linetype = "dashed",colour = "grey50") +
    geom_smooth( method = "lm") +
    geom_point(size = 2, colour = "steelblue") +
#    ggrepel::geom_text_repel(aes(label = Location),size = 3) +
    coord_cartesian(xlim = c(h_xmin, h_xmax),ylim = c(0, h_ymax)) +
    labs(x = "Empirical survey",y = "Citizen science",title = "(a)") +
    annotate("text",x = -Inf,y = Inf,label = cor_text_h,size = 2.5,
        hjust = h, vjust = v) +
    mytheme
print(p)

#  Hill diversity at the same sample coverage
abundance_list <- split(df$number,paste(df$Location, df$O_month, df$Type, sep = "_"))
# q = 0: richness; q = 1: Shannon diversity; q = 2: Simpson diversity
hill <- estimateD(abundance_list,q = c(0, 1, 2),
                  datatype = "abundance",
                  base = "coverage",level = 0.95)
# Assemblage names have the form Location_month_Type, e.g. Foping_7_CS.
hill$Location <- sub("_[0-9]+_(CS|EM)$", "", hill$Assemblage)
hill$O_month <- as.integer(sub("^.*_([0-9]+)_(CS|EM)$", "\\1", hill$Assemblage))
hill$Type <- sub(".*_", "", hill$Assemblage)

# Put the CS and EM estimates for each location-month in the same row
hill_wide <- dcast(hill, Location + O_month + Order.q ~ Type, value.var = "qD")
# Keep only paired, finite estimates for comparisons.
hill_wide <- hill_wide[is.finite(hill_wide$EM) & is.finite(hill_wide$CS), ]

hill_q0 <- hill_wide[hill_wide$Order.q == 0, ]
hill_q1 <- hill_wide[hill_wide$Order.q == 1, ]
hill_q2 <- hill_wide[hill_wide$Order.q == 2, ]

#Correlation between two diversity index
r_q0 <- cor.test(hill_q0$EM,hill_q0$CS,method = "pearson")
s_q0 <- cor.test(hill_q0$EM,hill_q0$CS,method = "spearman",exact = FALSE)
cor_text_q0 <- sprintf(    "r = %.2f, P = %.2f\n\u03c1 = %.2f, P = %.2f",
                        unname(r_q0$estimate), r_q0$p.value, unname(s_q0$estimate), s_q0$p.value)

r_q1 <- cor.test(hill_q1$EM,hill_q1$CS,method = "pearson")
s_q1 <- cor.test(hill_q1$EM,hill_q1$CS,method = "spearman",exact = FALSE)
cor_text_q1 <- sprintf(    "r = %.2f, P = %.2f\n\u03c1 = %.2f, P = %.2f",
                        unname(r_q1$estimate), r_q1$p.value, unname(s_q1$estimate), s_q1$p.value)

r_q2 <- cor.test(hill_q2$EM,hill_q2$CS,method = "pearson")
s_q2 <- cor.test(hill_q2$EM,hill_q2$CS,method = "spearman",exact = FALSE)
cor_text_q2 <- sprintf(    "r = %.2f, P = %.2f\n\u03c1 = %.2f, P = %.2f",
                        unname(r_q2$estimate), r_q2$p.value, unname(s_q2$estimate), s_q2$p.value)


#Plot the results
q0_xmax <- max(hill_q0$EM, na.rm = TRUE) * 1.1
q0_ymax <- max(hill_q0$CS, na.rm = TRUE) * 1.22
q1_xmax <- max(hill_q1$EM, na.rm = TRUE) * 1.1
q1_ymax <- max(hill_q1$CS, na.rm = TRUE) * 1.22
q2_xmax <- max(hill_q2$EM, na.rm = TRUE) * 1.1
q2_ymax <- max(hill_q2$CS, na.rm = TRUE) * 1.22

q0 <- ggplot(hill_q0, aes(x = EM, y = CS)) +
    geom_abline(slope = 1, intercept = 0,linetype = "dashed", colour = "grey50") +
    geom_smooth(method = "lm", formula = y ~ x) +
    geom_point(size = 2, colour = "steelblue") +
    coord_cartesian(xlim = c(0, q0_xmax), ylim = c(0, q0_ymax)) +
    labs(x = "Empirical survey", y = "Citizen science",
         title = "(b)") +
    annotate("text",x = -Inf,y = Inf,label = cor_text_q0,size = 2.5,
             hjust = h, vjust = v) +
    mytheme

q1 <- ggplot(hill_q1, aes(x = EM, y = CS)) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "grey50") +
    geom_smooth(method = "lm", formula = y ~ x) +
    geom_point(size = 2, colour = "steelblue") +
    coord_cartesian(xlim = c(0, q1_xmax), ylim = c(0, q1_ymax)) +
    labs(x = "Empirical survey", y = "Citizen science",
         title = "(c)") +
    annotate("text",x = -Inf,y = Inf,label = cor_text_q1,size = 2.5,
             hjust = h, vjust = v) +
    mytheme

q2 <- ggplot(hill_q2, aes(x = EM, y = CS)) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "grey50") +
    geom_smooth(method = "lm", formula = y ~ x) +
    geom_point(size = 2, colour = "steelblue") +
    coord_cartesian(xlim = c(0, q2_xmax), ylim = c(0, q2_ymax)) +
    labs(x = "Empirical survey", y = "Citizen science",
         title = "(d) ") +
    annotate("text",x = -Inf,y = Inf,label = cor_text_q2,size = 2.5,
             hjust = h, vjust = v) +
    mytheme

Appendix1 <- print((p|q0)/(q1|q2))
print(Appendix1)
ggsave("Appendix1.jpg", Appendix1,
       width = 5, height = 5, dpi = 300)
