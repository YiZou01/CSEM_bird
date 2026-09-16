#Data analysis and visualization for "Citizen-science and empirical field data yield limited agreement about avian diversity"  
#Last update: 2026-Sep-16

#Load library
library(readxl)
library(reshape2)
library(ggplot2)
library(patchwork)
library(iNEXT)
library(vegan)

# Read data
#CS <- read("CSdata5km.xlsx",sheet="Details")
CS <- read.csv("CSdata5km.csv")
CS$O_month <- as.integer(substr(CS$month, 6, 7))
#Sum number of individuals by species and location 
CSag <- aggregate(number~Location+IOC_14.1, data = CS, FUN=sum)
CSag$Type <- "CS"

EM <- read.csv("EMdata.csv")
EM$O_month <- as.integer(format(EM$Date, "%m"))
#Sum number of individuals by species and location 
EMag <- aggregate(number~Location+IOC_14.1, data = EM, FUN=sum)
EMag$Type <- "EM"

#Data combined
df <- rbind(CSag,EMag)

# Analysis 1. Diversity indices between two data sources
shannon_fun <- function(x) {p <- x / sum(x);-sum(p * log(p))}

shannon_long <- aggregate(number ~ Location + Type,data = df,FUN = shannon_fun)
names(shannon_long)[3] <- "Shannon"

# Put the CS and EM values for each location in the same row
shannon <- dcast(shannon_long, Location ~ Type, value.var = "Shannon")
shannon <- shannon[order(shannon$Location), ]

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
                 panel.background=element_rect(fill = "white"))
h <- 1.1
v <- -0.2

# Plot: empirical diversity vs citizen-science diversity
p <- ggplot(shannon, aes(x = EM, y = CS)) +
    geom_abline(slope = 1,intercept = 0,linetype = "dashed",colour = "grey50") +
    geom_smooth( method = "lm") +
    geom_point(size = 2, colour = "steelblue") +
    coord_equal(xlim = c(1.5, 5),ylim = c(1.5, 5)) +
    labs(x = "Empirical survey",y = "Citizen science",title = "(a)") +
    annotate("text",x = Inf,y = -Inf,label = cor_text_h,size = 2.5,
        hjust = h, vjust = v) +
    mytheme
print(p)

#  Hill diversity at the same sample coverage
abundance_list <- split(df$number,paste(df$Location, df$Type, sep = "_"))
# q = 0: richness; q = 1: Shannon diversity; q = 2: Simpson diversity
hill <- estimateD(abundance_list,q = c(0, 1, 2),
                  datatype = "abundance",
                  base = "coverage",level = 0.95)
hill$Location <- sub("_.*", "", hill$Assemblage)
hill$Type <- sub(".*_", "", hill$Assemblage)

# Put the CS and EM estimates for each location in the same row
hill_wide <- dcast(hill, Location + Order.q ~ Type, value.var = "qD")

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
q0_max <- max(hill_q0$EM, hill_q0$CS,na.rm = T) *1.1
q1_max <- max(hill_q1$EM, hill_q1$CS,na.rm = T) *1.1
q2_max <- max(hill_q2$EM, hill_q2$CS,na.rm = T) *1.1

q0 <- ggplot(hill_q0, aes(x = EM, y = CS)) +
    geom_abline(slope = 1, intercept = 0,linetype = "dashed", colour = "grey50") +
    geom_smooth(method = "lm", formula = y ~ x) +
    geom_point(size = 2, colour = "steelblue") +
    coord_equal(xlim = c(0, q0_max), ylim = c(0, q0_max)) +
    labs(x = "Empirical survey", y = "Citizen science",
         title = "(b)") +
    annotate("text",x = Inf,y = -Inf,label = cor_text_q0,size = 2.5,
             hjust = h, vjust = v) +
    mytheme

q1 <- ggplot(hill_q1, aes(x = EM, y = CS)) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "grey50") +
    geom_smooth(method = "lm", formula = y ~ x) +
    geom_point(size = 2, colour = "steelblue") +
    coord_equal(xlim = c(0, q1_max), ylim = c(0, q1_max)) +
    labs(x = "Empirical survey", y = "Citizen science",
         title = "(c)") +
    annotate("text",x = Inf,y = -Inf,label = cor_text_q1,size = 2.5,
             hjust = h, vjust = v) +
    mytheme

q2 <- ggplot(hill_q2, aes(x = EM, y = CS)) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "grey50") +
    geom_smooth(method = "lm", formula = y ~ x) +
    geom_point(size = 2, colour = "steelblue") +
    coord_equal(xlim = c(0, q2_max), ylim = c(0, q2_max)) +
    labs(x = "Empirical survey", y = "Citizen scienc",
         title = "(d) ") +
    annotate("text",x = Inf,y = -Inf,label = cor_text_q2,size = 2.5,
             hjust = h, vjust = v) +
    mytheme

p1 <- print((p|q0)/(q1|q2))
ggsave("Figure 2.jpg", p1,
       width = 5, height = 5, dpi = 300)


# Analysis 2, species composition
# Show the top 10 most abundant species separately for CS and EM
species_abundance <- aggregate(number ~ IOC_14.1 + Type,data = df,FUN = sum)
species_abundance$Proportion <- with(species_abundance, number / ave(number, Type, FUN = sum))

nth = 10
# CS data
top10_CS <- species_abundance[species_abundance$Type == "CS", ]
top10_CS <- head(top10_CS[order(-top10_CS$Proportion), ], nth)
top10_CS$IOC_14.1 <- factor(top10_CS$IOC_14.1,
                            levels = top10_CS$IOC_14.1)

# EM data
top10_EM <- species_abundance[species_abundance$Type == "EM", ]
top10_EM <- head(top10_EM[order(-top10_EM$Proportion), ], nth)
top10_EM$IOC_14.1 <- factor(top10_EM$IOC_14.1,
                            levels = top10_EM$IOC_14.1)
#common species
common_species <- intersect(as.character(top10_CS$IOC_14.1),
    as.character(top10_EM$IOC_14.1))
add_star <- function(x) {ifelse(x %in% common_species, paste0("*",x), x)}

# Plot the top 10 species
y_max <- max(top10_CS$Proportion, top10_EM$Proportion) * 1.1

barplot_CS <- ggplot(top10_CS, aes(x = IOC_14.1, y = Proportion)) +
    geom_col(fill = "#3B82B8", width = 0.7) +
    geom_text(aes(label = scales::percent(Proportion, accuracy = 0.1)),
              vjust = -0.4, size = 3) +
    scale_y_continuous(labels = scales::label_percent(),
                       limits = c(0, y_max)) +
    scale_x_discrete(labels = add_star)+
    labs(x = NULL, y = "Proportion of total abundance", title = "(a)") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                     face = "italic"))

barplot_EM <- ggplot(top10_EM, aes(x = IOC_14.1, y = Proportion)) +
    geom_col(fill = "#E68A3F", width = 0.7) +
    geom_text(aes(label = scales::percent(Proportion, accuracy = 0.1)),
              vjust = -0.4, size = 3) +
    scale_y_continuous(labels = scales::label_percent(),
                       limits = c(0, y_max)) +
    scale_x_discrete(labels = add_star)+
    labs(x = NULL, y = "", title = "(b)") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                     face = "italic"))

# Dissimilarity and NMDS
community <- dcast(df, Location + Type ~ IOC_14.1,
                   value.var = "number", fun.aggregate = sum, fill = 0)

community_info <- community[, c("Location", "Type")]
community_matrix <- community[, -(1:2)]
rownames(community_matrix) <- paste(community_info$Location,
                                    community_info$Type, sep = "_")

# Relative abundance reduces the influence of unequal total sample sizes
community_matrix <- decostand(community_matrix, method = "total")


nmds <- metaMDS(community_matrix, distance = "bray", k = 2,
                trymax = 100, autotransform = FALSE, trace = FALSE)

nmds_scores <- as.data.frame(scores(nmds, display = "sites"))
nmds_scores$Location <- community_info$Location
nmds_scores$Type <- community_info$Type

NMDSplot <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
    stat_ellipse(aes(colour = Type), geom = "polygon",
                 alpha = 0.12, linewidth = 0.7,show.legend = FALSE) +
    geom_line(aes(group = Location), colour = "grey70", linewidth = 0.5) +
    geom_point(aes(colour = Type), size = 3) +
    scale_colour_manual(values = c(CS = "#3B82B8", EM = "#E68A3F"),
                        labels = c(CS = "Citizen science",EM = "Empirical survey")) +
    labs(x = "NMDS1", y = "NMDS2",colour = "",title = "(c)") +
    mytheme+ 
theme(
    legend.position = c(0.15,1.10),
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.justification = "left",
    legend.box.margin = margin(0, 0, -2, 0),
    aspect.ratio = 0.80)+
    guides(colour = guide_legend(nrow = 1))
nmds$stress
print(NMDSplot)


# Analysis 3. Record frequency versus abundance
CS_freq <- aggregate(number ~ IOC_14.1, data = CS, FUN = sum)
names(CS_freq)[2] <- "Abundance"
CS_freq$Frequency <- as.numeric(table(CS$IOC_14.1)[CS_freq$IOC_14.1])
CS_freq$Source <- "CS"

EM_freq <- aggregate(number ~ IOC_14.1, data = EM, FUN = sum)
names(EM_freq)[2] <- "Abundance"
EM_freq$Frequency <- as.numeric(table(EM$IOC_14.1)[EM_freq$IOC_14.1])
EM_freq$Source <- "EM"

freq_abund <- rbind(CS_freq, EM_freq)

AF_plot <- ggplot(freq_abund,
             aes(y = Abundance, x = Frequency, colour = Source)) +
    geom_point(alpha = 0.45, size = 2) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 1) +
    scale_x_log10(labels = scales::label_comma()) +
    scale_y_log10(labels = scales::label_comma()) +
    scale_colour_manual(values = c(CS = "#3B82B8", EM = "#E68A3F"),
                        labels = c(CS = "Citizen science",EM = "Empirical survey")) +
    labs(y = "Total abundance",
         x = "Record frequency",title = "(d)",
         colour = " ") +
    theme_classic(base_size = 12) +
    theme(
        legend.position = c(0.05,0.9),
        legend.background = element_blank(),
        legend.justification = "left",
        legend.box.margin = margin(0, 0, -2, 0),
        aspect.ratio = 0.80)
AF_plot

p2 <- print((barplot_CS|barplot_EM)/(NMDSplot|AF_plot))
p2
ggsave("Figure 3.jpg", p2,
       width = 8, height = 8, dpi = 300)

#Model, frequency vs abundance data
mod <- lm(log10(Abundance) ~ log10(Frequency) * Source,data = freq_abund)
summary(mod) # 
mod1 <- lm(log10(Abundance) ~ log10(Frequency) ,data = subset(freq_abund,Source == "CS"))
summary(mod1)
mod2 <- lm(log10(Abundance) ~ log10(Frequency) ,data = subset(freq_abund,Source == "EM"))
summary(mod2)
