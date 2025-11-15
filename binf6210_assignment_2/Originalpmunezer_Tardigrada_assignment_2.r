##*******************************
## BINF*6210 - BIOINFORMATICS TOOLS
##
## Instructor: Prof. Karl Cottenie
##
## ASSIGNMENT 1: UNEVEN GLOBAL SAMPLING IN BOLD AFFECTS TARDIGRADA GENETIC DIVERSITY
## 
## Pierre Celestin Munezero
## 
## 2025-10-13
##
##*******************************
##
## _ Packages to be used ---------------------

library(stats)
library(tidyverse)
library(viridis)
library(vegan)
library(VennDiagram)
library(ggpubr)

# + scale_color/fill_viridis(discrete = T/F)
theme_set(theme_light())
conflicted::conflict_prefer("filter", "dplyr")

## _ Load Tardigrada_BOLD_data.tsv ---------

# Check the current directory
getwd()

## Load the tsv file
df_tardi <- read_tsv("../data/Tardigrada_BOLD_data.tsv")
df_tardi

# Check statistical summaries
class(df_tardi)
dim(df_tardi) #==> 5441_____93
names(df_tardi)

# In this study, I will need 3 variables:  bin_uri, species, country/ocean
# Let's edit the name of a variable with division operator, country/ocean as a vector

names.original <- as.vector(names(df_tardi))
names.original

# Replace "/" with "_"
names.edited <- str_replace(string = names.original, pattern = "/", replacement = "_")
names.edited

# Let's check the summary
summary(df_tardi) 

# Let's change the names in the original data frame
names(df_tardi) <- names.edited
names(df_tardi)

# Create a subset of variable needed for this project

df_tardi.sub <- df_tardi[, c("bin_uri", "species", "country_ocean")]
df_tardi.sub

# Let's check the new dataset
dim(df_tardi.sub)
#==> 5441______3

## _ Objective 1: Quantify specimen records, species, and BINs per country ----

## __ Descriptive analysis of species richness per country ------

# Create new object for descriptive analysis per country with no missing values

df_tardi.sub_Desc <- df_tardi.sub %>% 
  filter(!is.na(df_tardi.sub$species) & !is.na(df_tardi.sub$country_ocean)) %>% 
  group_by(country_ocean) %>% 
  summarize(records = n(), species_richness = length(unique(species)), bin_richness = length(unique(bin_uri))) %>% 
  arrange(desc(records)) %>% 
  print()

# Check how many row are remaining
nrow(df_tardi.sub_Desc)
#==> 51 countries/regions. At this stage I decided to keep "Unrecoverable" value because it will help us understand that the metadata lacks some important data

total_records <- sum(df_tardi.sub_Desc$records)
total_records
#==> Higher number of records belong to "Unrecoverable". Maybe these data were extracted from different database like GenBank where study sites were not specified.

# Let's build a ggplt to visualize the data
# Reorder countries base on species richness per country using fct_reorder() function to reorder the level of all values in country_ocean variable

ggplot(data = df_tardi.sub_Desc) +
  geom_col(mapping = aes(x = fct_reorder(country_ocean, species_richness), y = species_richness),  fill = "steelblue") +
  labs(title = "Tardigrada species-richness per country", x = "Country sampled", y = "Number of species") +
  coord_flip()

# Tardigrada species richness for the top 20 countries
# Arrange() function will sort row, this time in descending order and head() will take only top 20
df_tardi.sub_Desc_top20 <- df_tardi.sub_Desc %>% 
  arrange(desc(species_richness)) %>% 
  head(n = 20)

# Let's now plot the top 20 countries/regions
ggplot(data = df_tardi.sub_Desc_top20) +
  geom_col(mapping = aes(x = fct_reorder(country_ocean, species_richness), y = species_richness),  fill = "steelblue") +
  labs(title = "Tardigrada species-richness for top 20 countries", x = "Country sampled", y = "Number of species") +
  coord_flip()

# Find BIN richness and species richness using the original data frame removing missing values

dim(df_tardi)
df_tardi.Desc <- df_tardi %>% 
  filter(!is.na(df_tardi$species) & !is.na(df_tardi$country_ocean)) %>% 
  group_by(country_ocean) %>% 
  summarize(records = n(), species_richness = length(unique(species)), bin_richness = length(unique(bin_uri))) %>% 
  arrange(desc(records)) %>% 
  print()

nrow(df_tardi.sub_Desc)

# Check both tibbles
all.equal(df_tardi.sub_Desc, df_tardi.Desc)

#==> Return TRUE

total_records <- sum(df_tardi.sub_Desc$records)
total_records

## __ Sampling effort and species richness -------------
# Let's start with rarefaction curve. The curve will show the number of BINs per level of sampling sizes. I want to see how the number of species increases as the sampling effort increases

# Let use the original data frame

names(df_tardi)
dim(df_tardi)

# To count the number of records in each unique BIN, let's group data per BINs

df_tardi.count.by.bin <- df_tardi %>% 
  group_by(bin_uri) %>% 
  count(bin_uri)

# Let's convert the data in community data object format by spreading the data frame

df_tardi.BINs.spread <- pivot_wider(data = df_tardi.count.by.bin, names_from = bin_uri, values_from = n)

# Let's build a rarefaction curve

x <- rarecurve(df_tardi.BINs.spread, xlab = "Individual barcoded", ylab = "BIN richness")
#==> The curve shows that the number of new BINs continue to increase. This means that many new BINs are still being discovered with each additional sample. It suggests that the sampling effort is not yet sufficient to capture the full diversity of Tardigrada.

# Let's now draw the Species Accumulation curve. This will check the BIN richness by country sampled instead of individual sample.

df_tardi

## Let's remove missing values before grouping
df_tardi.BINs.by.country <- df_tardi %>% 
  filter(!is.na(country_ocean), !is.na(bin_uri)) %>% 
  group_by(country_ocean, bin_uri) %>% 
  count(country_ocean, bin_uri)


# Let's convert the data frame into community format

df_tardi.BINs.spread.by.country <- pivot_wider(data = df_tardi.BINs.by.country, names_from = bin_uri, values_from = n)

# Make sure the country name are stored as strings not factors
df_tardi.BINs.spread.by.country %>% 
  mutate(country_ocean = trimws(as.character(country_ocean))) %>%
  column_to_rownames(var = "country_ocean")

# check the data frame
class(df_tardi.BINs.by.country)
str(df_tardi.BINs.by.country)

# For downstream analysis, let's convert NA, which represent the absence of barcode in a given country, into 0 for no barcode recorded
df_tardi.BINs.spread.by.country[is.na(df_tardi.BINs.spread.by.country)] <- 0

# Let's build the species accumulation curve (SAC)

# Check any duplicated country
any(duplicated(df_tardi.BINs.spread.by.country$country_ocean))

#==> FALSE

# Check the type of df_tardi.BINs.spread.by.country
class(df_tardi.BINs.spread.by.country)

# Let's convert the tible into a dat frame and change the "country_ocean" into row name

df_tardi.BINs.spread.by.country_com <- df_tardi.BINs.spread.by.country %>% 
  as.data.frame() %>% 
  column_to_rownames(var = "country_ocean")

# Let's build the species accumulation curve
AccumCurve <- specaccum(df_tardi.BINs.spread.by.country_com)

# Plot the curve

plot(AccumCurve, main = "Tardigrada Species Accumulation Curve", xlab = "Number of countries sampled", ylab = "BIN richness")

#==> The SAC keeps rising steadly without reaching a plateau. This means that additional sampling is necessary as new species are still being discovered with each additional country sampled. This SAC suggests that the sampling effort has not fully captured the global diversity of Tardigrada

## _ Objective 2: To compare well_sampled countries to underrepresented countries -----------------
## Here I want to test if well-sampled countries have higher species richness

# To avoid arbitrary threshold, let's use a threshold of 75th percentile, which is a practical way to classify countries based on sampling effort. 75th percentile means that the top 25% of the sampled countries with the highest number of records will be considered as well-sampled countries in this dataset. Reference is made to the study by Skoulikidis et al. 2025 (https://doi.org/10.1016/j.ecolind.2025.114071), which use the 75th percentile of nutrient concentration to set data-driven threshold for classifying sites as pristine.

# Let's prepare a new dataset excluding "unrecoverable" as one value in the variable country_ocean. This is to avoid its influence on the classification of well-sampled countries.

head(df_tardi)

df_tardi.countries <- df_tardi %>% 
  filter(country_ocean != "Unrecoverable")

# Let's check how many rows remained
dim(df_tardi.countries)
#==> 3012----93

## __ Classification of countries as well-sampled and underrepresented countries ---------

# Let's create a new dataset of number of records, species, and BINs per country

names(df_tardi.countries) # To check the data set variables
df_tardi.countries.sub <- df_tardi.countries %>% 
  filter(!is.na(df_tardi.countries$species) & !is.na(df_tardi.countries$country_ocean)) %>% 
  group_by(country_ocean) %>% 
  summarize(records = n(), species_richness = length(unique(species)), bin_richness = length(unique(bin_uri))) %>% 
  arrange(desc(records)) %>% 
  print()

# Let's calculate the threshold. quantile() function will found the values above which 25% of the data fall.

threshold <- quantile(df_tardi.countries.sub$records, 0.75, na.rm = TRUE)
threshold
#==> The threshold is 17 records. All countries with >= 17 records will be classified as well-sampled countries and those with < 17 as underrepresented regions.

# Let's classify the countries in well-sampled and underrepresented regions

df_tardi.countries.sub_class <- df_tardi.countries.sub %>% 
  mutate(region_class = ifelse(records >= 17, "Well-sampled", "Underrepresented"))

# Let's see the summary

summary(df_tardi.countries.sub_class)
names(df_tardi.countries.sub_class)
length(names(df_tardi.countries.sub_class))

# library(ggpubr) will help display the statistical significance on the figure 

ggplot(data = df_tardi.countries.sub_class, aes(x = region_class, y = species_richness)) +
  geom_boxplot(outlier.color = "red", outlier.fill = "red", outlier.size = 2) +
  stat_compare_means(method = "wilcox.test", label = "p.signif" )+
  labs(title = "Species richness in well-sampled vs underrepresented countries", x = "Region classification", y = "Species richness") 

# Let's run ggplot without statistical test 
ggplot(data = df_tardi.countries.sub_class) +
  geom_boxplot(mapping = aes(x = region_class, y = species_richness),  outlier.color = "red", outlier.fill = "red", outlier.size = 2) +
  labs(title = "Species richness in well-sampled vs underrepresented countries", x = "Region classification", y = "Species richness")
## ggplot2 can't add directly wilcox.test()

# Wilcoxon rank-sum test is a non-paramatric test for comparison. With exact = FALSE, R will estimate the p-value using normal approximation.
?wilcox.test()
wilcox_test_result <- wilcox.test(species_richness ~ region_class, data = df_tardi.countries.sub_class, exact = FALSE)
wilcox_test_result
#==> W = 33, p-value = 1.59e-06 (p < 0.001)

# Let's check both codes
all.equal(ggplot(data = df_tardi.countries.sub_class, aes(x = region_class, y = species_richness)), ggplot(data = df_tardi.countries.sub_class))
#==> TRUE

## __ Shannon Diversity Index ---------------

# Let's calculate the shannon diversity index. This will help calculate BIN diversity per country (alpha-diversity). Shannon index score indicates both species richness and evenness per each country. The argument .groups = "drop" tells R not to keep grouping after summarization.

df_tardi.countries_shannon <- df_tardi.countries %>% 
  filter(!is.na(country_ocean) & !is.na(bin_uri)) %>% 
  group_by(country_ocean, bin_uri) %>% 
  summarize(count = n(), .groups = "drop") %>% 
  group_by(country_ocean) %>% 
  summarize(shannon = diversity(count, index = "shannon")) %>% 
  arrange(desc(shannon)) %>%
  print()

# Let's plot these results
ggplot(data = df_tardi.countries_shannon) +
  geom_point(mapping = aes(x = fct_reorder(country_ocean, shannon), y = shannon),  size = 3) +
  labs(title = "Shannon diversity of Tardigrada across countries", x = "Country sampled", y = "Shannon Diversity Index ") +
  coord_flip()
#==> Shannon diversity index doesn't clearly give information on how BINs are evenly distributed.

# Let's look at the shannon diversity and the evenness together. This will show how evenly distributed, normalized by total richness. Evenness come to complement shannon showing how uniformly BINS are presented
df_tardi.countries_shannon_evenness <- df_tardi.countries %>% 
  filter(!is.na(country_ocean) & !is.na(bin_uri)) %>% 
  group_by(country_ocean, bin_uri) %>% 
  summarize(count = n(), .groups = "drop") %>% 
  group_by(country_ocean) %>% 
  summarize(shannon = diversity(count, index = "shannon"), richness = n_distinct(bin_uri), evenness = ifelse(richness > 1, shannon / log(richness), NA), total_records = sum(count)) %>% 
  arrange(desc(shannon)) %>%
  print()

#==> All countries with shannon diversity = 0 returned evenness to NA. For visualization, these rows are first removed.

df_tardi.countries_shannon_evenness.clean <- df_tardi.countries_shannon_evenness %>% 
  filter(!is.na(evenness), !is.na(shannon))

# Before visualization of relationship, let's check if Shannon and evenness data are normally distributed

# Histogram for shannon scores
hist(df_tardi.countries_shannon_evenness.clean$shannon, main = "Histogram presentation of Shannon Diversity", xlab = "Shannon Index", col = "lightblue")
# Or let's use density plot to see well the skewness of the data
ggplot(data = df_tardi.countries_shannon_evenness.clean) +
  geom_density(mapping = aes(x = shannon, fill = "skyblue", alpha = 0.5)) +
  labs(title = "Density plot of Shannon Diversity")

#==> Not normally distributed

# Histogram for evennes 
hist(df_tardi.countries_shannon_evenness.clean$evenness, main = "Histogram presentation of Shannon Diversity", xlab = "Evenness", col = "steelblue")
# Let's see skewnwss with density plot 
ggplot(data = df_tardi.countries_shannon_evenness.clean) +
  geom_density(mapping = aes(x = evenness, fill = "skyblue", alpha = 0.5)) +
  labs(title = "Density plot of Evenness")
#==> Not normally distributed

# Let's confirm normality using a statistical test, Shapiro-Wilk normality test

shapiro.test(df_tardi.countries_shannon_evenness.clean$shannon)
shapiro.test(df_tardi.countries_shannon_evenness.clean$evenness)

#==> Both Shannon and Evenness are statistically deviated from normality (p < 0.05). Non-parametric tests are preferred.

# Visualize the relationship between shannon diversity and sampling effort (no of records per country). Here I want to see if higher diversity reflect more sampling. The geom_smooth() will run specific relationship between x and y axis. Because I have numerical variable I chose to use linear regression, so the argument method = lm means linear model to give a straight line of best fit: y = a + bx. Argument se = TRUE come to add the confidence interval around the fitted line. 

ggplot(data = df_tardi.countries_shannon_evenness.clean, aes(x = total_records, y = shannon)) +
  geom_point(color = "steelblue", size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "darkorange", linetype = "dashed") +
  labs(title = "Relationship between sampling effort and Shannon diversity", x = "Total records per country", y = "Shannon Diversity Index")

#==> Positive correlation

# Let's visualize the relationship between evenness and sampling effort
ggplot(data = df_tardi.countries_shannon_evenness.clean, aes(x = total_records, y = evenness)) +
  geom_point(color = "steelblue", size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "darkorange", linetype = "dashed") +
  labs(title = "Relationship between sampling effort and data distribution", x = "Total records per country", y = "Data distribution")
  
#==> Negative correlation

# Let's see a comparison between shannon diversity and evenness. 

# Because the data are not normally distributed, let's use non-paremetric test, Spearman correlation.

cor_test <- cor.test(df_tardi.countries_shannon_evenness.clean$shannon, df_tardi.countries_shannon_evenness.clean$evenness, method = "spearman")
cor_test

# Let's visualize this Spearman correlation
ggplot(data = df_tardi.countries_shannon_evenness.clean, aes(x = shannon, y = evenness)) +
  geom_point(color = "steelblue", size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "darkorange", linetype = "dashed") +
  labs(title = "Relationship between Shannon diversity and data distribution", x = "Total records per country", y = "Data distribution")

# Let's add shannon and evenness results to class dataset (df_tardi.countries.sub_class)
df_tardi.countries.sub_class_shannon <- df_tardi.countries.sub_class %>% 
  left_join(df_tardi.countries_shannon_evenness.clean, by = c("country_ocean"))
df_tardi.countries.sub_class_shannon

# Exclude missing values before visualization.
df_tardi.countries.sub_class_shannon.clean <- df_tardi.countries.sub_class_shannon %>% 
  filter(!is.na(evenness), !is.na(shannon))

# Let's now visualize the mean comparison of Shannon index between well-sampled and underrepresented regions 
ggplot(data = df_tardi.countries.sub_class_shannon.clean) +
  geom_boxplot(mapping = aes(x = region_class, y = shannon, fill = region_class)) +
  labs(title = "Shannon diversity in well-sampled vs underrepresented regions", x = "Region class", y = "Shannon Diversity Index")

# Let's find the p-value
wilcox.test(shannon ~ region_class, data = df_tardi.countries.sub_class_shannon.clean)

#==> W = 80, p-value = 0.001475

# Visualize the mean comparison Shannon diversity between well-sampled and underrepresented regions with the p-value
ggplot(data = df_tardi.countries.sub_class_shannon.clean, aes(x = region_class, y = shannon, fill = region_class)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test", label = "p.signif" )+
  labs(title = "Alpha-diversity in well-sampled vs underrepresented regions", x = "Region class", y = "Shannon Diversity Index")

## _ Objective 3: Identify potential hidden genetic diversity -------------

## Here I'm going to use the operator "%in%" because I want to compare something with a set of values not with one value where "==" should be the best choice

# Check the tibble to use

head(df_tardi.countries)

# Prepare a new dataset with no missing values

df_tardi.countries_clean <- df_tardi.countries %>% 
  filter(!is.na(country_ocean) & !is.na(bin_uri)) %>%
  print()

# Let's check dimensions
dim(df_tardi.countries_clean)

## __ Number of new unique BINs in underrepresented regions ----------

# Let's find unique BINs present in well-sampled countries

df_tardi.countries.sub_well <- df_tardi.countries_clean %>% 
  filter(country_ocean %in% df_tardi.countries.sub_class$country_ocean[df_tardi.countries.sub_class$region_class == "Well-sampled"]) %>% 
  pull(bin_uri) %>% 
  unique() %>% 
  print()
  
# Let's find unique BINs present in underrepresented countries
df_tardi.countries.sub_under <- df_tardi.countries_clean %>% 
  filter(country_ocean %in% df_tardi.countries.sub_class$country_ocean[df_tardi.countries.sub_class$region_class == "Underrepresented"]) %>% 
  pull(bin_uri) %>% 
  unique() %>% 
  print()
  

#Let's find new unique BINS in underrepresented regions that are not find in well-sampled regions. setdiff(x, y) identifies elements that are in x but not in y.

unique.BINs_under <- setdiff(df_tardi.countries.sub_under, df_tardi.countries.sub_well)
cat("Number of unique BINs in underrepresented regions:", length(unique.BINs_under), "\n")

#==> Number of unique BINs in underrepresented regions: 112 

## __ Venn Diagram ---------

## Venn Diagram will shoow shared BINs and unique BINs in each class.

# Let's create the Venn diagram as a file

# x = list(Underrepresented = df_tardi.countries.sub_under, Well_sampled = df_tardi.countries.sub_well, category_names = c("Underrepresented", "Well-sampled"), filename = "BINs_VennDiagram.png")

# Let's save the Venn Diagram as a grid object and not as a file

df_tardi.venn_plot <- venn.diagram(x = list(Underrepresented = df_tardi.countries.sub_under, Well_sampled = df_tardi.countries.sub_well), category_names = c("Underrepresented", "Well-sampled"), filename = NULL, fill = c("orange", "skyblue"), alpha = 0.5, lty = 1, col = c("orange", "skyblue"), cex = 1.5, cat.cex = 1.2, cat.pos = c(180, 0))
df_tardi.venn_plot

# Draw the diagram
# Let's first create a new page for the giagram
grid.newpage()

# Let's dram the diagram using the grid object
grid.draw(df_tardi.venn_plot)
grid.text("Overlap and unique BINs in well-sampled and underrepresented regions", x = unit(0.5, "npc"), y = unit(0.95, "npc"), gp = gpar(fontosize = 40))





