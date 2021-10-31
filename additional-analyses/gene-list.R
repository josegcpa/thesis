library(tidyverse)
library(ggpubr)
library(cowplot)

translate <- function(x) {
  ifelse(
    grepl("\\(",x),
    ifelse(grepl("controls",x),
           "Healthy individuals",
           ifelse(grepl("AML",x),"Progressed to AML","Therapy-related")),
    "Healthy individuals"
  ) %>%
    return
}

sample_sizes <- c(
  `Jaiswal 2014` = 15801,`Genovese 2014` = 12380,`Young 2016` = 20,
  `McKerrel 2015` = 4067,`Zink 2017` = 11262,
  `Acuna-Hidalgo 2017` = 2006,`Coombs 2017 (therapy-related)` = 5649,
  `Desai 2018 (controls)` = 181,`Desai 2018 (prog. to AML)` = 188,
  `Young 2019 (controls)` = 69,`Young 2019 (prog. to AML)` = 34,
  `Bolton 2020 (therapy-related)` = 5978)

df <- read_csv("../csv/gene-list-proportions.csv") %>%
  gather(key = "Study",value = "Proportion",-name,-abbreviation,-"function") %>%
  mutate(Proportion = as.numeric(Proportion)) %>% 
  group_by(Study) %>% 
  mutate(R = rank(-Proportion,ties.method = c("first"))) %>%
  mutate(R = ifelse(is.na(Proportion),NA,R)) %>%
  ungroup %>% 
  mutate(N = sample_sizes[Study],ind = translate(Study)) %>%
  rowwise() %>% 
  mutate(Study = paste(str_split(Study,' ')[[1]][1:2],collapse = " ")) %>% 
  group_by(abbreviation) %>%
  mutate(average_rank = mean(R,na.rm=T)) %>%
  group_by(Study) %>%
  mutate(N = sum(unique(N))) %>%
  group_by(abbreviation) %>% 
  mutate(Minimum = min(Proportion,na.rm = T)) %>%
  mutate(Maximum = max(Proportion,na.rm = T)) %>%
  mutate(AverageProportion = mean(Proportion,na.rm=T))

not_studied <- df %>% 
  subset(is.na(Proportion))

heatmap_plot <- df %>%
  subset(Proportion > 0) %>% 
  ggplot(aes(x = reorder(abbreviation,average_rank),
             y = reorder(Study,N),size = Proportion,
             colour = R,shape = ind,group = ind)) + 
  geom_point() +
  geom_tile(data = not_studied,aes(x = abbreviation,y = Study), 
            inherit.aes = F,fill = "grey93") + 
  theme_classic(base_size = 6) + 
  rotate_x_text(face = "italic") + 
  scale_radius(
    trans = 'log10',
    breaks = c(1e-4,1e-2,1),labels = function(x) sprintf("%s",x),
    range = c(1,6),name = "Mutations\nper individual") +
  scale_colour_gradient(low = "yellow2",high = "red4",guide = FALSE) + 
  theme(axis.text = element_text(size = 6),
        legend.text = element_text(size = 6)) +
  theme(legend.position = "bottom",legend.key.height = unit(0.5,"cm"),
        strip.text.y = element_blank(),
        strip.background = element_blank(),axis.title.x = element_blank()) + 
  scale_shape_manual(name = NULL,values = c(16,10,5)) +
  guides(shape = guide_legend(nrow = 3),size = guide_legend(ncol = 1)) + 
  ylab("Study")

freq_plot <- df %>% 
  ungroup %>% 
  select(Study,N) %>%
  distinct %>%
  group_by(Study) %>%
  summarise(N = sum(N)) %>%
  ggplot(aes(x = reorder(Study,N),y = N)) + 
  geom_bar(stat = "identity",fill = "grey60") + 
  theme_classic(base_size = 6) + 
  coord_flip() + 
  scale_y_continuous(trans = 'log10',expand = c(0,0,0.1,0)) +
  theme(axis.text = element_text(size = 6),
        legend.text = element_text(size = 6)) +
  theme(axis.text.y = element_blank(),axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  ylab("Study size")

combined_plot <- plot_grid(heatmap_plot,freq_plot,nrow = 1,align = "hv",axis = "bt",
          rel_widths = c(1,0.4)) 

ggsave(plot = combined_plot,"ch-gene-proportions.pdf",height = 3,width = 6)

chromosome_locations <- c(
  DNMT3A = "2p23.3",TET2 = "4q24",ASXL1 = "20q11.21",CTCF = "16q22.1",
  KRAS = "12p12.1",JAK2 = "9p24.1",PTPN11 = "12q24.13",SF3B1 = "2q33.1",
  SRSF2 = "17q25.1",U2AF1 = "21q22.3",PPM1D = "17q23.2",BRCC3 = "Xq28",
  IDH1 = "2q34",IDH2 = "15q26.1",TP53 = "17p13.1",GNB1 = "1p36.33", CBL = "11q23.3"
)

df %>%
  select(name,abbreviation,f = "function",AverageProportion,Minimum,Maximum) %>% 
  distinct %>% 
  mutate(f = ifelse(grepl("Cell cycle regula",f),"Cell cycle regulation",f)) %>%
  mutate(f = ifelse(grepl("DNA dama",f),"DNA damage regulation",f)) %>%
  group_by(f) %>%
  mutate(NU = length(unique(name))) %>%
  arrange(-NU,f,-AverageProportion) %>%
  transmute(name = name,ab = abbreviation,f,p = sprintf("%.4f (%.4f-%.4f)",AverageProportion,Minimum,Maximum)) %>% 
  mutate(chr = chromosome_locations[ab]) %>%
  write.csv("../csv/gene-list-proportions-thesis.csv",quote = F,row.names = F) 
