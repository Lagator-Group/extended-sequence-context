### PREAMBLE ####
setwd("~/Documents/PhD/mutational_spectra/R_mutational_spectra")

library(seqinr)
library(stringi)
library(tidyverse)
library(readxl)
library(magrittr)
library(cowplot)
library(scales)
library(ggtext)

groups = c('proofreading(+) MMR(+)', 'proofreading(+) MMR(-)', 'proofreading(-) MMR(+)', 'proofreading(-) MMR(-)') %>% as_factor()
spectrum = c('AT>GC', 'GC>AT', 'AT>CG', 'GC>TA', 'AT>TA', 'GC>CG') %>% as_factor()

pal_spec = c('#7E3F8F', '#EE81EC', '#6E9BF8', '#53B74C', '#D9893B', '#B3A033')
names(pal_spec) = c('AT>GC', 'GC>AT', 'AT>CG', 'GC>TA', 'AT>TA', 'GC>CG')

pal_bases = c('#B2CB94', '#DDA97C', '#7DA1D0', '#F4DB85')
names(pal_bases) = c('A', 'T', 'G', 'C')


# imports genome sequences sequence
mg1655 = read.fasta('genomes/mg1655_NC_000913.3.fasta') %>%
  unlist() %>% unname() %>% toupper() %>% factor(levels = c('A', 'T', 'G', 'C'))

gen.length = length(mg1655)


# creates function that reverse complements sequences
comp = c('A' = 'T', 'T' = 'A', 'G' = 'C', 'C' = 'G')

rev.comp = function(sequence) {
  if (is.factor(sequence)) {
    original_levels = levels(sequence)
    out = comp[as.character(rev(sequence))]
    return(factor(out, levels = original_levels))
  }
  
  if (is.character(sequence)) {
    if (length(sequence) == 1) {
      bases = strsplit(sequence, '')[[1]]
      return(paste0(rev(comp[bases]), collapse = ''))
    } else {
      return(unname(comp[rev(sequence)]))
    }
  }
  
  warning("Input must be a character vector or factor of bases.")
  return(NULL)
}


# creates function that converts numbers to pretty sci notation for plots
sci_label <- function(x, digits = 1) {
  if (all(is.na(x))) return(NA_character_)
  sapply(x, function(val) {
    if (is.na(val)) return(NA_character_)
    if (val == 0) return("")
    # if (val %in% c(1.5e-10, 0.5e-10, 15e-9, 5e-9)) return('')
    # if (val < 10) return(paste(paste0(rep('&nbsp;', 4 + digits), collapse = '') , as.character(val)))  
    sci <- format(val, scientific = TRUE)    # e.g. "4e+05"
    if (grepl("e", sci)) {
      parts <- strsplit(sci, "e")[[1]]
      base <- formatC(as.numeric(parts[1]), format = "f", digits = digits)
      exp  <- as.integer(parts[2])
      paste0(base, "×10<sup>", exp, "</sup>")
    } else {
      formatC(val, format = "f", digits = digits)
    }
  })
}



#### other datasets ####

# imports Lee 2012 data
all_data = tibble(source = 'lee_2012', genome = 'mg1655', read_xlsx('mut_accumulation_data/Lee 2012.xlsx')) %>%
  mutate(line = as.character(line))

# Foster 2015
data = tibble()
x = excel_sheets('mut_accumulation_data/Foster 2015.xlsx')
for (i in 1:length(x)) {
  data = rbind(data, tibble(source = 'foster_2015', strain = x[i],
                            read_excel("mut_accumulation_data/Foster 2015.xlsx", sheet = x[i])))
}


data %<>%
  filter(! strain %in% c('WT_LB')) %>% # filters out WT (already in dataset from Lee 2012)
  filter(! strain %in% c('ED1a', 'IAI1')) %>% # also gets rid of non-MG1655 strains
  mutate(genome = case_when(
    strain %in% c('ED1a', 'IAI1') ~ tolower(strain),
    T ~ 'mg1655'),
    strain = case_when(strain == 'WT_MM' ~ 'WT', T ~ strain),
    condition = case_when(strain == 'WT' ~ 'MM', T ~ 'LB')) %>%
  filter(consensus == ref) %>%
  # select(-consensus, -line) %>%
  rename(pos = position) %>%
  mutate(code = paste(strain, '-', condition, sep = ''))

all_data %<>%
  bind_rows(data)

# foster 2018
codes = c('mutL', 'mutL', 'mutS_mutL', 'mutS', 'mutS', 'mutS', 'mutH', 'mutS_mutL_mutH', 'mutS_mutL_mutH', 'uvrD', 'mutL_umuDC_dinB', 'mutL_mutY', 'mutL_mfd', 'mutS_mfd', 'mutL_ndk', 'mutL-MM', 'mutS-MM', 'mutL_mutY-MM', 'mutL-buffLB', 'mutS-dilLB', 'mutS-suppMM', 'mutS-28dC')
names(codes) = as.character(c(144, 288, 304, 342, 343, 555, 197, 567, 568, 191, 118, 137, 294, 368, 666, 'mm5', 'mm343', 'mm137', 'blb144', 'dlb343', 'smm343', '28dc342'))

data = read_csv('mut_accumulation_data/Foster 2018.csv') %>%
  mutate(genome = 'mg1655', source = 'foster_2018') %>%
  # filter(consensus == ref) %>%
  separate(code, into = c('code', 'line'), sep = '-') %>%
  mutate(strain = codes[code]) %>%
  separate(strain, into = c('strain', 'condition'), sep = '-') %>%
  mutate(condition = case_when(is.na(condition) ~ 'LB', T ~ condition))

data %<>% filter(strain != 'mutL_ndk') # take this out cause does match mg1655 ref seq??
# also ndk messes with the spectrum anyways


all_data %<>%
  bind_rows(data)

# Long 2016
data = rbind(tibble(strain = 'WT', read_xlsx('mut_accumulation_data/Long 2016.xlsx', sheet = 'S2', skip = 1, col_types = c('text'))),
             tibble(strain = 'WT', read_xlsx('mut_accumulation_data/Long 2016.xlsx', sheet = 'S3', skip = 1, col_types = c('text'))),
             tibble(strain = 'mutS', read_xlsx('mut_accumulation_data/Long 2016.xlsx', sheet = 'S6', skip = 1, col_types = c('text'))),
             tibble(strain = 'mutS', read_xlsx('mut_accumulation_data/Long 2016.xlsx', sheet = 'S7', skip = 1, col_types = c('text'))),
             tibble(strain = 'mutY', read_xlsx('mut_accumulation_data/Long 2016.xlsx', sheet = 'S15', skip = 1, col_types = c('text'))),
             tibble(strain = 'mutY', read_xlsx('mut_accumulation_data/Long 2016.xlsx', sheet = 'S16', skip = 1, col_types = c('text'))))

data %<>%
  mutate(source = 'long_2016',
         genome = 'long_2016',
         pos = sub('NC_000913:', '', Position) %>% as.numeric(),
         ref = substring(Mutation, 1, 1),
         alt = substring(Mutation, 3, 3),
         alt = case_when(
           substring(ref, 1, 1) == '+' ~ 'XX',
           substring(ref, 1, 1) == '-' ~ 'X',
           T ~ alt
         ),
         ref = case_when(
           substring(ref, 1, 1) == '+' ~ 'X',
           substring(ref, 1, 1) == '-' ~ 'XX',
           T ~ ref
         )) %>%
  rename(condition = norofloxacin_ng_mL,
         line = Lines) %>%
  mutate(condition = case_when(condition == '0' ~ 'LB', T ~ condition)) %>%
  select(source, strain, genome, pos, ref, alt, line, condition) %>%
  mutate(code = paste(strain, '-', condition, '-long', sep = ''))

data %<>% filter(!condition %in% c('3.125')) # excludes categories not present in every strain

all_data %<>% bind_rows(data)

## Niccum 2018 (DNA pol proofreading knockouts)
codes = c('D5_WT', rep('D5_mutL', 3), 'D5_dinB', rep('D5_dinB_umuDC', 2), 'D5_lexA3', 'D5_WT-MM', 'D5_mutL-MM')
names(codes) = as.character(c(163, 165, 397, 399, 479, 515, 517, 686, '163m', '165m'))

data = read_xlsx('mut_accumulation_data/Niccum 2018.xlsx') %>%
  mutate(genome = 'mg1655', source = 'niccum_2018') %>%
  # filter(consensus == ref) %>%
  separate(line, into = c('code', 'line'), sep = '-') %>%
  mutate(strain = codes[code]) %>%
  separate(strain, into = c('strain', 'condition'), sep = '-') %>%
  mutate(condition = ifelse(is.na(condition), 'LB', condition))

all_data %<>% bind_rows(data)


# ## Sane 2023
# data = read_xlsx('mut_accumulation_data/Sane 2023.xlsx') %>%
#   mutate(genome = 'long_2016', source = 'sane_2023') %>%
#   filter(strain == 'wt') %>%
#   mutate(line = as.character(line)) %>%
#   rename(pos = genomic_posn,
#          ref = orig_base,
#          alt = mut_base) %>%
#   mutate(condition = 'LB',
#          code = 'sane_2023_wt',
#          strain = 'WT') %>%
#   select(source, genome, code, strain, condition, line, pos, ref, alt)
# 
# all_data %<>% bind_rows(data)


## Tincher 2017
data = read_xlsx('mut_accumulation_data/Tincher 2017.xlsx', skip = 1) %>%
  mutate(genome = 'long_2016', source = 'tincher_2017') %>%
  mutate(samples = sub('Sample_', '', samples),
         strain = case_when(
           substr(samples, 1, 1) == 'A' ~ 'WT',
           substr(samples, 1, 2) == 'SA' ~ 'mutS',
           TRUE ~ NA_character_
         )) %>%
  filter(!is.na(strain)) %>%
  separate(substitution, into = c('ref', 'alt'), sep = '>') %>%
  mutate(pos = sub('NC_000913:', '', positions) %>% as.numeric()) %>%
  rename(line = samples) %>%
  mutate(condition = 'LB',
         code = ifelse(strain == 'WT', 'tincher_2017_wt', 'tincher_2017_mutS')) %>%
  select(source, genome, code, strain, condition, line, pos, ref, alt)

all_data %<>% bind_rows(data)



# classifies mutations into different BPS types
all_data %<>%
  mutate(snp = paste(ref, alt, sep = '')) %>%
  mutate(snp = case_when(
    snp %in% c('AG', 'TC') ~ spectrum[1],
    snp %in% c('GA', 'CT') ~ spectrum[2],
    snp %in% c('AC', 'TG') ~ spectrum[3],
    snp %in% c('GT', 'CA') ~ spectrum[4],
    snp %in% c('AT', 'TA') ~ spectrum[5],
    snp %in% c('GC', 'CG') ~ spectrum[6],
    TRUE ~ NA_character_
  ))

all_data$snp %<>% factor(levels = spectrum)
all_data %<>% filter(!is.na(snp)) # excludes insertions and deletions


## determines replichore
# Coordinates for OriC and Ter in E. coli K-12

OriC <- 3925696
Ter <- 1640202

# Annotate mutations with replichore
all_data %<>%
  rowwise() %>%
  mutate(
    # Replichore assignment
    replichore = case_when(
      !is.na(replichore) ~ tolower(replichore),
      TRUE ~ case_when(
        pos >= Ter & pos < OriC ~ 'left',
        TRUE ~ 'right'
      )
    )
  ) %>%
  ungroup()


## updates NC_000913.2 (mg1655) positions to NC_000913.3 (long_2016)
all_data %<>%
  mutate(pos = case_when(
    genome == 'mg1655' ~ case_when(
          pos <= 257901 ~ pos,
          pos >= 257902 & pos <= 547831 ~ pos + 776,
          pos >= 547832 & pos <= 1298719 ~ pos + 776 + 1,
          pos >= 1298720 & pos <= 2171384 ~ pos + 776 + 1 + 1199,
          pos >= 2171385 & pos <= 3558477 ~ pos + 776 + 1 + 1199 + 2,
          pos >= 3558479 ~ pos + 776 + 1 + 1199 + 2 - 1 
  ), TRUE ~ pos))



all_data %<>% select(-consensus)


# checks which experiments have mismatches between position and reference sequence
filter(all_data, ref != mg1655[pos]) %>%
  group_by(source, strain, code) %>%
  summarise(count = n())


#### makes the curated dataset ####
my_data = filter(all_data, condition == 'LB' & !grepl('subtilis', genome)
            & genome %in% c('mg1655', 'long_2016')) %>%
  mutate(group = case_when(
    strain %in% c('WT', 'ED1a', 'IAI1') ~ groups[1],
    strain %in% c('uvrA', 'alkA_tagA', 'ada_ogt', 'nfi', 'umuDC_dinB', 'umuDC_dinB_polB') ~ groups[1],
    strain %in% c('mutS', 'mutL', 'mutH', 'mutS_mutL', 'mutS_mutL_mutH') ~ groups[2],
    strain %in% c('mutL_umuDC_dinB', 'mutS_mfd', 'mutL_mfd') ~ groups[2],
    strain %in% c('D5_WT', 'D5_dinB', 'D5_dinB_umuDC') ~ groups[3],
    strain %in% c('D5_mutL') ~ groups[4],
    TRUE ~ NA
  )) %>%
  mutate(group = factor(group, levels = groups)) %>%
  filter(!is.na(group)) %>%
  select(-genome, -condition)

# identifies the -10 to +10 motif, used for later analyses
my_data %<>%
  rowwise() %>%
  mutate(motif = case_when(
    ref %in% c('A', 'G') ~ mg1655[(pos-10):(pos+10)] %>% paste(collapse = ''),
    TRUE ~ rev.comp(mg1655[(pos-10):(pos+10)]) %>% paste(collapse = '')
  )) %>%
  ungroup()

# determines if purine is on LDST or LGST, used for later analyses
my_data %<>%
  mutate(strand = case_when(
    (replichore == 'left' & ref %in% c('A', 'G')) | (replichore == 'right' & ref %in% c('T', 'C')) ~ 'LDST',
    TRUE ~ 'LGST'
  ))



# creates wrapped genome for motif analyses
wrap.genome = as.character(mg1655) %>% paste0(collapse = '')
wrap.genome = paste0(substr(wrap.genome, gen.length - 9, gen.length), wrap.genome, substr(wrap.genome, 1, 10))

# creates genome motif table for later analyses
genome.motifs = tibble(
  pos = 1:length(mg1655),
  ref = as.character(mg1655)
) %>%
  mutate(motif = substring(wrap.genome, pos, pos + 20)) %>%
  mutate(motif = ifelse(ref %in% c('T', 'C'), Vectorize(rev.comp)(motif), motif)) %>%
  mutate(focus = substr(motif, 11, 11))



#### Table 1 (collated data and mutational spectrum) ####
generations = read_csv('mut_accumulation_data/summary.csv') %>%
  filter(strain %in% unique(my_data$strain))

rates_groups <- my_data %>%
  group_by(group, strain, code, snp) %>%
  summarise(count = n(), .groups = "drop") %>%
  left_join(generations %>% select(code, generations)) %>%
  group_by(group, snp) %>%
  summarise(
    count = sum(count),
    generations = sum(generations),
    mean_rate = count / generations / gen.length) %>%
  mutate(
    lower = 0.5 * qchisq(0.025, 2 * count) / generations / gen.length,
    upper = 0.5 * qchisq(0.975, 2 * (count + 1)) / generations / gen.length,
  )


rates_summary = generations %>%
  group_by(group) %>%
  summarise(mean_rate = sum(BPS) / sum(generations) / gen.length,
            generations = sum(generations),
            lower = 0.5 * qchisq(0.025, 2 * sum(BPS)) / sum(generations) / gen.length,
            upper = 0.5 * qchisq(0.975, 2 * (sum(BPS) + 1)) / sum(generations) / gen.length,
            .groups = "drop")


breaks = list(c(3e-11, 6e-11, 9e-11), c(1e-8, 2e-8, 3e-8), c(3e-7, 6e-7, 9e-7),  c(3e-7, 6e-7, 9e-7))

limits = tibble(group = groups,
                lim = c(9e-11, 3e-8, 9e-7, 9e-7))

plot_list = list() ; for (i in 1:4) {
  plot_data = filter(rates_groups, group == groups[i])
  limits2 = filter(limits, group == groups[i])
  
  plot_list %<>% append(list(
    ggplot(plot_data) +
      geom_bar(aes(x = snp, y = mean_rate, fill = snp),
               stat = 'identity', position = 'dodge') +
      geom_errorbar(aes(x = snp, ymax = upper,
                        ymin = lower),
                    width = 0.6) +
      # geom_label(data = stats, aes(x = snp, y = ymax + (0.1 * max), label = ratio), size = 2,
      #            label.size = 0, label.padding = unit(0, 'points')) +
      geom_hline(data = limits2, aes(yintercept = lim), linewidth = 0) +
      scale_x_discrete(labels = function(x) sub('>', '→', x)) +
      scale_y_continuous(labels = function(x) sci_label(x, digits = 0),
                         breaks = breaks[[i]]) +
      scale_fill_manual(values = pal_spec, guide = NULL) +
      facet_wrap(~group, scales = 'free_y', nrow = 1) +
      labs(x = NULL, y = 'Mutation rate<br>(per bp per generation)', fill = 'Base pair sites where purine templates the:') +
      theme_minimal() +
      theme(panel.grid.major.x = element_blank(),
            strip.background = element_rect(fill = NA, color = NA),
            strip.text = element_blank(), axis.title.y = element_blank(),
            plot.background = element_blank(),
            axis.text.y.left = element_markdown(),
            axis.text.x = element_markdown(angle = 90, vjust = 0.3,
                                           margin = margin(t = 2),
                                           family = 'Aptos Mono',
                                           size = 9))
  ))
}

avg_rate_labels = rates_summary$mean_rate %>% sci_label(digits = 1)

avg_rate_labels = paste0(substr(avg_rate_labels, 1, 3),
       ' (± 0.0', substr(as.character(rates_summary$upper - rates_summary$mean_rate), 1, 1), ')<br>',
       substr(avg_rate_labels, 4, 20))

table = tibble(x = c(1, 2, 3, 4) %>% rep(5),
               y = c(6, 5.2, 4.4, 3.8, 3.1) %>% rep(each = 4),
               text = c(groups %>% gsub(' ', '<br>', .) %>% paste0('**', ., '**'),
                        c('PFM2 (WT),<br>*ΔuvrA, Δnfi,<br>Δada+ogt<br>ΔalkA+tagA<br>ΔdinB+umuDC<br>ΔdinB+umuDC+polB*',
                          '*ΔmutS, ΔmutL,<br>ΔmutH, ΔmutSL<br>ΔmutSLH,<br>ΔmutS+mfd,<br>ΔmutL+mfd,<br>ΔmutL+dinB+umuDC*',
                          '*dnaQ*-T15I,<br>*dnaQ*-T15I (*ΔdinB* ),<br>*dnaQ*-T15I<br>(*ΔdinB+umuDC*)',
                          '*dnaQ*-T15I (*ΔmutL*)'),
                        my_data$group %>% table() %>% unname() %>% prettyNum(big.mark = ','), # no. mutations
                        rates_summary$generations %>% prettyNum(big.mark = ','), # no. generations
                        avg_rate_labels # mean mutation rate
               ))

headings = tibble(x = 0, y = c(unique(table$y), 1.7),
                  text = c('DNA repair<br>group',
                           'Strains<br>used',
                           'Total #<br>of BPS',
                           'Total # of<br>generations',
                           'Avg. BPS rate<br><sub>(per base pair<br>per generation)</sub>',
                           'Mutational<br>spectrum'))

a = ggplot() +
  geom_richtext(data = table, aes(x = x, y = y, label = text, size = as.character(y)),
                label.size = 0, fill = NA) +
  scale_size_manual(values = c(5, 5, 5, 3.7, 4), guide = NULL)+
  geom_richtext(data = headings, aes(x = x, y = y, label = text),
                fontface = 'bold', size = 5, label.size = 0, fill = NA) +
  geom_hline(yintercept = c(2.7, 3.5, 4.1, 4.7, 5.7), color = 'grey') +
  geom_vline(xintercept = c(0.5, 1.5, 2.5, 3.5), color = 'grey') +
  coord_cartesian(xlim = c(-0.2, 4.3), ylim = c(1, 6)) +
  theme_void()


ggdraw() +
  draw_plot(a, x = 0, y = 0, width = 1, height = 1) +
  # draw_plot(b, x = 0.2, y = 0, width = 0.8, height = 0.3) +
  draw_plot(plot_list[[1]], x = 0.19, y = 0, height = 0.35, width = 0.2) +
  draw_plot(plot_list[[2]], x = 0.39, y = 0, height = 0.35, width = 0.2) +
  draw_plot(plot_list[[3]], x = 0.59, y = 0, height = 0.35, width = 0.2) +
  draw_plot(plot_list[[4]], x = 0.794, y = 0, height = 0.35, width = 0.2)


ggsave('Table 1.svg', width = 9, height = 7)



## stats
rates_groups %>%
  select(group, snp, count) %>%
  pivot_wider(names_from = snp, values_from = count) %>%
  column_to_rownames("group") %>%
  as.matrix() %>%
  chisq.test()


plot_data = rates_groups %>%
  filter(group != groups[1]) %>%
  left_join(rates_groups %>%
              filter(group == groups[1]) %>%
              rename(base_rate = mean_rate) %>%
              select(snp, base_rate)) %>%
  mutate(fold = mean_rate / base_rate)


ggplot(plot_data, aes(x = snp, y = fold, fill = snp)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = pal_spec) +
  labs(x = NULL, y = 'Fold increase in mutation rate') +
  facet_wrap(~group, scales = 'free_y') +
  theme_bw()






#### Figure S1 (Mutational spectrum of strains in each repair group) ####

plot_data <- my_data %>%
  group_by(group, strain, code, snp) %>%
  summarise(count = n(), .groups = "drop") %>%
  left_join(generations %>% select(code, generations)) %>%
  group_by(group, strain, snp) %>%
  summarise(
    count = sum(count),
    generations = sum(generations),
    mean_rate = count / generations / gen.length) %>%
  mutate(
    lower = 0.5 * qchisq(0.025, 2 * count) / generations / gen.length,
    upper = 0.5 * qchisq(0.975, 2 * (count + 1)) / generations / gen.length,
  )


plot_data %<>%
  mutate(strain = case_when(
    strain == 'WT' ~ 'WT',
    strain == 'D5_WT' ~ 'dnaQ-T15I',
    strain == 'D5_dinB' ~ 'dnaQ-T15I (ΔdinB)',
    strain == 'D5_dinB_umuDC' ~ 'dnaQ-T15I (ΔdinB+umuDC)',
    TRUE ~ paste0('Δ', strain)
  )) %>%
  mutate(strain = sub('_', '+', strain),
         strain = sub('_', '\n+', strain))

plot_data$strain %<>% factor(
  levels = c('WT', 'ΔuvrA', 'Δnfi', 'Δada+ogt', 'ΔalkA+tagA', 'ΔumuDC+dinB', 'ΔumuDC+dinB\n+polB',
             'ΔmutS', 'ΔmutL', 'ΔmutH', 'ΔmutS+mutL', 'ΔmutS+mutL\n+mutH', 'ΔmutL+umuDC\n+dinB', 'ΔmutS+mfd', 'ΔmutL+mfd',
             'dnaQ-T15I', 'dnaQ-T15I (ΔdinB)', 'dnaQ-T15I (ΔdinB+umuDC)', 'dnaQ-T15I (ΔmutL)'))


ggplot(plot_data %>% filter(group != groups[4]), aes(x = strain, y = mean_rate, fill = snp)) +
  geom_bar(stat = 'identity', position = position_dodge(0.9)) +
  geom_errorbar(aes(ymax = upper,
                    ymin = lower),
                position = position_dodge(0.9), width = 0.5) +
  scale_x_discrete(labels = function(x) sub('>', '→', x)) +
  scale_y_continuous(labels = function(x) sci_label(x, digits = 1)) +
  scale_fill_manual(values = pal_spec, guide = NULL) +
  facet_wrap(~group, scales = 'free', ncol = 1) +
  labs(x = NULL, y = 'Mutation rate (per bp per generation)', fill = 'Base pair sites where purine templates the:') +
  theme_bw() +
  theme(axis.text.y = element_markdown(),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_markdown(),
        axis.ticks = element_blank(),
        strip.text = element_text(size = 12),
        strip.background = element_blank())

ggsave('figS1.svg', width = 8, height = 8)


#### Table S1 ####


strains = c('PFM2 (WT)', 'ΔuvrA', 'Δnfi', 'Δada+ogt', 'ΔalkA+tagA', 'ΔdinB+umuDC', 'ΔdinB+umuDC+polB',
            'ΔmutS', 'ΔmutL', 'ΔmutH', 'ΔmutSL', 'ΔmutSLH', 'ΔmutS+mfd', 'ΔmutL+mfd', 'ΔmutL+dinB+umuDC',
            'dnaQ-T15I', 'dnaQ-T15I (ΔdinB)', 'dnaQ-T15I (ΔdinB+umuDC)',
            'dnaQ-T15I (ΔmutL)')

names(strains) = c('WT', 'uvrA', 'nfi', 'ada_ogt', 'alkA_tagA', 'umuDC_dinB', 'umuDC_dinB_polB',
                   'mutS', 'mutL', 'mutH', 'mutS_mutL', 'mutS_mutL_mutH', 'mutS_mfd', 'mutL_mfd', 'mutL_umuDC_dinB',
                   'D5_WT', 'D5_dinB', 'D5_dinB_umuDC',
                   'D5_mutL')

my_data %>%
  group_by(group, source, strain, code) %>%
  summarise(BPS = n(), lines = length(unique(line)), .groups = "drop") %>%
  left_join(generations %>% select(code, generations)) %>%
  group_by(group, source, strain) %>%
  mutate(experiment = dense_rank(code),
         strain = strains[strain],
         strain = factor(strain, levels = strains),
         source = factor(source, levels = c('lee_2012', 'foster_2015', 'long_2016', 'tincher_2017', 'foster_2018', 'niccum_2018'))) %>%
  ungroup() %>%
  select(group, source, strain, experiment, lines, BPS, generations, -code) %>%
  arrange(group, source, strain, experiment, lines) %>%
  mutate(rate = BPS / generations / gen.length,
         lower = 0.5 * qchisq(0.025, 2 * BPS) / generations / gen.length,
         upper = 0.5 * qchisq(0.975, 2 * (BPS + 1)) / generations / gen.length) %>%
  write_csv('Table S1.csv')





#### Figure S2 (Average genomic sequence context) ####


genome_context = genome.motifs %>%
  mutate(base = list(5:17)) %>%
  unnest(base) %>%
  mutate(value = substr(motif, base, base)) %>%
  group_by(focus, base, value) %>%
  summarise(genome.count = n(), .groups = "drop")

genome_context %<>%
  left_join(genome_context %>%
              group_by(focus, base) %>%
              summarise(total = sum(genome.count))) %>%
  mutate(perc = genome.count / total,
         base = base - 11,
         value = factor(value, levels = c('A', 'T', 'G', 'C'))) %>%
  select(-total)


ggplot(genome_context, aes(x = base, y = perc, fill = value)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = pal_bases, name = 'Nucleotide') +
  scale_y_continuous(breaks = c(0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = seq(-6, 6, 2), labels = function(x) ifelse(x > 0, paste0('+', x), x)) +
  geom_hline(yintercept = c(0.5), linewidth = 0.5) +
  geom_hline(yintercept = c(0.25, 0.75), linewidth = 0.25) +
  facet_wrap(~focus, labeller = as_labeller(function(x) paste0('Context around every ', x)), ncol = 2) +
  coord_cartesian(xlim = c(-6, 6)) +
  labs(y = 'Frequency', x = "Context Position", fill = "Value") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  guides(fill = guide_legend(position = 'top', order = 1))


ggsave('figS2.svg', width = 6, height = 3)

#### Figure 1 (sequence context nucleotide frequencies) ####

bases = c(1:6)
bases = c(rev(bases)*-1, 0, bases)

ctx_freq = my_data %>%
  mutate(base = list(bases)) %>%
  unnest(base) %>%
  mutate(value = substr(motif, 11 + base, 11 + base))

ctx_freq$value %<>% factor(levels = c('A', 'T', 'G', 'C'))

plot_data <- ctx_freq %>%
  group_by(snp, base, value, group) %>%
  summarise(count = n(), .groups = "drop") %>%
  complete(snp, base, value, group, fill = list(n = 0)) # Fill missing combinations

plot_data %<>%
  mutate(focus = case_when(grepl('AT>', snp) ~ 'A', TRUE ~ 'G')) %>% # adds col for focus base (A or G)
  left_join(genome_context) %>% # adds info for the genome_context bias at each position
  mutate(weight.count = (0.25 / perc) * count) # weights by genome structure

plot_data %<>%
  left_join(
    plot_data %>%
      group_by(snp, base, group) %>%
      summarise(total = sum(weight.count))
  ) %>%
  mutate(perc = weight.count/total,
         perc.sem = sqrt(perc * (1 - perc) / total)) %>%
  ungroup() %>%
  select(-total)


# chi squared test to determine if significant bias
stats = plot_data %>%
  group_by(snp, base, focus, group) %>%
  summarise(n = sum(count))
  
p = tibble()
for (i in 1:nrow(stats)) {
  if (stats[[i,'base']] == 0) {p %<>% bind_rows(c(chi = NA, p = NA)) ; next}
  x = semi_join(ctx_freq, stats[i,], by = join_by(snp, base, group))$value %>% table()
  
  expected = semi_join(genome_context, stats[i,], by = join_by(base, focus)) %>%
    with(setNames(perc, value))
  
  y = chisq.test(x, p = expected)[c('statistic', 'p.value')] %>% unlist()
  y %<>% c(sum(abs(x/sum(x) - expected)))
  names(y) = c('chi', 'p', 'bias')
  
  p %<>% bind_rows(y)
}

p %<>% mutate(p.adj = p.adjust(p, method = 'fdr')) 


stats %<>% bind_cols(p)


plot_data %<>%
  left_join(stats)


# stops base 0 from skewing results with maximum bias
plot_data$bias[plot_data$base == 0] = 
  max(plot_data$bias[plot_data$base != 0])
plot_data$perc[plot_data$base == 0 & !is.na(plot_data$count)] = 1


a = ggplot(plot_data, aes(x = base, y = perc, fill = value, alpha = bias)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = pal_bases, name = 'Nucleotide') +
  scale_alpha_continuous(
    range = c(0, 1),
    name = 'Bias',
    limits = c(0, 0.5),
    breaks = c(0, 0.25, 0.5),
    labels = c('0', '0.25', '> 0.5')) +
  scale_y_continuous(breaks = c(0.25, 0.5, 0.75), sec.axis = dup_axis(name = 'Weighted frequency')) +
  scale_x_continuous(breaks = seq(-6, 6, 2), labels = c('-6', '-4', '-2', '0', '+2', '+4', '+6')) +
  geom_hline(yintercept = c(0.5), linewidth = 0.5) +
  geom_hline(yintercept = c(0.25, 0.75), linewidth = 0.25) +
  # facet_wrap(~snp, labeller = as_labeller(labels)) +
  facet_grid(rows = vars(snp), cols = vars(group), switch = 'y',
             labeller = as_labeller(function(x) sub(">", "→", x))
  ) +
  coord_cartesian(xlim = c(-6, 6)) +
  labs(y = NULL, x = "Context Position", fill = "Value") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y.left = element_text(family = 'Aptos'), 
        axis.text.y.left = element_blank(),
        legend.text = element_text(size = 10)) + 
  guides(fill = guide_legend(position = 'top', order = 1),
         alpha = guide_legend(position = 'top', order = 2,
                              override.aes = list(color = 'black',
                                                  linewidth = 0.1,
                                                  fill = '#B6C999')))



## plot for bias at each position, first estimates min significant bias
counts = my_data %>%
  group_by(group, snp) %>%
  summarise(count = n(), .groups = 'drop') %>%
  arrange(group, snp)

# splits a number (n) into four random numbers that sum to n, normally distributed around 0.25 split
split_four <- function(n, sd = 0.07) {
  if (n == 0) return(rep(0, 4))
  
  probs <- rnorm(4, mean = 0.25, sd = sd)
  probs[probs < 0] <- 0  # truncate negatives
  probs <- probs / sum(probs)  # renormalize to sum to 1
  
  counts <- floor(probs * n)
  remainder <- n - sum(counts)
  
  # randomly distribute the remainder to make total = n
  if (remainder > 0) {
    add <- sample(1:4, remainder, replace = TRUE)
    for (i in add) counts[i] <- counts[i] + 1
  }
  
  return(counts)
}

many_stats = tibble()
for (i in 1:100) {
  stats = plot_data %>% distinct(snp, base, focus, group) %>%
    filter(base != 0)
  
  p = tibble()
  for (i in 1:nrow(stats)) {
    x = semi_join(counts, stats[i,], by = join_by(snp, group))$count %>% split_four()
    
    expected = semi_join(genome_context, stats[i,], by = join_by(base, focus)) %>%
      with(setNames(perc, value))
    
    p %<>% bind_rows(c(bias = abs((x / sum(x)) - expected) %>% sum(),
                       p = chisq.test(x, p = expected)$p.value))
  }
  
  p$p %<>% p.adjust(method = 'fdr')
  stats %<>% bind_cols(p = p)
  
  many_stats %<>% bind_rows(stats)
}

ggplot(many_stats) +
  geom_hline(yintercept = 0.01) +
  geom_point(aes(x = bias, y = p, color = snp)) +
  scale_y_log10() +
  coord_cartesian(ylim = c(1e-100, 1)) +
  facet_grid(rows = vars(snp), cols = vars(group))

many_stats %<>% filter(p < 0.05) %>%
  group_by(group, snp) %>%
  summarise(thresh = min(bias), .groups = 'drop')

b = ggplot(filter(plot_data, base != 0), aes(x = base, y = bias, color = snp)) +
  geom_hline(data = many_stats, aes(yintercept = thresh), linetype = 'dashed', alpha = 0.5) +
  geom_vline(xintercept = 0, color = 'grey') +
  geom_point(aes(alpha = p <= 0.05)) +
  geom_label(data = counts, aes(label = paste0('n=', count)),
             x = -6.4 , y = 1, color = 'black', fill = 'white', size = 3,
             label.size = 0, hjust = 0) +
  # geom_label(data = counts, aes(label = loss),
  # x = 6.4 , y = 1, color = 'black', fill = 'white', size = 3,
  # label.size = 0, hjust = 1) +
  facet_grid(rows = vars(snp), cols = vars(group), switch = 'y',
             labeller = as_labeller(function(x) sub(">", "→", x))) +
  scale_x_continuous(breaks = seq(-6, 6, 2), labels = c('-6', '-4', '-2', '0', '+2', '+4', '+6')) +
  scale_y_continuous(sec.axis = dup_axis(name = 'Absolute difference from expected frequencies'),
                     breaks = c(0, 0.5, 1)) +
  scale_color_manual(values = pal_spec, guide = NULL) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.1), guide = NULL) + 
  coord_cartesian(ylim = c(0, 1.2)) +
  labs(y = NULL, x = "Context Position", color = "SNP") +
  theme_bw() +
  theme(axis.text.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        strip.text.y.left = element_text(family = 'Aptos'),
        strip.background = element_rect(fill = NA, color = NA))


cowplot::plot_grid(plotlist = list(a, b), ncol = 1, labels = c('B', 'C'))
ggsave('fig1bc.svg', width = 7, height = 9)


#### Figure S3 (Dam and Dcm hotspots) ####
plot_data = my_data %>%
  mutate(type = case_when(
    substr(motif, 10, 13) == 'GATC' ~ 'G**A**TC',
    substr(motif, 8, 12) %in% c('CCTGG', 'CCAGG') ~ 'C**C**WGG',
    TRUE ~ 'other'
  ))


plot_data %<>%
  group_by(group, snp, type) %>%
  summarise(count = n()) %>%
  left_join(plot_data %>%
              group_by(group, snp) %>%
              summarise(total = n())) %>%
  mutate(prop = count/total,
         type = factor(type, levels = c('G**A**TC', 'C**C**WGG', 'other')))


expected = tibble(
  group = rep('Whole\ngenome', 3),
  type = c('G**A**TC', 'C**C**WGG', 'other'),
  gen.count = c(38240, 24090, gen.length - 24090 - 38240)
) %>% mutate(prop = gen.count / gen.length,
             type = factor(type, levels = c('G**A**TC', 'C**C**WGG', 'other')))



p1 = plot_data %>%
  mutate(type = factor(type, levels = c('G**A**TC', 'C**C**WGG', 'other'))) %>%
  ggplot(aes(x = snp, y = prop, fill = type)) +
  geom_bar(stat = 'identity') +
  scale_x_discrete(labels = function(x) sub('>', '→', x)) +
  scale_fill_manual(values = c('#DBA237', '#469C76', 'lightgrey'), guide = NULL) +
  facet_wrap(~group, ncol = 1) +
  labs(x = NULL, y = 'Proportion of all mutations', fill = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(family = 'Aptos'))

p2 = expected %>%
  mutate(type = factor(type, levels = c('G**A**TC', 'C**C**WGG', 'other'))) %>%
  ggplot(aes(x = '', y = prop, fill = type)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c('#DBA237', '#469C76', 'lightgrey'), guide = NULL) +
  facet_wrap(~group, ncol = 1) +
  labs(x = NULL, y = 'Proportion of sites in genome', fill = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_blank())


plot_data %<>%
  left_join(expected %>% select(-prop, -group)) %>%
  left_join(rates_summary %>% select(group, generations)) %>%
  mutate(rate = count / gen.count / generations)


p3 = ggplot(plot_data %>% mutate(group = factor(group, levels = groups)),
            aes(x = snp, y = rate, fill = type)) +
  geom_bar(stat = 'identity', position = 'dodge', color = 'black', linewidth = 0.3) +
  scale_fill_manual(values = c('#DBA237', '#469C76', 'lightgrey')) +
  scale_x_discrete(labels = function(x) sub('>', '→', x)) +
  scale_y_continuous(labels = function(x ) sci_label(x, digits = 1)) +
  facet_wrap(~group, ncol = 1, scales = 'free_y') +
  labs(x = NULL, y = 'Mutation rate (per site in genome per generation)', fill = NULL) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_markdown(),
        axis.text.x = element_markdown(family = 'Aptos'),
        legend.text = element_markdown()) +
  guides(fill = guide_legend(position = 'top'))

plot_grid(p3, p2, p1, rel_widths = c(1, 0.5, 1), ncol = 3, labels = c('A', 'B'))

ggsave('figS3.png', width = 8.5, height = 6, dpi = 1000)



#### Figure S4 (Run of 3’ C hotspots) ####

c.runs.genome = genome.motifs %>%
  mutate(run.length = attr(regexpr("^C*", substring(motif, 12, 21), perl = TRUE), "match.length")) %>%
  group_by(run.length, focus) %>%
  summarise(gen.count = n(), .groups = 'drop')


c.runs.genome %<>%
  left_join(c.runs.genome %>%
              group_by(focus) %>%
              summarise(total = sum(gen.count))) %>%
  mutate(perc.gen = gen.count / total)


plot_data = my_data %>%
  mutate(run.length = attr(regexpr("^C*", substring(motif, 12, 21), perl = TRUE), "match.length")) %>%
  group_by(group, snp, run.length) %>%
  summarise(count = n(), .groups = 'drop')


plot_data %<>%
  left_join(plot_data %>%
              group_by(group, snp) %>%
              summarise(total = sum(count))) %>%
  mutate(perc = count / total,
         focus = substr(snp, 1, 1))



# join expected proportions
plot_data %<>%
  left_join(c.runs.genome %>% select(focus, run.length, perc.gen)) %>%
  rowwise() %>%
  mutate(pval = binom.test(count, total, p = perc.gen)$p.value) %>%
  ungroup() %>%
  mutate(pval = p.adjust(pval, method = 'fdr')) %>%
  mutate(sig = case_when(
    pval < 0.001 ~ "***",
    pval < 0.01  ~ "**",
    pval < 0.05  ~ "*",
    TRUE         ~ ""
  ))


plot_list = list()
for (mut in c('AT>CG', 'GC>CG')) {
  plot_list %<>% append(list(
    ggplot(c.runs.genome %>% filter(focus == substr(mut, 1, 1)),
           aes(x = run.length, y = perc.gen)) +
      geom_point(color = 'red') + geom_line(color = 'red') +
      geom_point(data = plot_data %>% filter(snp == mut),
                 aes(x = run.length, y = perc), color = 'black') +
      geom_text(data = plot_data %>% filter(snp == mut),
                aes(x = run.length, y = perc * 2.3, label = sig),  # adjust multiplier to position stars
                size = 3.5, vjust = 0) +
      scale_y_log10(labels = c(1, '10<sup>-2</sup>', '10<sup>-4</sup>', '10<sup>-6</sup>'),
                    breaks = c(1, 1e-2, 1e-4, 1e-6)) +
      scale_x_continuous(breaks = 0:10) +
      coord_cartesian(xlim = c(0, 10), ylim = c(1e-7, 10)) +
      labs(x = 'Run length of Cs 3\' of mutation site',
           y = paste0('Proportion of ', sub('>', '→', mut), ' mutations')) +
      facet_wrap(~group, ncol = 1, scales = 'free_x') +
      theme_bw() +
      theme(strip.background = element_blank(),
            axis.text.y = element_markdown(),
            axis.title = element_markdown(family = 'Aptos'),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank())
  ))
}

plot_grid(plotlist = plot_list, ncol = 2, labels = 'AUTO')

ggsave('figS4.png', height = 7, width = 7, , dpi = 1000)



#### Figure S5 (Fig 1 without hotspots) ####

filtered_data = my_data %>%
  mutate(run.length = attr(regexpr("^C*", substring(motif, 12, 21), perl = TRUE), "match.length"),
         hotspot = case_when(
           substr(motif, 10, 13) == 'GATC' ~ 'GATC',
           substr(motif, 8, 12) %in% c('CCTGG', 'CCAGG') ~ 'CCWGG',
           run.length >= 3 & snp %in% c('AT>CG', 'GC>CG') ~ 'C3+',
           TRUE ~ NA_character_
         )) %>%
  filter(is.na(hotspot))


bases = c(1:6)
bases = c(rev(bases)*-1, 0, bases)

ctx_freq = filtered_data %>%
  mutate(base = list(bases)) %>%
  unnest(base) %>%
  mutate(value = substr(motif, 11 + base, 11 + base))

ctx_freq$value %<>% factor(levels = c('A', 'T', 'G', 'C'))

plot_data <- ctx_freq %>%
  group_by(snp, base, value, group) %>%
  summarise(count = n(), .groups = "drop") %>%
  complete(snp, base, value, group, fill = list(n = 0)) # Fill missing combinations

plot_data %<>%
  mutate(focus = case_when(grepl('AT>', snp) ~ 'A', TRUE ~ 'G')) %>% # adds col for focus base (A or G)
  left_join(genome_context) %>% # adds info for the genome_context bias at each position
  mutate(weight.count = (0.25 / perc) * count) # weights by genome structure

plot_data %<>%
  left_join(
    plot_data %>%
      group_by(snp, base, group) %>%
      summarise(total = sum(weight.count))
  ) %>%
  mutate(perc = weight.count/total,
         perc.sem = sqrt(perc * (1 - perc) / total)) %>%
  ungroup() %>%
  select(-total)


# chi squared test to determine if significant bias
stats = plot_data %>%
  group_by(snp, base, focus, group) %>%
  summarise(n = sum(count))

p = tibble()
for (i in 1:nrow(stats)) {
  if (stats[[i,'base']] == 0) {p %<>% bind_rows(c(chi = NA, p = NA)) ; next}
  x = semi_join(ctx_freq, stats[i,], by = join_by(snp, base, group))$value %>% table()
  
  expected = semi_join(genome_context, stats[i,], by = join_by(base, focus)) %>%
    with(setNames(perc, value))
  
  y = chisq.test(x, p = expected)[c('statistic', 'p.value')] %>% unlist()
  y %<>% c(sum(abs(x/sum(x) - expected)))
  names(y) = c('chi', 'p', 'bias')
  
  p %<>% bind_rows(y)
}

p %<>% mutate(p.adj = p.adjust(p, method = 'fdr')) 


stats %<>% bind_cols(p)


plot_data %<>%
  left_join(stats)


# stops base 0 from skewing results with maximum bias
plot_data$bias[plot_data$base == 0] = 
  max(plot_data$bias[plot_data$base != 0])
plot_data$perc[plot_data$base == 0 & !is.na(plot_data$count)] = 1


a = ggplot(plot_data, aes(x = base, y = perc, fill = value, alpha = bias)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = pal_bases, name = 'Nucleotide') +
  scale_alpha_continuous(
    range = c(0, 1),
    name = 'Bias',
    limits = c(0, 0.5),
    breaks = c(0, 0.25, 0.5),
    labels = c('0', '0.25', '> 0.5')) +
  scale_y_continuous(breaks = c(0.25, 0.5, 0.75), sec.axis = dup_axis(name = 'Weighted frequency')) +
  scale_x_continuous(breaks = seq(-6, 6, 2), labels = c('-6', '-4', '-2', '0', '+2', '+4', '+6')) +
  geom_hline(yintercept = c(0.5), linewidth = 0.5) +
  geom_hline(yintercept = c(0.25, 0.75), linewidth = 0.25) +
  # facet_wrap(~snp, labeller = as_labeller(labels)) +
  facet_grid(rows = vars(snp), cols = vars(group), switch = 'y',
             labeller = as_labeller(function(x) sub(">", "→", x))
  ) +
  coord_cartesian(xlim = c(-6, 6)) +
  labs(y = NULL, x = "Context Position", fill = "Value") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y.left = element_text(family = 'Aptos'), 
        axis.text.y.left = element_blank(),
        legend.text = element_text(size = 10)) + 
  guides(fill = guide_legend(position = 'top', order = 1),
         alpha = guide_legend(position = 'top', order = 2,
                              override.aes = list(color = 'black',
                                                  linewidth = 0.1,
                                                  fill = '#B6C999')))



## plot for bias at each position, first estimates min significant bias
og.counts = my_data %>%
  group_by(group, snp) %>%
  summarise(og.count = n(), .groups = 'drop')

filt.counts = filtered_data %>%
  group_by(group, snp) %>%
  summarise(count = n(), .groups = 'drop') %>%
  arrange(group, snp)

filt.counts %<>%
  left_join(og.counts) %>%
  mutate(difference = og.count - count) %>%
  mutate(loss = round(difference / og.count, 2) * 100) %>%
  mutate(difference = paste0('-', difference)) %>%
  mutate(loss = paste0('-', loss, '%'))

# splits a number (n) into four random numbers that sum to n, normally distributed around 0.25 split
split_four <- function(n, sd = 0.07) {
  if (n == 0) return(rep(0, 4))
  
  probs <- rnorm(4, mean = 0.25, sd = sd)
  probs[probs < 0] <- 0  # truncate negatives
  probs <- probs / sum(probs)  # renormalize to sum to 1
  
  counts <- floor(probs * n)
  remainder <- n - sum(counts)
  
  # randomly distribute the remainder to make total = n
  if (remainder > 0) {
    add <- sample(1:4, remainder, replace = TRUE)
    for (i in add) counts[i] <- counts[i] + 1
  }
  
  return(counts)
}

many_stats = tibble()
for (i in 1:100) {
  stats = plot_data %>% distinct(snp, base, focus, group) %>%
    filter(base != 0)
  
  p = tibble()
  for (i in 1:nrow(stats)) {
    x = semi_join(filt.counts, stats[i,], by = join_by(snp, group))$count %>% split_four()
    
    expected = semi_join(genome_context, stats[i,], by = join_by(base, focus)) %>%
      with(setNames(perc, value))
    
    p %<>% bind_rows(c(bias = abs((x / sum(x)) - expected) %>% sum(),
                       p = chisq.test(x, p = expected)$p.value))
  }
  
  p$p %<>% p.adjust(method = 'fdr')
  stats %<>% bind_cols(p = p)
  
  many_stats %<>% bind_rows(stats)
}

ggplot(many_stats) +
  geom_hline(yintercept = 0.01) +
  geom_point(aes(x = bias, y = p, color = snp)) +
  scale_y_log10() +
  coord_cartesian(ylim = c(1e-100, 1)) +
  facet_grid(rows = vars(snp), cols = vars(group))

many_stats %<>% filter(p < 0.05) %>%
  group_by(group, snp) %>%
  summarise(thresh = min(bias), .groups = 'drop')


b = ggplot(filter(plot_data, base != 0), aes(x = base, y = bias, color = snp)) +
  geom_hline(data = many_stats, aes(yintercept = thresh), linetype = 'dashed', alpha = 0.5) +
  geom_vline(xintercept = 0, color = 'grey') +
  geom_point(aes(alpha = p <= 0.05)) +
  geom_label(data = filt.counts, aes(label = paste0('n=', count)),
             x = -6.4 , y = 1, color = 'black', fill = 'white', size = 3,
             label.size = 0, hjust = 0) +
  geom_label(data = filt.counts, aes(label = loss),
  x = 6.4 , y = 1, color = 'black', fill = 'white', size = 3,
  label.size = 0, hjust = 1) +
  facet_grid(rows = vars(snp), cols = vars(group), switch = 'y',
             labeller = as_labeller(function(x) sub(">", "→", x))) +
  scale_x_continuous(breaks = seq(-6, 6, 2), labels = c('-6', '-4', '-2', '0', '+2', '+4', '+6')) +
  scale_y_continuous(sec.axis = dup_axis(name = 'Absolute difference from expected frequencies'),
                     breaks = c(0, 0.5, 1)) +
  scale_color_manual(values = pal_spec, guide = NULL) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.1), guide = NULL) + 
  coord_cartesian(ylim = c(0, 1.2)) +
  labs(y = NULL, x = "Context Position", color = "SNP") +
  theme_bw() +
  theme(axis.text.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        strip.text.y.left = element_text(family = 'Aptos'),
        strip.background = element_rect(fill = NA, color = NA))


cowplot::plot_grid(plotlist = list(a, b), ncol = 1, labels = c('A', 'B'))
ggsave('figS5.svg', width = 7, height = 8.5)


#### Figure 2 () ####

snp.conv = rep(spectrum, each = 2)
names(snp.conv) = c('AG', 'TC', 'GA', 'CT', 'AC', 'TG', 'GT', 'CA', 'AT', 'TA', 'GC', 'CG')


all.runs = genome.motifs %>%
  # finds runs 3' of purine and pyrimidine
  mutate(
    base.3 = substr(motif, 12, 12),
    base.5 = substr(motif, 10, 10),
    run.length.3 = attr(regexpr("^(.)\\1*", substring(motif, 12, 21), perl = TRUE), "match.length"),
    run.length.5 = attr(regexpr("^(.)\\1*", substring(motif, 1, 10) %>% stri_reverse(), perl = TRUE), "match.length")
  ) %>%
  
  # extends run length if focal base matches the run base
  mutate(
    run.length.3 = case_when(
      focus == base.3 & focus == base.5 ~ run.length.3 + run.length.5 + 1,
      focus == base.3 ~ run.length.3 + 1,
      TRUE ~ run.length.3
    ),
    run.length.5 = case_when(
      focus == base.3 & focus == base.5 ~ run.length.3, # only works cause run.length.3 updates first
      focus == base.5 ~ run.length.5 + 1,
      TRUE ~ run.length.5
    )
  ) %>%
  
  # identifies run sites where transient misalignment (nascent or template) could happen
  mutate(
    cat.pur = case_when(
      run.length.3 == 1 ~ 'no',
      focus == base.3 & focus == base.5 ~ 'no',
      focus != base.3 ~ 'nascent',
      TRUE ~ 'template'),
    cat.pyr = case_when(
      run.length.5 == 1 ~ 'no',
      focus == base.3 & focus == base.5 ~ 'no',
      focus != base.5 ~ 'nascent',
      TRUE ~ 'template'
    )) %>%
  
  # does final update to give terminal run bases the full run length in 5' direction too
  mutate(
    run.length.3 = case_when(
      cat.pur == 'no' & focus == base.5 & focus != base.3 ~ run.length.5,
      TRUE ~ run.length.3
    ),
    run.length.5 = case_when(
      cat.pyr == 'no' & focus == base.3 & focus != base.5 ~ run.length.3,
      TRUE ~ run.length.5
    )
  ) %>%
  
  # identifies which snp would result from transient misalignment
  mutate(
    snp.pur = case_when(
      cat.pur == 'nascent' ~ snp.conv[paste0(focus, base.3)],
      cat.pur == 'template' ~ snp.conv[paste0(focus, base.5)],
      TRUE ~ NA_character_
    ),
    snp.pyr = case_when(
      cat.pyr == 'nascent' ~ snp.conv[paste0(focus, base.5)],
      cat.pyr == 'template' ~ snp.conv[paste0(focus, base.3)],
      TRUE ~ NA_character_
    )
  )

# diagnostic tibble
all.runs %>% select(motif, base.3, run.length.3, cat.pur, snp.pur, base.5, run.length.5, cat.pyr, snp.pyr)  %>%
  mutate(motif = paste0(substr(motif, 1, 10), '*', substr(motif, 11, 11), '*', substr(motif, 12, 21))) %>%
  print(n=20)


# pivot longer
all.runs %<>%
  rename(base.pur = base.3, base.pyr = base.5,
         run.length.pur = run.length.3, run.length.pyr = run.length.5) %>%
  pivot_longer(
    cols = c(base.pur, base.pyr,
             run.length.pur, run.length.pyr,
             cat.pur, cat.pyr,
             snp.pur, snp.pyr),
    names_to = c(".value", "template"),
    names_pattern = "(.*)\\.(pur|pyr)"
  )


# summarise counts
all.runs %<>%
  group_by(focus, template, run.length, cat, snp) %>%
  summarise(gen.count = n(), .groups = 'drop') %>%
  mutate(snp = factor(snp, levels = spectrum),
         cat = factor(cat, levels = c('no', 'template', 'nascent')))

# diagnostic plot
ggplot(all.runs, aes(x = run.length, y = gen.count, shape = cat, color = snp)) +
  # geom_bar(stat = 'identity', position = 'dodge') +
  geom_point() +
  scale_color_manual(values = pal_spec) +
  scale_y_log10() +
  facet_grid(cols = vars(template), rows = vars(snp)) +
  theme_bw()


# combines runs >6 into 6+ category
max.length = 6
all.runs %<>%
  filter(run.length < max.length) %>%
  bind_rows(filter(all.runs, run.length >= max.length) %>%
              group_by(focus, template, cat, snp) %>%
              summarise(gen.count = sum(gen.count), .groups = 'drop') %>%
              bind_cols(run.length = max.length))


### does the same for the mutation data

plot_data = my_data %>%
  mutate(focus = substr(snp, 1, 1)) %>%
  # finds runs 3' of purine and pyrimidine
  mutate(
    base.3 = substr(motif, 12, 12),
    base.5 = substr(motif, 10, 10),
    run.length.3 = attr(regexpr("^(.)\\1*", substring(motif, 12, 21), perl = TRUE), "match.length"),
    run.length.5 = attr(regexpr("^(.)\\1*", substring(motif, 1, 10) %>% stri_reverse(), perl = TRUE), "match.length")
  ) %>%
  
  # extends run length if focal base matches the run base
  mutate(
    run.length.3 = case_when(
      focus == base.3 & focus == base.5 ~ run.length.3 + run.length.5 + 1,
      focus == base.3 ~ run.length.3 + 1,
      TRUE ~ run.length.3
    ),
    run.length.5 = case_when(
      focus == base.3 & focus == base.5 ~ run.length.3, # only works cause run.length.3 updates first
      focus == base.5 ~ run.length.5 + 1,
      TRUE ~ run.length.5
    )
  ) %>%
  
  # identifies run sites where transient misalignment (nascent or template) could happen
  mutate(
    cat.pur = case_when(
      run.length.3 == 1 ~ 'no',
      focus == base.3 & focus == base.5 ~ 'no',
      focus != base.3 ~ 'nascent',
      TRUE ~ 'template'),
    cat.pyr = case_when(
      run.length.5 == 1 ~ 'no',
      focus == base.3 & focus == base.5 ~ 'no',
      focus != base.5 ~ 'nascent',
      TRUE ~ 'template'
    )) %>%
  
  # does final update to give run-internal bases the full run length
  mutate(
    run.length.3 = case_when(
      cat.pur == 'no' & focus == base.5 & focus != base.3 ~ run.length.5,
      TRUE ~ run.length.3
    ),
    run.length.5 = case_when(
      cat.pyr == 'no' & focus == base.3 & focus != base.5 ~ run.length.3,
      TRUE ~ run.length.5
    )
  ) %>%
  
  # identifies which snp would result from transient misalignment
  mutate(
    snpTM.pur = case_when(
      cat.pur == 'nascent' ~ snp.conv[paste0(focus, base.3)],
      cat.pur == 'template' ~ snp.conv[paste0(focus, base.5)],
      TRUE ~ NA_character_
    ),
    snpTM.pyr = case_when(
      cat.pyr == 'nascent' ~ snp.conv[paste0(focus, base.5)],
      cat.pyr == 'template' ~ snp.conv[paste0(focus, base.3)],
      TRUE ~ NA_character_
    )
  )


plot_data %<>%
  rename(base.pur = base.3, base.pyr = base.5,
         run.length.pur = run.length.3, run.length.pyr = run.length.5) %>%
  pivot_longer(
    cols = c(base.pur, base.pyr,
             run.length.pur, run.length.pyr,
             cat.pur, cat.pyr,
             snpTM.pur, snpTM.pyr),
    names_to = c(".value", "template"),
    names_pattern = "(.*)\\.(pur|pyr)"
  )


plot_data %<>%
  mutate(cat = case_when(
    cat == 'no' ~ cat,
    snp ==  snpTM ~ cat,
    TRUE ~ 'no'
  )) %>%
  group_by(group, snp, focus, template, run.length, cat) %>%
  summarise(count = n(), .groups = 'drop') %>%
  left_join(rates_summary %>% select(group, generations))


# combines runs >6 into 6+ category
plot_data %<>% filter(run.length < max.length) %>%
  bind_rows(filter(plot_data, run.length >= max.length) %>%
              group_by(group, snp, focus, template, cat, generations) %>%
              summarise(count = sum(count), .groups = 'drop') %>%
              bind_cols(run.length = max.length))


one_data = plot_data %>%
  filter(run.length <= 1) %>%
  group_by(group, snp, focus, template, generations) %>% summarise(count = sum(count), .groups = 'drop') %>%
  left_join(all.runs %>%
              filter(run.length <= 1) %>%
              group_by(focus, template) %>% summarise(gen.count = sum(gen.count), .groups = 'drop') %>%
              select(focus, template, gen.count)) %>%
  mutate(one.rate = count / generations / gen.count)


plot_data = bind_rows(
  plot_data %>%
    filter(cat == 'no' & run.length > 1) %>%
    left_join(all.runs %>% filter(cat == 'no') %>% select(-snp)),
  plot_data %>%
    filter(cat != 'no' & run.length > 1) %>%
    left_join(all.runs %>% filter(cat != 'no'))
) %>%
  mutate(rate = count / generations / gen.count)


plot_data %<>% left_join(one_data %>% select(group, snp, template, one.rate)) %>%
  mutate(rel.rate = rate / one.rate)


plot_data$cat %<>% factor(levels = c('nascent', 'template', 'no'))
plot_data$group %<>% factor(levels = groups)

plot_list = list() ; for (temp in c('pur', 'pyr')) {
  plot_list %<>% append(
    list(
      ggplot(plot_data %>% filter(template == temp & rel.rate >= 1),
             aes(x = run.length, y = rel.rate, fill = cat)) +
        geom_bar(stat = 'identity', position = position_dodge2(preserve = 'single')) +
        geom_hline(yintercept = 1) +
        scale_x_continuous(breaks = c(2:6), labels = c(2, 3, 4, 5, '6+')) +
        scale_y_log10(breaks = c(10, 100, 1000), labels = c('10<sup>1</sup>', '10<sup>2</sup>', '10<sup>3</sup>'),
                      sec.axis = dup_axis(name = ifelse(temp == 'pyr', 'Fold increase in mutation rate relative to no-run sites', ''))) +
        coord_cartesian(ylim = c(1, 5e3)) +
        scale_fill_manual(values = c('#56BCC2', '#E77D72', 'darkgrey'),
                          guide = NULL) +
        labs(x = 'Run length',
             y = '',
             title = paste('Assuming', c(pur = 'purine', pyr = 'pyrimidine')[temp], 'templated mispair')) +
        facet_grid(rows = vars(snp), cols = vars(group), switch = 'y',
                   labeller = as_labeller(function(x) sub('proofreading', 'Pr', x) %>% sub('>', '→', .))) +
        theme_bw() +
        theme(axis.text.y.right = if(temp == 'pur') {element_blank()} else{element_markdown()},
              axis.text.y.left = element_blank(),
              strip.text.y = if(temp == 'pur') {element_text(family = 'Aptos')} else {element_blank()},
              strip.background = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.border = element_blank(),
              axis.ticks.y = element_blank())
        # guides(fill = guide_legend(position = 'top', title = 'Transient\nmisalignment?'))
    )
  )
}



plot_grid(plotlist = plot_list, ncol = 2, labels = c('B', 'C'))

ggsave('fig2bc.svg', width = 9, height = 5)



#### Figure 3 / Figure S6 () ####

# run all this code with either chosen.snp = 'GC>CG' or 'AT>CG
chosen.snp = 'GC>CG'

c.runs.genome = genome.motifs %>%
  mutate(run.length = attr(regexpr("^C*", substring(motif, 12, 21), perl = TRUE), "match.length"),
         five.base = substr(motif, 10, 10)) %>%
  group_by(focus, five.base, run.length) %>%
  summarise(gen.count = n(), .groups = 'drop') %>%
  mutate(five.base = factor(five.base, levels = c('A', 'T', 'G', 'C'))) %>%
  filter(focus == substr(chosen.snp, 1, 1))


# fills dataframe with missing combinations, gives them NA value
c.runs.genome %<>% complete(focus, five.base, run.length)


if (chosen.snp == 'GC>CG') {breaks = list(c(4e5), c(1.5e5), c(5e4), c(8e3), c(1.5e3), c(3e2), c(40), c(3), c(1), c(NA), c(1))}

if (chosen.snp == 'AT>CG') {breaks = list(c(5e5), c(1.2e5), c(4e4), c(6e3), c(1e3), c(1.5e2), c(25), c(4), c(1))}


# makes y axis labels for counts
labels <- prettyNum(format(breaks, scientific = FALSE), big.mark = ",")
pad <- function(x, width) {
  npad <- width - nchar(x)
  paste0(strrep("\u2007", npad), x)
}
labels %<>% pad(max(nchar(labels)))



plot_list = list() ; for (i in unique(c.runs.genome$run.length)) {
  plot_data = filter(c.runs.genome, run.length == i)
  
  plot_list %<>% append(list(
    ggplot(plot_data, aes(x = five.base, y = gen.count)) +
      geom_hline(aes(yintercept = 0), color = 'grey') +
      geom_bar(aes(fill = five.base),
               position = position_dodge2(preserve = 'single'),
               # color = 'black', size = 0.3,
               stat = "identity") +
      scale_fill_manual(values = pal_bases, guide = NULL) +
      scale_y_continuous(breaks = breaks[[i+1]], labels = labels[[i+1]]) +
      facet_wrap(~run.length, ncol = 1, scales = 'free_y', strip.position = 'left') +
      labs(x = i, y = 'Count in genome') +
      theme_minimal() +
      theme(panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            strip.text = element_blank(),
            plot.background = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_text(color = '#4D4D4D', size = 9), 
            axis.text.y.left = element_text(family = 'Aptos Mono'))
  ))
}

c = plot_grid(plotlist = plot_list, nrow = 2)




plot_data = my_data %>%
  filter(snp == chosen.snp) %>%
  mutate(run.length = attr(regexpr("^C*", substring(motif, 12, 21), perl = TRUE), "match.length"),
         five.base = substr(motif, 10, 10)) %>%
  group_by(group, snp, five.base, run.length) %>%
  summarise(count = n(), .groups = 'drop')


plot_data %<>%
  left_join(rates_summary %>% select(group, generations)) %>%
  mutate(focus = substr(snp, 1, 1)) %>%
  left_join(c.runs.genome) %>%
  mutate(rate = count / generations / gen.count,
         group = factor(group, levels = groups),
         five.base = factor(five.base, levels = c('A', 'T', 'G', 'C')))


zero_data = plot_data %>%
  filter(run.length == 0)

plot_data %<>%
  filter(run.length > 0) %>%
  left_join(zero_data %>%
              select(group, snp, five.base, rate) %>%
              rename(zero.rate = rate)) %>%
  mutate(rel.rate = rate / zero.rate)


# calculates n= labels
labels = plot_data %>%
  group_by(group, snp, run.length) %>%
  summarise(count = sum(count), height = max(rel.rate) * 8, .groups = 'drop')


library(ggtext)
b = ggplot(plot_data, aes(x = run.length, y = rel.rate)) +
  geom_bar(aes(fill = five.base), stat = "identity",
           position = position_dodge2(preserve = "single")) +
  geom_richtext(data = labels, aes(x = run.length, y = height,
                                   label = paste0('<i>', count, '</i>')),
                inherit.aes = F, size = 3, label.size = 0, fill = NA) +
  facet_wrap(~group,
             # scales = 'free_x',
             ncol = 1) +
  scale_fill_manual(values = pal_bases) +
  scale_x_continuous(breaks = c(1:10)) +
  scale_y_log10(breaks = c(1e1, 1e2, 1e3, 1e4, 1e5),
                labels = paste0(rep(10, 5), '<sup>', c(1:5), '</sup>')) +
  coord_cartesian(xlim = c(0.8, max(plot_data$run.length))) +
  geom_hline(yintercept = 1) +
  labs(x = paste0('Number of consecutive Cs 3\' of focal ',
                  substr(chosen.snp, 1, 1), ' nucleotide'),
       y = 'Fold increase in mutation rate') +
  theme_bw() +
  guides(fill = guide_legend(position = 'top', title = '5\' nucleotide')) +
  theme(axis.text.y = element_markdown(),
        axis.title.x.bottom = element_text(hjust = 0),
        legend.text = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.ticks.x = element_blank(),
        # strip.text = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        legend.justification = 'left')


zero_data %<>%
  mutate(group = gsub('proofreading', '', group)) %>%
  mutate(group = gsub('MMR', '', group)) %>%
  mutate(group = factor(group, levels = c('(+) (+)', '(+) (-)', '(-) (+)', '(-) (-)')))


limits <- if (chosen.snp == "GC>CG") {
  c(7e-11, 3.9e-10, 7e-9, 7e-9)
} else {
  c(8e-11, 3e-10, 4e-8, 12e-9)
}

labels <- zero_data %>%
  group_by(group, snp, run.length) %>%
  summarise(count = sum(count),
            height = max(rate) * 1.1,
            .groups = "drop") %>%
  bind_cols(limit.up = limits) %>%
  mutate(
    height = limit.up,
    limit.down = -0.1 * limit.up
  )


a = ggplot(zero_data, aes(x = run.length, y = rate)) +
  geom_hline(data = labels, aes(yintercept = limit.up), size = 0) +
  geom_hline(data = labels, aes(yintercept = limit.down), size = 0) +
  geom_hline(aes(yintercept = 0), color = 'grey') +
  geom_bar(aes(fill = five.base), position = 'dodge', stat = "identity",
           # color = 'black', size = 0.3
  ) +
  geom_richtext(data = labels, aes(x = run.length, y = height,
                                   label = paste0('*n = ', count, '*')),
                size = 3, label.size = 0, label.padding = unit(0, 'pt'),
                vjust = 0.9) +
  facet_wrap(~group,
             scales = 'free_y',
             ncol = 1) +
  coord_cartesian(clip = "on", expand = FALSE) +
  # facet_wrap(~snp + group, ncol = 4, scales = 'free_y') +
  scale_fill_manual(values = pal_bases) +
  scale_x_continuous(breaks = c(0:6), labels = c(paste(0:5), '6+')) +
  scale_y_continuous(labels = function(x) sci_label(x, digits = 0)) +
  labs(x = '', y = 'Mutation rate') +
  theme_bw() +
  guides(fill = guide_legend(position = 'top', title = NULL, override.aes = list(color = NA, fill = NA))) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        # axis.ticks.x = element_blank(),
        strip.text = element_text(color = 'white'),
        strip.background = element_blank(), 
        legend.text = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_markdown(),
        axis.text.y.left = element_markdown())


final_plot <- ggdraw() +
  draw_plot(a, x = 0, y = 0.205, width = 0.2, height = 0.8) +
  draw_plot(b, x = 0.2, y = 0.205, width = 0.8, height = 0.8) +
  draw_plot(c, x = 0.05, y = 0.015, width = 0.95, height = 0.18) +
  # Add labels manually
  draw_plot_label(c("A", "B", "C"), 
                  x = c(0, 0.23, 0), 
                  y = c(0.96, 0.96, 0.2), 
                  size = 14) +
  draw_plot_label(c('Count in genome', paste0('Number of consecutive Cs 3\' of focal ', substr(chosen.snp, 1, 1), ' nucleotide')),
                  x = c(0.04, 0.5),
                  y = c(0.1, 0.005),
                  fontface = 'plain',
                  size = 11,
                  angle = c(90, 0),
                  hjust = 0.5,
                  vjust = 0)



final_plot
ggsave('fig3.svg', height = 10, width = 7)
ggsave('figS6.svg', height = 10, width = 7)




#### Figure 4 () ####

### Fig 4A

genome.strands = c(
  c(mg1655[OriC : length(mg1655)],
    mg1655[1 : (Ter - 1)]) %>% as.character(),
  comp[mg1655[Ter : (OriC - 1)]]
) %>%
  table() %>%
  as_tibble() %>%
  mutate(strand = ifelse(. %in% c('A', 'G'), 'LDST', 'LGST')) %>%
  group_by(strand) %>% summarise(genome.count = sum(n))


plot_data = my_data %>%
  group_by(group, snp, strand) %>%
  summarise(count = n()) %>%
  left_join(genome.strands) %>%
  left_join(rates_summary %>% select(group, generations)) %>%
  mutate(rate = count / genome.count / generations / 2) %>% # divides by 2 to make equivalent with non-split graph
  mutate(group = factor(group, levels = groups)) %>%
  mutate(strand = factor(strand, levels = c('LDST', 'LGST')))

stats = plot_data %>%
  select(group, snp, strand, rate) %>%
  pivot_wider(names_from = strand, values_from = rate) %>%
  mutate(
    total = LDST + LGST,
    LDST_pct = round(LDST / total * 100),
    LGST_pct = 100 - LDST_pct,
    ratio = paste0(LDST_pct, ":", LGST_pct),
    ymax = max(LDST, LGST)
  ) %>%
  select(group, snp, ratio, ymax)

stats %<>%
  left_join(stats %>%
              group_by(group) %>%
              summarise(max = max(ymax)))

limits = tibble(
  group = groups,
  lim = c(6e-11, 2e-8, 6e-7, 6e-7)
)

ggplot(plot_data) +
  geom_bar(aes(x = snp, y = rate, fill = strand, color = snp),
           stat = 'identity', position = 'dodge') +
  geom_label(data = stats, aes(x = snp, y = ymax + (0.2 * max), label = ratio), size = 2,
             label.size = 0, label.padding = unit(0, 'points')) +
  geom_hline(data = limits, aes(yintercept = lim), size = 0) +
  scale_x_discrete(labels = function(x) sub('>', '→', x)) +
  # scale_y_continuous(labels = function(x) sci_label(x, digits = 1)) +
  scale_fill_manual(values = c('white', 'darkgrey'), labels = c('Leading strand', 'Lagging strand')) +
  scale_color_manual(values = pal_spec, guide = NULL) +
  facet_wrap(~group, scales = 'free_y', nrow = 1) +
  labs(x = NULL, y = 'Mutation rate', fill = 'Mutation rate at base pairs where purine templates the:') +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = NA, color = NA),
        axis.title.y = element_markdown(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, family = 'Aptos'),
        axis.text.y = element_markdown()) +
  guides(fill = guide_legend(position = 'top', order = 1, 
                             override.aes = list(color = 'black')))

breaks = list(c(6e-11, 3e-11, 0), c(2e-8, 1e-8, 0), c(6e-7, 3e-7, 0), c(6e-7, 3e-7, 0))

plot_list = list() ; for(i in 1:4) {
  plot_data2 = filter(plot_data, group == groups[i])
  stats2 = filter(stats, group == groups[i])
  
  plot_list = append(plot_list, list(
    ggplot(plot_data2) +
      geom_bar(aes(x = snp, y = rate, fill = strand, color = snp),
               stat = 'identity', position = 'dodge') +
      geom_label(data = stats2, aes(x = snp, label = ratio), y = stats2$ymax + 0.1 * max(breaks[[i]]), size = 2,
                 label.size = 0, label.padding = unit(0, 'points')) +
      scale_x_discrete(labels = function(x) sub('>', '→', x)) +
      scale_y_continuous(labels = function(x) sci_label(x, digits = 0),
                         breaks = breaks[[i]]) +
      scale_fill_manual(values = c('white', 'darkgrey'),
                        labels = c('Leading strand', 'Lagging strand'),
                        guide = NULL) +
      scale_color_manual(values = pal_spec, guide = NULL) +
      coord_cartesian(ylim = c(NA, limits$lim[i])) +
      facet_wrap(~group, scales = 'free_y', nrow = 1) +
      labs(x = NULL, y = 'Mutation rate', fill = 'Mutation rate at base pairs where purine templates the:') +
      theme_bw() +
      theme(panel.grid.major.x = element_blank(),
            strip.background = element_rect(fill = NA, color = NA),
            strip.text = element_text(size = 8),
            panel.border = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            axis.text.x = element_markdown(angle = 90, vjust = 0.4, family = 'Aptos Mono', size = 10,
                                           margin = margin(t = 3)),
            axis.text.y = element_markdown())
  ))
}


legend = ggplot(plot_data, aes(x = snp, y = rate, fill = strand)) +
  geom_bar(stat = 'identity', color = 'black') +
  scale_fill_manual(values = c('white', 'darkgrey'),
                    labels = c('Leading strand', 'Lagging strand'),
                    name = 'Mutation rate at base pairs where purine templates the:',
                    guide = guide_legend(nrow = 1, title.position = 'left'))

legend = get_legend(legend)


a = ggdraw() +
  draw_plot(legend,
            x = 0, y = 0.9, width = 1, height = 0.1) +
  draw_plot(plot_grid(plotlist = plot_list, nrow = 1),
            x = 0.02, y = 0.05, width = 0.98, height = 0.82) +
  draw_plot_label('Mutation rate', x = 0, y = 0.35, angle = 90, hjust = 0,
                  fontface = 'plain', size = 12)




### Fig 4B
plot_data = my_data %>%
  mutate(base = list(bases)) %>%
  unnest(base) %>%
  mutate(value = substr(motif, 11 + base, 11 + base))

plot_data$value %<>% factor(levels = c('A', 'T', 'G', 'C'))

plot_data %<>%
  group_by(group, snp, strand, base, value) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(!is.na(count)) %>%
  complete(group, snp, strand, base, value, fill = list(count = 0)) # Fill missing combinations


plot_data %<>%
  mutate(focus = substr(snp, 1, 1)) %>% # adds col for focus base (A or G)
  left_join(genome_context %>% rename(gen.perc = perc)) %>% # adds info for the genome_context bias at each position
  mutate(weight.count = (0.25 / gen.perc) * count) # weights by genome structure


plot_data %<>% left_join(
  plot_data %>%
    group_by(group, snp, strand, base) %>%
    summarise(weight.total = sum(weight.count), .groups = 'drop')
) %>%
  mutate(weight.perc = weight.count / weight.total)

plot_data %<>% mutate(weight.perc = ifelse(base == 0 & value == 'G', 1, weight.perc))


b = ggplot(filter(plot_data, snp == 'GC>CG'), aes(x = base, y = weight.perc, fill = value)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = pal_bases, name = 'Base', guide = NULL) +
  scale_y_continuous(breaks = c(0.25, 0.5, 0.75), sec.axis = dup_axis(name = 'Weighted frequency')) +
  scale_x_continuous(breaks = seq(-6, 6, 2), labels = c('-6', '-4', '-2', '0', '+2', '+4', '+6')) +
  geom_hline(yintercept = c(0.5), size = 0.5) +
  geom_hline(yintercept = c(0.25, 0.75), size = 0.25) +
  facet_grid(rows = vars(strand), cols = vars(group), switch = 'y',
             labeller = labeller(strand = c(
               'LDST' = 'Leading strand',
               'LGST' = 'Lagging strand'
             ))) +
  coord_cartesian(xlim = c(-6, 6)) +
  labs(y = 'G:C→C:G sites\nwhere G templates the:', x = NULL, fill = "Value") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y.left = element_blank(),
        axis.title.y.left = element_text(family = 'Aptos'))


### Fig 4C

plot_data %<>%
  filter(base != 0) %>%
  select(group, snp, base, value, strand, weight.perc) %>%
  pivot_wider(names_from = strand, values_from = weight.perc)


regression_data <- plot_data %>%
  group_by(snp, group) %>%
  summarise({
    model <- lm(LGST ~ LDST)
    coefs <- summary(model)$coefficients
    
    # car::linearHypothesis test: is slope == 1?
    lh <- car::linearHypothesis(model, "LDST = 1")
    
    tibble(
      intercept = coefs["(Intercept)", "Estimate"],
      intercept_p = coefs["(Intercept)", "Pr(>|t|)"],
      slope = coefs["LDST", "Estimate"],
      slope_p = coefs["LDST", "Pr(>|t|)"],
      r_squared = summary(model)$r.squared,
      slope_test_F = lh$F[2],       # F statistic from car output
      slope_test_p = lh$`Pr(>F)`[2] # p-value for H0: slope = 1
    )
  }, .groups = "drop")

regression_data %<>%
  mutate(slope_test_p = p.adjust(slope_test_p, method = 'fdr'))

# Merge regression data with the original data
plot_data <- plot_data %>%
  left_join(regression_data, by = c("snp", "group"))

# plot_data %<>% filter(snp == 'GC>CG')

# Plot with manual annotations
c = ggplot(plot_data %>% filter(snp == 'GC>CG'),
           aes(x = LDST, y = LGST, fill = value)) +
  geom_hline(yintercept = 0.25) +
  geom_vline(xintercept = 0.25) +
  geom_abline(slope = 1, intercept = 0, color = 'darkgrey', alpha = 0.5) +
  geom_abline(aes(intercept = intercept, slope = slope), color = "black", linetype = "dashed", size = 0.75) +
  geom_point(color = 'black', shape = 21, stroke = 0.1, size = 2.5) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  coord_cartesian(xlim = c(0, max(c(plot_data$LDST, plot_data$LGST))), ylim = c(0, max(c(plot_data$LDST, plot_data$LGST)))) +
  scale_fill_manual(values = pal_bases, guide = NULL) +
  facet_grid(rows = vars(snp), cols = vars(group),
             labeller = as_labeller(function(x) sub(">", "→", x))) +
  labs(x = 'Nucleotide frequency across G:C→C:G sites where G templates the... leading strand',
       y = '... lagging strand') +
  theme_minimal() +
  theme(
    strip.text.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(family = 'Aptos'),
    axis.title.x = element_text(family = 'Aptos'))


### Fig 4D

d = ggplot(regression_data) +
  geom_bar(aes(x = snp, y = slope, fill = snp), stat = 'identity') +
  geom_hline(yintercept = 1, size = 0.5, color = 'black') +
  geom_text(aes(x = snp, y = slope - 0.12,
                label = ifelse(slope_test_p < 0.001, '***',
                               ifelse(slope_test_p < 0.01, '**',
                                      ifelse(slope_test_p < 0.05, '*', '')))),
            color = 'white'
  ) +
  scale_fill_manual(values = pal_spec, guide = NULL) +
  scale_x_discrete(labels = function(x) sub('>', '→', x)) +
  scale_y_continuous(name = 'Slope', breaks = c(0, 0.5, 1)) +
  facet_wrap(~group, nrow = 1) +
  labs(x = NULL, fill = NULL) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        legend.text = element_markdown(size = 12, family = 'Aptos'),
        axis.title.y = element_markdown(),
        axis.text.x = element_markdown(angle = 90,
                                       vjust = 0.4,
                                       family = 'Aptos Mono',
                                       size = 10))




plot_grid(a, b, c, d, ncol = 1, rel_heights = c(1, 1, 0.87, 0.6), labels = 'AUTO')

scale = 1
ggsave('fig4.svg', width = 8*scale, height = 10*scale)



#### Supp. Figure 7 (Strand-comparison dot plots for all BPSs) ####

# use plot_data created in above section (Fig 4C) 


ggplot(plot_data, aes(x = LDST, y = LGST, fill = value)) +
  geom_hline(yintercept = 0.25) +
  geom_vline(xintercept = 0.25) +
  geom_abline(slope = 1, intercept = 0, color = 'darkgrey', alpha = 0.5) +
  geom_abline(aes(intercept = intercept, slope = slope), color = "black", linetype = "dashed", size = 0.75) +
  geom_point(color = 'black', shape = 21, stroke = 0.1, size = 2.5) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  coord_cartesian(xlim = c(0, max(c(plot_data$LDST, plot_data$LGST))), ylim = c(0, max(c(plot_data$LDST, plot_data$LGST)))) +
  scale_fill_manual(values = pal_bases) +
  facet_grid(rows = vars(snp), cols = vars(group),
             labeller = as_labeller(function(x) sub(">", "→", x))) +
  labs(x = 'Nucleotide frequency across sites where purine templates the leading strand', y = 'Nucleotide frequency across sites where purine templates the lagging strand') +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    strip.text.y = element_text(family = 'Aptos'),
  ) +
  guides(fill = guide_legend(position = 'left', title = 'Nucleotide'))

ggsave('figS7.svg', width = 7*1.1, height = 8.6*1.1)


### Supplementary Table 2

regression_data %>%
  select(group, snp, slope, slope_test_p, r_squared) %>%
  arrange(group, snp) %>%
  mutate(snp = sub('>', '→', snp)) %>%
  write_csv('Table S2.csv')
  



#### Supp. Figure 8 (Strand bias of A:T→C:G and G:C→C:G hotspots) ####


OriC <- 3925696
Ter <- 1640202


c.runs.genome = genome.motifs %>%
  mutate(replichore = ifelse(pos >= Ter & pos < OriC, 'left', 'right'),
         strand = case_when(
           (replichore == 'left' & ref %in% c('A', 'G'))
           | (replichore == 'right' & ref %in% c('T', 'C')) ~ 'LDST',
           TRUE ~ 'LGST'),
         run.length = attr(regexpr("^C*", substring(motif, 12, 21), perl = TRUE), "match.length")) %>%
  group_by(focus, strand, run.length) %>%
  summarise(gen.count = n(), .groups = 'drop')


plot_data = my_data %>%
  mutate(run.length = attr(regexpr("^C*", substring(motif, 12, 21), perl = TRUE), "match.length"),
         focus = substr(snp, 1, 1)) %>%
  group_by(group, snp, focus, strand, run.length) %>%
  summarise(count = n(), .groups = 'drop')


plot_data %<>%
  left_join(rates_summary %>% select(group, generations)) %>%
  left_join(c.runs.genome) %>%
  mutate(rate = count / generations / gen.count,
         group = factor(group, levels = groups))


zero_data = plot_data %>%
  filter(run.length == 0)

plot_data %<>%
  left_join(zero_data %>% rename(zero.rate = rate) %>% select(group, snp, strand, zero.rate)) %>%
  mutate(rel.rate = rate / zero.rate)

zero_data %<>%
  mutate(group = gsub('proofreading', '', group)) %>%
  mutate(group = gsub('MMR', '', group)) %>%
  mutate(group = factor(group, levels = c('(+) (+)', '(+) (-)', '(-) (+)', '(-) (-)')))



plot_list = list() ; for (chosen.snp in c('GC>CG', 'AT>CG')) {
  
  plot_data2 = plot_data %>% filter(snp == chosen.snp)
  zero_data2 = zero_data %>% filter(snp == chosen.snp)
  
  labels <- zero_data2 %>%
    group_by(group, snp, run.length) %>%
    summarise(count = sum(count), height = max(rate) * 1.15, .groups = "drop") %>%
    mutate(limit.up = c(8e-11, 3e-10, 3.1e-8, 8e-9))
  
  
  a = ggplot(zero_data2, aes(x = run.length, y = rate)) +
    geom_hline(data = labels, aes(yintercept = limit.up), size = 0) +
    geom_hline(aes(yintercept = 0), color = 'grey') +
    geom_bar(aes(fill = strand), position = 'dodge', stat = "identity",
             color = 'black') +
    # geom_richtext(data = labels, aes(x = run.length, y = height,
    #                                  label = paste0('*n = ', count, '*')),
    #               size = 3, label.size = 0, label.padding = unit(0, 'pt'),
    #               vjust = 0.9) +
    facet_wrap(~group,
               scales = 'free',
               ncol = 1) +
    coord_cartesian(clip = "on", expand = FALSE) +
    # facet_wrap(~snp + group, ncol = 4, scales = 'free_y') +
    scale_fill_manual(values = c('white', 'grey')) +
    scale_x_continuous(breaks = c(0:6), labels = c(paste(0:5), '6+')) +
    scale_y_continuous(labels = function(x) sci_label(x, digits = 0)) +
    labs(x = '', y = 'Mutation rate') +
    theme_bw() +
    guides(fill = guide_legend(position = 'top', title = NULL, override.aes = list(color = NA, fill = NA))) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          # axis.ticks.x = element_blank(),
          strip.text = element_text(color = 'white'),
          strip.background = element_blank(), 
          legend.text = element_blank(),
          panel.border = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_markdown(),
          axis.text.y.left = element_markdown())
  
  # calculates n= labels
  labels = plot_data2 %>%
    group_by(group, snp, run.length) %>%
    summarise(count = sum(count), height = max(rel.rate) * 8, .groups = 'drop')
  
  
  b = ggplot(plot_data2, aes(x = run.length, y = rel.rate)) +
    geom_bar(aes(fill = strand), stat = "identity",
             position = position_dodge2(preserve = "single"),
             color = 'black') +
    # geom_richtext(data = labels, aes(x = run.length, y = height,
    #                                  label = paste0('<i>', count, '</i>')),
    #               inherit.aes = F, size = 3, label.size = 0, fill = NA) +
    facet_wrap(~group,
               scales = 'free_x',
               ncol = 1) +
    scale_fill_manual(values = c('white', 'grey'), labels = c('Leading strand', 'Lagging strand')) +
    scale_x_continuous(breaks = c(1:10)) +
    scale_y_log10(breaks = c(1e1, 1e2, 1e3, 1e4, 1e5),
                  labels = paste0(rep(10, 5), '<sup>', c(1:5), '</sup>')) +
    coord_cartesian(xlim = c(0.8, 8), ylim = c(0.5, 1e5)) +
    geom_hline(yintercept = 1) +
    labs(x = paste0('Number of consecutive Cs 3\' of focal ',
                    substr(chosen.snp, 1, 1)),
         y = 'Fold increase in mutation rate') +
    theme_bw() +
    guides(fill = guide_legend(position = 'top',
                               title = paste(sub('>', '→', chosen.snp),
                                             'sites where',
                                             substr(chosen.snp, 1, 1),
                                             'templates the: '))) +
    theme(axis.text.y = element_markdown(),
          axis.title.x.bottom = element_text(hjust = 0),
          legend.text = element_text(size = 10),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.y = element_blank(),
          # axis.ticks.x = element_blank(),
          # strip.text = element_blank(),
          panel.border = element_blank(),
          strip.background = element_blank(),
          legend.title = element_text(family = 'Aptos'),
          legend.justification = 'left')
  
  plot_list %<>% append(list(plot_grid(a, b, rel_widths = c(0.35, 1))))
}

plot_grid(plotlist = plot_list, labels = 'AUTO')

ggsave('figS8.svg', width = 9, height = 8)






#### Figure 5 / Figure S10 () ####

# GC function for considering context with purine as focal nucleotide (20 bp sliding window)
window = 20
GC_context <- function(sequence, position, context, base) {
  if (base %in% c('A', 'G')) {
    if (context > 0) {
      indices = (((position+context) : (position+context+window-1)) - 1) %% length(sequence) + 1
    } else {
      indices = (((position+context-window+1) : (position+context)) - 1) %% length(sequence) + 1
    }
  } else {
    if (context > 0) {
      indices = (((position-context-window+1) : (position-context)) - 1) %% length(sequence) + 1
    } else {
      indices = (((position-context) : (position-context+window-1)) - 1) %% length(sequence) + 1
    }
  }
  # sequence[indices]
  length(which(sequence[indices] %in% c('G', 'C'))) / (window)
}


contexts <- c(1:1000)
contexts %<>% c(-contexts)

plot_data = tibble()
for (i in 1:length(contexts)) {
  plot_data %<>%
    bind_rows(
      my_data %>%
        mutate(context = contexts[i]) %>%
        rowwise() %>%
        mutate(GC = GC_context(mg1655, pos, context, ref)) %>%
        ungroup() %>%
        group_by(group, snp, context) %>%
        summarise(GC_mean = mean(GC), GC_sd = sd(GC), n = n(), .groups = "drop")
    )
  if (i %% 10 == 0) {
    print(i)
  }
}

window1000 = plot_data
write_csv(window1000, 'sliding_window.csv')

plot_data = read_csv('sliding_window.csv')

plot_data %<>%
  mutate(ts.tv = case_when(
    snp %in% c('AT>GC', 'GC>AT') ~ 1,
    snp %in% c('AT>CG', 'GC>TA') ~ 2,
    snp %in% c('AT>TA', 'GC>CG') ~ 2
  ))

custom_label <- function(x, digits = 0) {
  if (all(is.na(x))) return(NA_character_)
  sapply(x, function(val) {
    if (is.na(val)) return(NA_character_)
    if (val == 0) return("0")
    if (val < 99) return(as.character(val))
    if (val %in% c(1.5e-10, 0.5e-10, 15e-9, 5e-9)) return('')
    # if (val < 10) return(paste(paste0(rep('&nbsp;', 4 + digits), collapse = '') , as.character(val)))  
    sci <- format(val, scientific = TRUE)    # e.g. "4e+05"
    if (grepl("e", sci)) {
      parts <- strsplit(sci, "e")[[1]]
      base <- formatC(as.numeric(parts[1]), format = "f", digits = digits)
      exp  <- as.integer(parts[2])
      paste0(base, "×10<sup>", exp, "</sup>")
    } else {
      formatC(val, format = "f", digits = digits)
    }
  })
}


plots = list()
for (j in 1:2) {
  for (i in 1:4) {
    plot_data2 = plot_data %>% filter(group == groups[i] & ts.tv == j)
    
    p1 = ggplot(plot_data2 %>% filter(context < 0),
                aes(x = context, y = GC_mean, group = snp, color = snp)) +
      geom_hline(yintercept = length(which(mg1655 %in% c('G', 'C'))) / length(mg1655),
                 color = 'black', linetype = 'dashed') +
      geom_smooth(method = 'loess', span = 1, se = F) +
      geom_ribbon(aes(ymin = GC_mean - GC_sd / sqrt(n), ymax = GC_mean + GC_sd / sqrt(n), fill = snp),
                  linetype = 0, alpha = 0.3) +
      scale_fill_manual(values = pal_spec,
                        labels = c('A:T→G:C', 'G:C→A:T', 'A:T→C:G', 'G:C→T:A', 'A:T→T:A', 'G:C→C:G'),
                        guide = NULL) +
      scale_color_manual(values = pal_spec, guide = NULL) +
      scale_x_log10(trans = pseudo_log_trans(base = 10), breaks= c(-1000, -100, -10, -1)) +
      coord_cartesian(ylim = c(0.425, 0.575), xlim = c(-1000, -1)) +
      labs(y = NULL, x = NULL) +
      theme_bw() +
      theme(panel.grid.minor = element_blank(), panel.background = element_blank(),
            panel.border = element_blank(), plot.background = element_blank(),
            axis.text.y.left = if (i == 1) element_text() else element_blank(),
            axis.ticks.y.left = if (i == 1) element_line() else element_blank(),
            axis.text.x.bottom = if (j == 2) element_markdown() else element_blank(),
            axis.ticks.x.bottom = if (j == 2) element_line() else element_blank(),
            axis.line.x = element_line(linewidth = 0.3),
            axis.line.y.left = element_line(linewidth = 0.3),
            axis.line.y.right = element_blank())
    
    p2 = ggplot(plot_data2 %>% filter(context > 0),
                aes(x = context, y = GC_mean, group = snp, color = snp)) +
      geom_hline(yintercept = length(which(mg1655 %in% c('G', 'C'))) / length(mg1655),
                 color = 'black', linetype = 'dashed') +
      geom_smooth(method = 'loess', span = 1, se = F) +
      geom_ribbon(aes(ymin = GC_mean - GC_sd / sqrt(n), ymax = GC_mean + GC_sd / sqrt(n), fill = snp),
                  linetype = 0, alpha = 0.3) +
      scale_fill_manual(values = pal_spec,
                        labels = c('A:T→G:C', 'G:C→A:T', 'A:T→C:G', 'G:C→T:A', 'A:T→T:A', 'G:C→C:G'),
                        guide = NULL) +
      scale_color_manual(values = pal_spec, guide = NULL) +
      scale_y_continuous(position = 'right') +
      scale_x_log10(trans = pseudo_log_trans(base = 10), breaks = c(1, 10, 100, 1000)) +
      coord_cartesian(ylim = c(0.425, 0.575), xlim = c(1, 1000)) +
      labs(y = NULL, x = NULL) +
      theme_bw() +
      theme(panel.grid.minor = element_blank(), panel.background = element_blank(),
            panel.border = element_blank(), plot.background = element_blank(),
            axis.text.y.right = if (i == 4) element_text() else element_blank(),
            axis.ticks.y.right = if (i == 4) element_line() else element_blank(),
            axis.text.x.bottom = if (j == 2) element_text() else element_blank(),
            axis.ticks.x.bottom = if (j == 2) element_line() else element_blank(),
            axis.line.x = element_line(linewidth = 0.3),
            axis.line.y.left = element_blank(),
            axis.line.y.right = element_line(linewidth = 0.3))
    
    if (i == 1) {my_plot = plot_grid(p1, p2, rel_widths = c(1.12, 0.88))}
    else if (i == 4) {my_plot = plot_grid(p1, p2, rel_widths = c(0.88, 1.12))}
    else {my_plot = plot_grid(p1, p2)}
    
    plots %<>% append(list(my_plot))
  }
}

plot = plot_grid(plotlist = plots, ncol = 4,
                 rel_widths = c(1.13, 1, 1, 1.13),
                 rel_heights = c(1, 1.05))


plots_spaced <- plots
for (i in seq_along(plots_spaced)) {
  plots_spaced[[i]] <- plots_spaced[[i]] +
    theme(plot.margin = margin(r = 4, l = 4))   # adjust as needed
}

plot = plot_grid(
  plotlist = plots_spaced,
  ncol = 4,
  rel_widths = c(1.13, 1, 1, 1.13),
  rel_heights = c(1, 1.09)
)

plot = ggdraw() +
  draw_plot(plot, x = 0, y = 0.035, width = 1, height = 0.94) +
  draw_text(groups, x = c(0.15, 0.38, 0.62, 0.85), y = 0.99, hjust = 0.5, size = 10) +
  draw_text('Starting position of 20 bp sliding context window', x = 0.5, y = 0.025, size = 12)


a = ggdraw() +
  draw_plot(plot, x = 0.025, y = 0, width = 0.975, height = 0.975) + 
  draw_label('Mean GC content of context window', x = 0.015, y = 0.5, angle = 90, size = 12)

ggsave('plot.png', width = 9, height = 4)


### Figure 5B

window = 100

wrap.genome.wide = as.character(mg1655) %>% paste0(collapse = '')

wrap.genome.wide = paste0(
  substr(wrap.genome.wide, gen.length - window + 1, gen.length),
  wrap.genome.wide,
  substr(wrap.genome.wide, 1, window)
)

genome.motifs = tibble(
  pos = 1:length(mg1655),
  ref = as.character(mg1655)
  ) %>%
  mutate(motif = substring(wrap.genome.wide, pos + 100 - window, pos + 100 + window),
         five.context = case_when(
           ref %in% c('A', 'G') ~ substr(motif, 1, window),
           TRUE ~ substr(motif, window+2, 2*window+1)),
         three.context = case_when(
           ref %in% c('A', 'G') ~ substr(motif, window+2, 2*window+1),
           TRUE ~ substr(motif, 1, window))
         ) %>%
  mutate(focus = ifelse(ref %in% c('A', 'T'), 'A', 'G'),
         five.GC = nchar(gsub('[^GC]', '', five.context)) / window,
         three.GC = nchar(gsub('[^GC]', '', three.context)) / window) %>%
  select(-motif, -five.context, -three.context)


ggplot(genome.motifs) +
  geom_histogram(aes(y = five.GC, x = -after_stat(density), fill = focus),
                 position = 'identity', binwidth = 0.01, alpha = 0.5,
                 color = 'black', linewidth = 0.3) +
  geom_histogram(aes(y = three.GC, x = after_stat(density), fill = focus),
                 position = 'identity', binwidth = 0.01, alpha = 0.5,
                 color = 'black', linewidth = 0.3) +
  geom_hline(yintercept = length(which(mg1655 %in% c('G', 'C'))) / length(mg1655),
             color = 'black', linetype = 'dashed', linewidth = 0.7) +
  scale_fill_manual(values = c('blue', 'orange')) +
  labs(x = 'Probability Density (left is 5\' context, right is 3\' context)',
       y = 'GC-content of ±100 bp context',
       fill = 'Focal nucleotide') +
  theme_minimal() +
  guides(fill = guide_legend(position = 'top'))


plot_data = my_data %>%
  mutate(motif = substring(wrap.genome.wide, pos + 100 - window, pos + 100 + window),
         five.context = case_when(
           ref %in% c('A', 'G') ~ substr(motif, 1, window),
           TRUE ~ substr(motif, window+2, 2*window+1)),
         three.context = case_when(
           ref %in% c('A', 'G') ~ substr(motif, window+2, 2*window+1),
           TRUE ~ substr(motif, 1, window))
  ) %>%
  mutate(focus = ifelse(ref %in% c('A', 'T'), 'A', 'G'),
         five.GC = nchar(gsub('[^GC]', '', five.context)) / window,
         three.GC = nchar(gsub('[^GC]', '', three.context)) / window) %>%
  select(-motif, -five.context, -three.context)


base_plot = bind_rows(
  genome.motifs %>% mutate(snp = ifelse(focus == 'A', 'AT>GC', 'GC>AT')),
  genome.motifs %>% mutate(snp = ifelse(focus == 'A', 'AT>CG', 'GC>TA')),
  genome.motifs %>% mutate(snp = ifelse(focus == 'A', 'AT>TA', 'GC>CG'))
) %>%
  mutate(snp = factor(snp, levels = spectrum))

plot_data %<>% pivot_longer(cols = c(five.GC, three.GC), names_to = 'context', values_to = 'GC') 
base_plot %<>% pivot_longer(cols = c(five.GC, three.GC), names_to = 'context', values_to = 'GC') 


stats = tibble()
for (i in 1:4) {
  for (j in 1:6) {
    for (k in 1:2) {
      data1 = filter(plot_data, group == groups[i] & snp == spectrum[j] & context == c('five.GC', 'three.GC')[k])
      data2 = filter(base_plot, snp == spectrum[j] & context == c('five.GC', 'three.GC')[k])
      
      ks.result = ks.test(data1$GC, data2$GC)
      
      stats %<>% bind_rows(tibble(
        group = groups[i], snp = spectrum[j], context = c('five.GC', 'three.GC')[k],
        p = ks.result$p.value, D = ks.result$statistic
      ))
    }
  } ; print(i)
}

stats %<>% mutate(p.adj = p.adjust(p, method = 'fdr'))

stats %<>% left_join(my_data %>%
                       group_by(group, snp) %>%
                       summarise(count = paste0('n =  ', n())) %>%
                       mutate(context = 'five.GC'))


# run only for main text figure
stats %<>%
  mutate(label = case_when(
    p.adj >= 0.05 ~ paste0('N.S.<br>', 'D = ', round(D, 2)),
    p.adj >= 0.01 ~ paste0('<i>p</i> < 0.05<br>D = ', round(D, 2)),
    p.adj >= 0.001 ~ paste0('<i>p</i> < 0.01<br>D = ', round(D, 2)),
    TRUE ~ paste0('<i>p</i> < 0.001<br>D = ', round(D, 2))
  ))


# plot_data %<>% mutate(snp = factor(snp, levels = spectrum[c(1,3,5,2,4,6)]))
# base_plot %<>% mutate(snp = factor(snp, levels = spectrum[c(1,3,5,2,4,6)]))


plot_data %<>% filter(snp %in% spectrum[c(1, 4, 6)] & group %in% groups[c(1:2)])
base_plot %<>% filter(snp %in% spectrum[c(1, 4, 6)])
stats %<>% filter(snp %in% spectrum[c(1, 4, 6)] & group %in% groups[c(1:2)])


binwidth = 0.03
b = ggplot() +
  geom_histogram(data = base_plot %>% filter(context == 'five.GC'),
                 aes(y = GC, x = -after_stat(density)),
                 binwidth = binwidth, fill = 'white', color = 'black', linewidth = 0.5) +
  geom_histogram(data = base_plot %>% filter(context == 'three.GC'),
                 aes(y = GC, x = after_stat(density)),
                 binwidth = binwidth, fill = 'white', color = 'black', linewidth = 0.5) +
  geom_histogram(data = plot_data %>% filter(context == 'five.GC'),
                 aes(y = GC, x = -after_stat(density), fill = snp),
                 binwidth = binwidth, alpha = 0.7) +
  geom_histogram(data = plot_data %>% filter(context == 'three.GC'),
                 aes(y = GC, x = after_stat(density), fill = snp),
                 binwidth = binwidth, alpha = 0.7) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  geom_hline(yintercept = length(which(mg1655 %in% c('G', 'C'))) / length(mg1655),
             linetype = 'dashed', linewidth = 0.4) +
  geom_richtext(data = stats,
                aes(x = ifelse(context == 'five.GC', -4.5, 4.5), label = label),
                y = 0.2, label.padding = unit(0, 'in'), label.size = 0, hjust = 0.5, size = 3.5) +
  geom_richtext(data = stats,
                aes(x = -0.3, y = 0.82, label = count),
                label.size = 0, hjust = 0.5, size = 3.7, label.padding = unit(3, 'points')) +
  scale_y_continuous(breaks = c(0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(-8, -4, 0, 4, 8), labels = c(8, 4, 0, 4, 8)) +
  scale_fill_manual(values = pal_spec, guide = NULL) +
  # facet_wrap(~snp, ncol = 3) +
  facet_grid(cols = vars(snp), rows = vars(group)) +
  labs(x = 'Probability Density (left is 5\' context, right is 3\' context)',
       y = 'GC-content of ±100 bp context',
       fill = NULL) +
  theme_minimal() +
  theme(strip.text.x = element_blank(),
        axis.text.x = element_blank())


plot_grid(a, b, ncol = 1, labels = 'AUTO')

ggsave('fig5.svg', width = 9, height = 9)



### supplementary Figure 10

# run only for supplementary figure
stats %<>% mutate(label = case_when(
  p.adj < 1e-10 ~ paste0('*p* < 10<sup>-10</sup><br>D = ', round(D, 2)),
  p.adj >= 0.05 ~ paste0('N.S.<br>D = ', round(D, 2)),
  TRUE ~ paste0('*p* = ', sci_label(p.adj, digits = 1), '<br>D = ', round(D, 2))
))


plot_data %<>% mutate(snp = factor(snp, levels = spectrum))
base_plot %<>% mutate(snp = factor(snp, levels = spectrum))

binwidth = 0.03
ggplot() +
  geom_histogram(data = base_plot %>% filter(context == 'five.GC'),
                 aes(y = GC, x = -after_stat(density)),
                 binwidth = binwidth, fill = 'white', color = 'black', linewidth = 0.5) +
  geom_histogram(data = base_plot %>% filter(context == 'three.GC'),
                 aes(y = GC, x = after_stat(density)),
                 binwidth = binwidth, fill = 'white', color = 'black', linewidth = 0.5) +
  geom_histogram(data = plot_data %>% filter(context == 'five.GC'),
                 aes(y = GC, x = -after_stat(density), fill = snp),
                 binwidth = binwidth, alpha = 0.7) +
  geom_histogram(data = plot_data %>% filter(context == 'three.GC'),
                 aes(y = GC, x = after_stat(density), fill = snp),
                 binwidth = binwidth, alpha = 0.7) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  geom_hline(yintercept = length(which(mg1655 %in% c('G', 'C'))) / length(mg1655),
             linetype = 'dashed', linewidth = 0.3) +
  geom_richtext(data = stats,
                aes(x = ifelse(context == 'five.GC', -5, 5), label = label),
                y = 0.2, label.padding = unit(0, 'in'), label.size = 0, hjust = 0.5, size = 3) +
  geom_richtext(data = stats, aes(x = -0.3, y = 0.82, label = count),
                label.size = 0, hjust = 0.5, size = 3.5, label.padding = unit(0, 'in')) +
  scale_y_continuous(breaks = c(0.25, 0.5, 0.75)) +
  scale_x_continuous(breaks = c(-5, 0, 5), labels = c(5, 0, 5)) +
  scale_fill_manual(values = pal_spec, labels = function(x) sub('>', '→', x)) +
  facet_grid(cols = vars(group), rows = vars(snp),
             labeller = as_labeller(function(x) sub('>', '→', x))) +
  labs(x = 'Probability Density (left is 5\' context, right is 3\' context)',
       y = 'GC-content of ±100 bp context',
       fill = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(family = 'Aptos'),
        legend.text = element_text(family = 'Aptos')) +
  guides(fill = guide_legend(position = 'top', nrow = 1))


ggsave('figS10.png', width = 7*1.1, height = 9*1.1)





#### Figure S9 (Bias towards GC-rich regions for G:C→C:G is not exclusive to GC3+ mutations) ####


plot_data = my_data %>%
  filter(snp == 'GC>CG') %>%
  mutate(run.length = attr(regexpr("^C*", substring(motif, 12, 21), perl = TRUE), "match.length"),
         hotspot = ifelse(run.length > 2, TRUE, FALSE)) %>%
  mutate(context = list(contexts)) %>%
  unnest(context) %>%
  rowwise() %>%
  mutate(GC = GC_context(mg1655, pos, context, ref)) %>%
  ungroup() %>%
  group_by(group, snp, hotspot, context) %>%
  summarise(GC_mean = mean(GC), GC_sd = sd(GC), n = n(), .groups = "drop")



plots = list()
for (j in c(TRUE, FALSE)) {
  for (i in 1:4) {
    plot_data2 = plot_data %>% filter(group == groups[i] & hotspot == j)

    a = ggplot(plot_data2 %>% filter(context < 0),
               aes(x = context, y = GC_mean)) +
      geom_hline(yintercept = length(which(mg1655 %in% c('G', 'C'))) / length(mg1655), color = 'black') +
      geom_ribbon(aes(ymin = GC_mean - GC_sd / sqrt(n), ymax = GC_mean + GC_sd / sqrt(n)),
                  fill = pal_spec[6], linetype = 0, alpha = 0.3) +
      geom_smooth(method = 'loess', se = F, color = pal_spec[6]) +
      scale_x_log10(trans = pseudo_log_trans(base = 10), breaks= c(-1000, -100, -10, -1)) +
      coord_cartesian(ylim = c(0.425, 0.575), xlim = c(-1000, -1)) +
      labs(y = NULL, x = NULL) +
      theme_bw() +
      theme(panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.text.y.left = if (i == 1) element_text() else element_blank(),
            axis.line.x = element_line(linewidth = 0.3),
            axis.line.y.left = element_line(linewidth = 0.3),
            axis.line.y.right = element_blank())
    
    b = ggplot(plot_data2 %>% filter(context > 0),
               aes(x = context, y = GC_mean)) +
      geom_hline(yintercept = length(which(mg1655 %in% c('G', 'C'))) / length(mg1655), color = 'black') +
      geom_ribbon(aes(ymin = GC_mean - GC_sd / sqrt(n), ymax = GC_mean + GC_sd / sqrt(n)),
                  fill = pal_spec[6], linetype = 0, alpha = 0.3) +
      geom_smooth(method = 'loess', se = F, color = pal_spec[6]) +
      scale_fill_manual(values = pal_spec,
                        labels = c('A:T→G:C', 'G:C→A:T', 'A:T→C:G', 'G:C→T:A', 'A:T→T:A', 'G:C→C:G'),
                        guide = NULL) +
      scale_color_manual(values = pal_spec, guide = NULL) +
      scale_y_continuous(position = 'right') +
      scale_x_log10(trans = pseudo_log_trans(base = 10), breaks = c(1, 10, 100, 1000)) +
      coord_cartesian(ylim = c(0.425, 0.575), xlim = c(1, 1000)) +
      labs(y = NULL, x = NULL) +
      theme_bw() +
      theme(panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.text.y.right = if (i == 4) element_text() else element_blank(),
            axis.line.x = element_line(linewidth = 0.3),
            axis.line.y.left = element_blank(),
            axis.line.y.right = element_line(linewidth = 0.3))
    
    if (i == 1) {my_plot = plot_grid(a, b, rel_widths = c(1.12, 0.88))}
    else if (i == 4) {my_plot = plot_grid(a, b, rel_widths = c(0.88, 1.12))}
    else {my_plot = plot_grid(a, b)}
    
    if (j == 1) {
      my_plot = ggdraw() +
        draw_plot(my_plot, x = 0, y = 0, width = 1, height = 0.95) +
        draw_label(groups[i], x = ifelse(i == 1, 0.55, ifelse(i == 4, 0.45, 0.5)), y = 0.98, hjust = 0.5, size = 12)
    }
    
    plots %<>% append(list(my_plot))
  }
}

plot = plot_grid(plotlist = plots, ncol = 4,
                 rel_widths = c(1.13, 1, 1, 1.13),
                 rel_heights = c(1.05, 1))

plot = ggdraw() +
  draw_plot(plot, x = 0, y = 0.025, width = 1, height = 0.975) +
  draw_label('Starting position of 20 bp sliding context window', x = 0.5, y = 0.015, size = 12)

ggdraw() +
  draw_plot(plot, x = 0.025, y = 0, width = 0.95, height = 0.975) + 
  draw_label('Mean GC content of context window', x = 1 - 0.015, y = 0.5, angle = 270, size = 12)


ggsave('figS9.svg', width = 10, height = 5)

my_data %>%
  filter(snp == 'GC>CG') %>%
  mutate(run.length = attr(regexpr("^C*", substring(motif, 12, 21), perl = TRUE), "match.length"),
         hotspot = ifelse(run.length > 2, TRUE, FALSE)) %>%
  group_by(group, hotspot) %>%
  summarise(count = n())


#### Supp. Figure 10 (Distribution of GC-content values for the ±1 to ±100 bp context regions) ####
#### 3' regions of AT>CG or GC>CG by strand ####

focal.base = 'G'
## calculates observed numbers of Cs at different positions
## downstream of every G in the genome

# adds 10 bases either side of genome to circular searching over the genome ends
genome = as.character(mg1655) %>% paste0(collapse = '')
genome = paste0(substr(genome, nchar(genome) - 9, nchar(genome)), genome, substr(genome, 1, 10))

# identifies postions of Cs to 3' of Gs across genome
swap <- c(A = "T", T = "A", G = "C", C = "G")
gvec_padded <- strsplit(genome, "")[[1]]
L <- length(gvec_padded) - 20  # original genome length
pos <- 1:L
center <- gvec_padded[pos + 10]

# Strand classification
strand <- ifelse(pos >= 1640202 & pos < 3925696, "leading", "lagging")

# --- G positions ---
g_pos <- pos[center == focal.base]
g_motif <- substring(genome, g_pos + 11, g_pos + 20)
g_strand <- strand[center == focal.base]

# --- C positions ---
c_pos <- pos[center == swap[focal.base]]
c_raw <- substring(genome, c_pos + 0, c_pos + 9)

c_motif <- vapply(
  strsplit(c_raw, ""),
  \(x) paste(rev(swap[x]), collapse = ""),
  character(1)
)

c_strand <- ifelse(c_pos >= 1640202 & c_pos < 3925696, "lagging", "leading")

# --- Combine, mask non-Cs, and count ---
run_length <- attr(regexpr("^C*", c(g_motif, c_motif), perl = TRUE), "match.length")

genome.Cs <- tibble(
  run_length = run_length,
  strand = c(g_strand, c_strand)
) %>%
  dplyr::count(run_length, strand, name = "expected")


## creates dataframe looking at 3' Cs up to 10 nucleotides away
my_data %<>% filter(snp == ifelse(focal.base == 'G', 'GC>CG', 'AT>CG'))

nucleotide = 'C'
start = -1
end = 12

my_data %<>%
  rowwise() %>%
  mutate(value = case_when(
    ref %in% c('A', 'G') ~ get(genome)[((pos + start) : (pos + end) - 1) %% length(get(genome)) +  1] %>% paste(collapse = ''),
    T ~ rev(rev.comp(get(genome)[((pos - start) : (pos - end) - 1) %% length(get(genome)) +  1])) %>% paste(collapse = '')
  )) %>%
  ungroup()


my_data %<>% mutate(start = substr(value, 1, 1),
               value = substr(value, 3, 10))


my_data %<>% mutate(value = gsub(paste0('[^', nucleotide, ']'), "-", value))



plot_data = my_data %>%
  mutate(strand = case_when(
    (ref %in% c('A', 'G') & replichore == 'left') | (ref %in% c('T', 'C') & replichore == 'right') ~ 'leading',
    TRUE ~ 'lagging'
  )) %>%
  group_by(group,
           strand,
           snp,
           value) %>%
  summarise(count = n(), .groups = 'drop') %>%
  mutate(run_length = attr(regexpr("^C*", value), "match.length"))


# summarise counts at different run lengths
plot_data %<>%
  group_by(group, strand, snp, run_length) %>%
  summarise(count = sum(count), .groups = 'drop')


# # groups larger run lengths if wanted ONLY DO IF WANTED
# num = 5
# genome.Cs %<>%
#   filter(run_length < num) %>%
#   bind_rows(filter(genome.Cs, run_length >= num) %>%
#               group_by(strand) %>%
#               summarise(expected = sum(expected)) %>%
#               bind_cols(run_length = num))
# 
# plot_data %<>%
#   filter(run_length < num) %>%
#   bind_rows(filter(plot_data, run_length >= num) %>%
#               group_by(group, strand, snp) %>%
#               summarise(count = sum(count)) %>%
#               bind_cols(run_length = num))


# calculates rate from number of sites in genome and generations
plot_data %<>%
  left_join(genome.Cs) %>%
  left_join(rates_summary %>% select(group, generations)) %>%
  mutate(expected = as.numeric(expected),
         rate = count / expected / generations,
         group = factor(group, levels = groups),
         strand = factor(strand, levels = c('leading', 'lagging')))


# calculates rate relative to having 0 Cs ONLY DO IF WANTED
zero_data = plot_data %>%
  filter(run_length == 0)

plot_data = zero_data %>%
  select(group, strand, snp, rate) %>%
  rename(base_rate = rate) %>%
  right_join(plot_data) %>%
  filter(run_length != 0) %>%
  mutate(rel_rate = rate / base_rate)

# # fills dataframe with missing combinations, gives them 0 value ONLY DO IF WANTED
# plot_data %<>%
#   complete(
#     group, start, snp, run_length,
#     fill = list(count = NA, expected = NA, generations = NA, rate = NA, rel_rate = NA, base_rate = NA))

# calculates n= labels
labels = plot_data %>%
  group_by(group, snp, run_length) %>%
  summarise(count = sum(count), height = max(rel_rate) * 8, .groups = 'drop')


# # # fills dataframe with missing combinations, gives them 0 value ONLY DO IF WANTED
# plot_data %<>%
#   complete(
#     group, strand, snp, run_length,
#     fill = list(count = NA, expected = NA, generations = NA, rate = NA, rel_rate = NA, base_rate = NA))



library(ggtext)
b = ggplot(plot_data, aes(x = run_length, y = rel_rate)) +
  geom_bar(aes(fill = strand), stat = "identity", color = 'black',
           position = position_dodge2(preserve = "single")) +
  geom_label(data = labels, aes(x = run_length, y = height, label = count),
             size = 3, label.size = 0, fill = NA) +
  facet_wrap(~group,
             scales = 'free_x',
             ncol = 1) +
  # facet_wrap(~snp + group, ncol = 4, scales = 'free_y') +
  scale_fill_manual(values = c('white', 'darkgrey'), labels = c('Leading strand', 'Lagging strand')) +
  scale_x_continuous(breaks = c(1:8)) +
  scale_y_log10(breaks = c(1e0, 1e1, 1e2, 1e3, 1e4, 1e5),
                labels = c('1', paste0(rep(10, 5), '<sup>', c(1:5), '</sup>'))) +
  coord_cartesian(xlim = c(0.8, max(plot_data$run_length) + 0.2)) +
  geom_hline(yintercept = 1) +
  labs(x = 'Run length of Cs 3\' of mutation site', y = 'Fold increase in mutation rate relative to zero 3\' Cs') +
  theme_bw() +
  # guides(fill = guide_legend(position = 'top', title = 'Leading strand template')) +
  guides(fill = guide_legend(position = 'top', title = paste0(focal.base, ' templates the:'))) +
  theme(axis.text.y.left = element_markdown(),
        axis.title.x = element_text(hjust = 0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        # strip.text = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        legend.justification = 'left')


zero_data %<>%
  mutate(group = gsub('proofreading', '', group)) %>%
  mutate(group = gsub('MMR', '', group)) %>%
  mutate(group = factor(group, levels = c('(+) (+)', '(+) (-)', '(-) (+)', '(-) (-)')))

labels = zero_data %>%
  group_by(group, snp, run_length) %>%
  summarise(count = sum(count), height = max(rate) * 1.1, .groups = 'drop')

limits = tibble(group = levels(zero_data$group) %>% as_factor(),
                value = c(7e-11, 2e-10, 2e-8, 7e-9))

a = ggplot(zero_data, aes(x = run_length, y = rate)) +
  geom_bar(aes(fill = strand), position = 'dodge', stat = "identity", color = 'black') +
  geom_label(data = labels, aes(x = run_length, y = height,
                                label = paste0('*n = ', count, '*')),
             size = 3, label.size = 0, label.padding = unit(0, 'pt')) +
  geom_hline(data = limits, aes(yintercept = value), size = 0) +
  facet_wrap(~group,
             scales = 'free',
             ncol = 1) +
  # facet_wrap(~snp + group, ncol = 4, scales = 'free_y') +
  scale_fill_manual(values = c('white', 'darkgrey'), labels = c('Leading strand', 'Lagging strand')) +
  scale_x_continuous(breaks = c(0:6), labels = c(paste(0:5), '6+')) +
  scale_y_continuous(labels = function(x) (sci_label(x, digits = 0))) +
  # scale_y_continuous(position = 'right') +
  labs(x = '', y = 'Mutation rate (per site per generation)') +
  theme_bw() +
  guides(fill = guide_legend(position = 'top', title = NULL, override.aes = list(color = NA, fill = NA))) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background = element_blank(), 
        legend.text = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_markdown(),
        axis.text.y.left = element_markdown(),
        axis.ticks.y.left = element_blank())


p2 <- ggdraw() +
  draw_plot(a, x = 0, y = 0, width = 0.3, height = 1) +
  draw_plot(b, x = 0.3, y = 0, width = 0.7, height = 1)


plot_grid(p1, p2, labels = 'AUTO')

ggsave('plot.svg', height = 10, width = 9)

#### mutation rate / spectra by strain (supplementary) ####

rel_freq = tibble()
for (i in 1:3) {
  df_group = filter(my_data, group == groups[i])
  
  strains = unique(df_group$strain)
  
  x = tibble()
  for (j in 1:length(strains)) {
    
    plot_data = filter(my_data, strain == strains[j]) %>%
      group_by(strain, snp) %>%
      summarise(count = n(), .groups = 'drop')
    
    plot_data %<>%
      group_by(strain) %>%
      summarise(total = sum(count)) %>%
      right_join(plot_data, join_by(strain)) %>%
      mutate(freq = count/total) %>%
      select(strain, snp, freq)
    
    x %<>% bind_rows(plot_data)
  }
  
  rel_freq %<>% bind_rows(x)
}

strain_rates <- rates %>%
  mutate(
    rate = BPS / generations / length(mg1655),
    # se = CI * 1e-10 / 1.96,
    w  = generations,
    rate_w = rate * w
  ) %>%
  group_by(group, strain) %>%
  summarise(
    rate = sum(rate_w) / sum(w),
    # se_w = sqrt(sum((w^2) * (se^2)) / (sum(w)^2)),
    # CI   = 1.96 * se_w,
    .groups = "drop"
  )


rel_freq %<>%
  left_join(strain_rates, by = 'strain', relationship = 'many-to-many') %>%
  mutate(rel_rate = rate * freq
         # , CI = CI * freq
         )


x = c('WT', 'uvrA', 'nfi', 'alkA+tagA', 'ada+ogt', 'dinB+<br>umuDC', 'dinB+<br>umuDC+<br>polB',
      'mutS', 'mutL', 'mutH', 'mutS+<br>mutL', 'mutS+<br>mutL+<br>mutH', 'mutL+<br>dinB+<br>umuDC', 'mutS+<br>mfd', 'mutL+<br>mfd',
      'mutD5', 'mutD5+dinB', 'mutD5+dinB+umuDC')
names(x) = c('WT', 'uvrA', 'nfi', 'alkA_tagA', 'ada_ogt', 'umuDC_dinB', 'umuDC_dinB_polB',
             'mutS', 'mutL', 'mutH', 'mutS_mutL', 'mutS_mutL_mutH', 'mutL_umuDC_dinB', 'mutS_mfd', 'mutL_mfd',
             'D5_WT', 'D5_dinB', 'D5_dinB_umuDC')

rel_freq %<>%
  mutate(group = factor(group, levels = groups),
         strain = factor(strain, levels = names(x)))


ggplot(rel_freq) +
  geom_bar(aes(x = strain, y = rel_rate, fill = snp),
           stat = 'identity', position = 'dodge') +
  geom_errorbar(aes(x = strain, ymin = rel_rate - CI, ymax = rel_rate + CI, group = snp),
                stat = 'identity', position = 'dodge') +
  scale_x_discrete(labels = x) +
  scale_y_continuous(labels = function(x) sci_label(x, digits = 1)) +
  scale_fill_manual(values = pal_spec,
                    labels = function(x) sub(">", "→", x)) + 
  scale_color_manual(values = pal_spec) + 
  facet_wrap(~group, scales = 'free', ncol = 1) +
  labs(x = NULL, y = 'Mutation rate per base pair per generation') +
  theme_bw() +
  guides(fill = guide_legend(position = 'top', title = NULL, nrow = 1)) +
  theme(strip.background = element_blank(),
        axis.text = element_markdown(),
        legend.text = element_text(family = 'Aptos', margin = margin(r = 6, l = 3),
                                   size = 11))

ggsave('plot.png', height = 8, width = 6)




#### filtering out specific motifs ####



my_data %>%
  group_by(motif) %>%
  summarise(count = n()) %>%
  arrange(desc(count))


my_data %<>% filter(substr(motif, 10, 13) == 'GATC')

my_data %<>% filter(substr(motif, 8, 12) %in% c('CCTGG', 'CCAGG'))




test = my_data %>%
  filter(snp == 'GC>CG') %>%
  mutate(C.run = pmax(attr(regexpr("^C+", substr(motif, 12, 21)), "match.length"), 0))


filter(test, group != groups[1] & C.run > 2) %>% nrow() %>% divide_by(
  (filter(test, group != groups[1]) %>% nrow())
)



my_data %<>%
  filter(! (substr(motif, 8, 12) %in% c('CCTGG', 'CCAGG')
            # & snp == 'GC>AT'
  )) %>%
  filter(! (substr(motif, 10, 13) == 'GATC'
            # & snp == 'AT>TA'
  )) %>%
  mutate(run.3 = substr(motif, 12, 21),
         run.5 = substr(motif, 1, 10) %>% stringi::stri_reverse()) %>%
  mutate(run_length.3 = attr(regexpr("^(.)(\\1*)", run.3, perl = TRUE), "match.length"),
         run_length.5 = attr(regexpr("^(.)(\\1*)", run.5, perl = TRUE), "match.length")) %>%
  filter(! (run_length.3 > 2 | run_length.5 > 2
            # & snp %in% c('AT>CG', 'GC>CG')
  ))




my_data %<>% anti_join(exclude)

my_data %>%
  # filter(substr(motif, 10, 13) == 'GATC') %>%
  filter(substr(motif, 8, 12) %in% c('CCTGG', 'CCAGG')) %>%
  select(snp) %>% table()





#### new three prime c ####


nucleotide = 'C'



plot_data = my_data %>% mutate(start = substr(motif, 10, 10),
               value = substr(motif, 12, 21))


plot_data %<>% mutate(value = gsub(paste0('[^', nucleotide, ']'), "-", value))

# count run of Cs
plot_data %<>%
  mutate(run_length = attr(regexpr("^C*", value), "match.length")) %>%
  mutate(end = substr(motif, 12 + run_length, 13 + run_length))


plot_data %<>%
  group_by(group,
           start,
           end,
           snp,
           run_length) %>%
  summarise(count = n(), .groups = 'drop')


plot_data %<>%
  filter(snp %in% c('AT>CG', 'GC>CG')) %>%
  filter(run_length > 2) %>%
  group_by(start, end, snp, group) %>%
  summarise(count = sum(count))

plot_data %<>%
  left_join(plot_data %>%
              group_by(group, snp) %>%
              summarise(total = sum(count))) %>%
  mutate(prop = count / total)


ggplot(plot_data, aes(y = end, x = start)) +
  geom_point(aes(size = count, color = count)) +
  labs(x = '5\' nucleotide', y = '3\' dinucleotide') +
  # facet_wrap(~snp) +
  facet_grid(cols = vars(snp), rows = vars(group)) +
  theme_bw() +
  theme(strip.background = element_blank())



# summarise counts at different run lengths
plot_data %<>%
  group_by(group, start, snp, run_length) %>%
  summarise(count = sum(count), .groups = 'drop')

#### model ####

my_data = df_save


model_df = my_data %>% mutate(orientation = case_when(
  replichore == 'left' & ref %in% c('A', 'G') | replichore == 'right' & ref %in% c('T', 'C') ~ 'leading',
  TRUE ~ 'lagging'
))


model_df %<>%
  select(group, snp, base, value, orientation, total) %>%
  rename(position = base, nucleotide = value, strand = orientation) %>%
  filter(position != 0)


# of mutations with that -6+6, # in genome with taht -6 to +6, -6 to +6 letters, strand, BPS type


model_data = my_data %>%
  mutate(motif = substr(motif, 5, 17)) %>%
  group_by(group, snp, strand, motif) %>%
  summarise(count = n(), .groups = 'drop')

for (i in c(1:6, 8:13)) {
  model_data %<>% mutate(!!paste0("base_", i) := substr(motif, i, i))
}

model_data %<>% select(-motif)


model1 = glm(count ~ group*snp*strand*base_1+
      group*snp*strand*base_2+
      group*snp*strand*base_3+
      group*snp*strand*base_4+
      group*snp*strand*base_5+
      group*snp*strand*base_6+
      group*snp*strand*base_8+
      group*snp*strand*base_9+
      group*snp*strand*base_10+
      group*snp*strand*base_11+
      group*snp*strand*base_12+
      group*snp*strand*base_13,
    data = model_data,
    family = poisson(link = "log")
)


library(randomForest)
randomForest::randomForest(count ~ ., importance = TRUE,
                           data = model_data %>% select(-motif) %>% filter(!snp %in% c('AT>GC', 'GC>AT')),
                           ntree = 5000)

randomForest(y = pull(model_data, count), x = #remove count from model_df ., importance = TRUE,
                           data = model_data %>% select(-motif) %>% filter(!snp %in% c('AT>GC', 'GC>AT')),
                           ntree = 5000, proximity = TRUE, importance = TRUE, keep.inbag = TRUE, keep.forest = TRUE)

#### Rok's di/tri nucleotide stuff ####

wrap.genome = as.character(mg1655) %>% paste0(collapse = '')
wrap.genome = paste0(substr(wrap.genome, gen.length - 9, gen.length), wrap.genome, substr(genome, 1, 10))

genome.motifs = tibble(
  pos = 1:length(mg1655),
  ref = as.character(mg1655)
  ) %>%
  mutate(motif = substring(wrap.genome, pos, pos + 20)) %>%
  mutate(motif = ifelse(ref %in% c('T', 'C'), Vectorize(rev.comp)(motif), motif)) %>%
  mutate(focus = substr(motif, 11, 11))


bases = c(1:5) # dinucs
bases = c(1:3) # trinucs

bases = c(rev(bases)*-1, bases)

exp_dinucs = genome.motifs %>%
  mutate(base = list(bases)) %>%
  unnest(base) %>%
  mutate(value = case_when(
    base < 0 ~ substr(motif, 11 + base - 2, 11 + base),
    TRUE ~ substr(motif, 11 + base, 11 + base + 2)
  ))

exp_dinucs %<>%
  group_by(focus, base, value) %>%
  summarise(count = n(), .groups = 'drop') %>%
  mutate(perc = count / gen.length) 


ggplot(exp_dinucs) +
  geom_tile(aes(x = as_factor(base), y = value, fill = count)) +
  geom_vline(xintercept = 5.5, linewidth = 2, color = 'grey') +
  facet_wrap(~focus)

di_nucs = my_data %>%
  mutate(base = list(bases)) %>%
  unnest(base) %>%
  mutate(value = case_when(
    base < 0 ~ substr(motif, 11 + base - 2, 11 + base),
    TRUE ~ substr(motif, 11 + base, 11 + base + 2)
  ))

di_nucs %<>%
  group_by(group, snp, base, value) %>%
  summarise(count = n(), .groups = 'drop')

di_nucs %<>%
  left_join(di_nucs %>% group_by(group, snp) %>%
              summarise(total = sum(count), .groups = 'drop')) %>%
  mutate(perc = count / total)


di_nucs %<>%
  mutate(focus = substr(snp, 1, 1)) %>%
  left_join(exp_dinucs %>% rename(exp.perc = perc) %>% select(-count))

ggplot(di_nucs) +
  geom_tile(aes(x = as_factor(base), y = value, fill = perc - exp.perc)) +
  geom_vline(xintercept = 3.5, linewidth = 1, color = 'grey') +
  facet_grid(cols = vars(group), rows = vars(snp)) +
  theme(axis.text.y = element_text(size = 6))









#### rowans hotspot stuff
plot_data = my_data %>%
  filter(group == groups[4]) %>%
  mutate(hotspot = case_when(
    substr(motif, 9, 13) == 'CCAGG' | substr(motif, 9, 13) == 'CCTGG' ~ 'CCWGG',
    substr(motif, 10, 13) == 'GATC' ~ 'GATC',
    TRUE ~ NA_character_
  )) %>%
  mutate(new_motif = paste0(substr(motif, 1, 10), substr(snp, 4, 4), substr(motif, 12, 21))) %>%
  mutate(new_hotspot = case_when(
    substr(new_motif, 9, 13) == 'CCAGG' | substr(new_motif, 9, 13) == 'CCTGG' ~ 'CCWGG',
    substr(new_motif, 10, 13) == 'GATC' ~ 'GATC',
    TRUE ~ NA_character_
  )) %>%
  select(-new_motif)





# adds 10 bases either side of genome to circular searching over the genome ends
wrap.genome = as.character(mg1655) %>% paste0(collapse = '')
wrap.genome = paste0(substr(wrap.genome, gen.length - 9, gen.length), wrap.genome, substr(genome, 1, 10))

genome.motifs = tibble(
  pos = 1:length(mg1655),
  ref = as.character(mg1655)
) %>%
  mutate(motif = substring(wrap.genome, pos, pos + 20)) %>%
  mutate(motif = ifelse(ref %in% c('T', 'C'), Vectorize(rev.comp)(motif), motif)) %>%
  mutate(focus = substr(motif, 11, 11))


genome.hotspots = genome.motifs %>%
  mutate(dam.og = str_count(motif, 'GATC'),
         dcm.og = str_count(motif, 'CCAGG') + str_count(motif, 'CCTGG'))


#you are rambunctious fo sure my little ggplot2mistresslittle baby ggplot 




#bug
