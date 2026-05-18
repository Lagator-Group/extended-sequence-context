### PREAMBLE ####

library(seqinr)
library(stringi)
library(tidyverse)
library(magrittr)
library(cowplot)
library(scales)
library(ggtext)


# helpful vectors and color palettes
groups = c('proofreading(+) MMR(+)', 'proofreading(+) MMR(-)', 'proofreading(-) MMR(+)', 'proofreading(-) MMR(-)') %>% as_factor()
spectrum = c('AT>GC', 'GC>AT', 'AT>CG', 'GC>TA', 'AT>TA', 'GC>CG') %>% as_factor()

pal_spec = c('#7E3F8F', '#EE81EC', '#6E9BF8', '#53B74C', '#D9893B', '#B3A033')
names(pal_spec) = spectrum

bases = c('A', 'T', 'G', 'C')
pal_bases = c('#B2CB94', '#DDA97C', '#7DA1D0', '#F4DB85')
names(pal_bases) = bases

# imports collated dataset
my_data = read_csv('supplementary_data_1.csv',
                    skip = 1) %>%
  mutate(group = factor(group, levels = groups),
         BPS = factor(BPS, levels = spectrum))

# imports genome sequence
mg1655 = read.fasta('mg1655_NC_000913.3.fasta') %>%
  unlist() %>% unname() %>% toupper() %>% factor(levels = c('A', 'T', 'G', 'C'))

gen.length = length(mg1655)

# this data is only in SI Appendix-Table S1 (i'm sorry) you can copy it from there 
sum_stats = read_csv('summary_stats.csv')

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
sci_label <- function(x, digits = 1, drop_one = TRUE) {
  if (all(is.na(x))) return(NA_character_)
  
  sapply(x, function(val) {
    if (is.na(val)) return(NA_character_)
    if (val == 0) return("")
    
    sci <- format(val, scientific = TRUE)
    
    if (grepl("e", sci)) {
      parts <- strsplit(sci, "e")[[1]]
      base  <- as.numeric(parts[1])
      exp   <- as.integer(parts[2])
      
      base_fmt <- formatC(base, format = "f", digits = digits)
      
      if (drop_one && abs(base - 1) < 10^(-digits)) {
        paste0("10<sup>", exp, "</sup>")
      } else {
        paste0(base_fmt, "×10<sup>", exp, "</sup>")
      }
    } else {
      formatC(val, format = "f", digits = digits)
    }
  })
}

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
  mutate(focus = substr(motif, 11, 11)) %>%
  mutate(ref = factor(ref, levels = bases))



#### Table 1 (collated data and mutational spectrum) ####
rates_groups <- my_data %>%
  group_by(group, strain, experiment, BPS) %>%
  summarise(count = n(), .groups = "drop") %>%
  left_join(sum_stats %>% select(experiment, generations)) %>%
  group_by(group, BPS) %>%
  summarise(
    count = sum(count),
    generations = sum(generations),
    mean_rate = count / generations / gen.length) %>%
  mutate(
    lower = 0.5 * qchisq(0.025, 2 * count) / generations / gen.length,
    upper = 0.5 * qchisq(0.975, 2 * (count + 1)) / generations / gen.length,
  )


rates_summary = sum_stats %>%
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
      geom_bar(aes(x = BPS, y = mean_rate, fill = BPS),
               stat = 'identity', position = 'dodge') +
      geom_errorbar(aes(x = BPS, ymax = upper,
                        ymin = lower),
                    width = 0.6) +
      # geom_label(data = stats, aes(x = BPS, y = ymax + (0.1 * max), label = ratio), size = 2,
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


ggsave('Table 1.png', width = 9, height = 7)

#### Figure S1 (Mutational spectrum of strains in each repair group) ####

plot_data <- my_data %>%
  group_by(group, strain, experiment, BPS) %>%
  summarise(count = n(), .groups = "drop") %>%
  left_join(sum_stats %>% select(experiment, generations)) %>%
  mutate(experiment = as.numeric(as_factor(experiment))) %>%
  group_by(group, strain, BPS) %>%
  summarise(
    count = sum(count),
    generations = sum(generations),
    mean_rate = count / generations / gen.length,
    .groups = 'drop') %>%
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


ggplot(plot_data %>% filter(group != groups[4]), aes(x = strain, y = mean_rate, fill = BPS)) +
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


ggsave('figS1.png', width = 8, height = 8)


#### Table S1 (summary of MA experiments) ####

# gives the strains "fancy" names
strains = c('PFM2 (WT)', 'ΔuvrA', 'Δnfi', 'Δada+ogt', 'ΔalkA+tagA', 'ΔdinB+umuDC', 'ΔdinB+umuDC+polB',
            'ΔmutS', 'ΔmutL', 'ΔmutH', 'ΔmutSL', 'ΔmutSLH', 'ΔmutS+mfd', 'ΔmutL+mfd', 'ΔmutL+dinB+umuDC',
            'dnaQ-T15I', 'dnaQ-T15I (ΔdinB)', 'dnaQ-T15I (ΔdinB+umuDC)',
            'dnaQ-T15I (ΔmutL)')

names(strains) = c('WT', 'uvrA', 'nfi', 'ada_ogt', 'alkA_tagA', 'umuDC_dinB', 'umuDC_dinB_polB',
                   'mutS', 'mutL', 'mutH', 'mutS_mutL', 'mutS_mutL_mutH', 'mutS_mfd', 'mutL_mfd', 'mutL_umuDC_dinB',
                   'D5_WT', 'D5_dinB', 'D5_dinB_umuDC',
                   'D5_mutL')
# write the table
my_data %>%
  group_by(group, source, strain, experiment) %>%
  summarise(BPS = n(), .groups = "drop") %>%
  left_join(sum_stats %>% select(experiment, generations, lines)) %>%
  group_by(group, source, strain) %>%
  arrange(group, source, strain, experiment, lines) %>%
  mutate(rate = BPS / generations / gen.length,
         lower = 0.5 * qchisq(0.025, 2 * BPS) / generations / gen.length,
         upper = 0.5 * qchisq(0.975, 2 * (BPS + 1)) / generations / gen.length) %>%
  write_csv('Table S1.csv')

#### Table S3 ####

plot_data <- my_data %>%
  group_by(group, strain, experiment, BPS) %>%
  summarise(count = n(), .groups = "drop") %>%
  left_join(sum_stats %>% select(experiment, generations)) %>%
  group_by(strain) %>%
  mutate(experiment = dense_rank(experiment)) %>%
  ungroup()

# which strains we have >1 experiment for
unique(plot_data %>% filter(experiment > 1) %>% pull(strain))

stats = tibble()

# calculates chi2 and cramers V across experiments for individual strains
for (str in unique(plot_data %>% filter(experiment > 1) %>% pull(strain))) {
  chi2 = plot_data %>%
    filter(strain == str) %>%
    select(experiment, BPS, count) %>%
    pivot_wider(names_from = BPS, values_from = count) %>%
    column_to_rownames('experiment')
  
  test = chisq.test(chi2)
  
  stats %<>% bind_rows(
    tibble(
      strain = str,
      exp = nrow(chi2),
      chi = test$statistic,
      df = test$parameter,
      pval = test$p.value,
      n = sum(chi2),
      k = min(nrow(chi2), ncol(chi2))
    ) %>%
      mutate(cramer = sqrt(chi / (n * (k-1))))
  )
}

# calculates chi2 and cramers V across experiments for each group
for (grp in groups) {
  chi2 = plot_data %>%
    filter(group == grp) %>%
    mutate(strain = paste0(strain, '_', experiment)) %>%
    group_by(strain, BPS) %>%
    summarise(count = sum(count), .groups = 'drop') %>%
    select(strain, BPS, count) %>%
    pivot_wider(names_from = BPS, values_from = count) %>%
    column_to_rownames('strain')
  
  test = chisq.test(chi2)
  
  stats %<>% bind_rows(
    tibble(
      strain = grp,
      exp = nrow(chi2),
      chi = test$statistic,
      df = test$parameter,
      pval = test$p.value,
      n = sum(chi2),
      k = min(nrow(chi2), ncol(chi2))
    ) %>%
      mutate(cramer = sqrt(chi / (n * (k-1))))
  )
}

write_csv(stats, 'Table S3.csv')


#### Figure S2 (Average genomic sequence context) ####

genome_context = genome.motifs %>%
  mutate(context = list(5:17)) %>%
  unnest(context) %>%
  mutate(value = substr(motif, context, context)) %>%
  group_by(focus, context, value) %>%
  summarise(genome.count = n(), .groups = "drop")

genome_context %<>%
  left_join(genome_context %>%
              group_by(focus, context) %>%
              summarise(total = sum(genome.count))) %>%
  mutate(perc = genome.count / total,
         context = context - 11,
         value = factor(value, levels = c('A', 'T', 'G', 'C'))) %>%
  select(-total)


ggplot(genome_context, aes(x = context, y = perc, fill = value)) +
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


ggsave('figS2.png', width = 6, height = 3)

#### Figure 1A (diagram) ####

ncol = 7
nrow = 8

plot_data = tibble(
  row = rep(1:nrow, ncol),
  col = rep(1:ncol, each = nrow)
)

# for custom base selection
my_bases = list(
  c('A', 'A', 'T', 'T', 'G', 'G', 'C', 'C'),
  c('A', 'A', 'T', 'A', 'G', 'G', 'C', 'G'),
  c('A', 'A', 'T', 'T', 'G', 'G', 'C', 'C'),
  rep('A', 8),
  c('A', 'C', 'T', 'C', 'G', 'C', 'C', 'C'),
  c('A', 'A', 'T', 'C', 'G', 'G', 'C', 'C'),
  c('A', 'A', 'T', 'T', 'G', 'G', 'C', 'C')
)

set.seed(11)
col_bases = c()
for (i in 1:length(my_bases)) {
  col_bases %<>% c(sample(my_bases[[i]], 8))
}

plot_data %<>% bind_cols(base = col_bases)

a = ggplot(plot_data, aes(x = col, y = row, fill = base)) +
  geom_tile(linewidth = 0) +
  geom_hline(yintercept = seq(0.5, nrow - 0.5), color = 'white', linewidth = 2) +
  # geom_rect(aes(xmin = 3.5, xmax = 4.5, ymin = row - 0.44, ymax = row + 0.44), fill = NA, color = 'black', linewidth = 1) +
  scale_fill_manual(values = pal_bases, guide = NULL) +
  scale_x_continuous(breaks = c(1:7), labels = c('-3', '-2', '-1', '0', '+1', '+2', '+3')) +
  labs(x = 'Context position') +
  theme_void() +
  theme(axis.text.x = element_text(margin = margin(t = -5), family = 'Aptos', size = 25),
        axis.title.x = element_text(size = 25))


b = plot_data %>%
  mutate(base = factor(base, levels = bases)) %>%
  arrange(col, base) %>%
  bind_cols(new_row = rep(1:nrow, ncol)) %>%
  ggplot(aes(x = col, y = new_row, fill = base, alpha = factor(col))) +
  geom_tile() +
  geom_vline(xintercept = seq(0.5, ncol - 0.5), color = 'white', linewidth = 2) +
  geom_hline(yintercept = seq(0.5, nrow + 0.5, length.out = 5)[-c(1, 5)], color = 'black', linewidth = 1) +
  scale_y_reverse(breaks = seq(0.5, nrow + 0.5, length.out = 5)[-c(1, 5)],
                  labels = c('25%', '50%', '75%'), sec.axis = dup_axis()) +
  scale_x_continuous(breaks = c(1:7), labels = c('-3', '-2', '-1', '0', '+1', '+2', '+3')) +
  scale_fill_manual(values = pal_bases, guide = NULL) +
  scale_alpha_manual(values = c(0.5, 1, 0.5, 1, 1, 1, 0.5), guide = NULL) +
  labs(x = 'Context position') +
  coord_cartesian(xlim = c(0.8, ncol+0.2)) +
  theme_void() +
  theme(axis.text.x = element_text(margin = margin(t = -5), family = 'Aptos', size = 25),
        axis.text.y.right = element_text(margin = margin(l = 7), family = 'Aptos', size = 25),
        axis.title.x = element_text(size = 25))

plot_grid(a, b)

ggsave('fig1a.png', width = 10, height = 4)

#### Figure 1B+C (sequence context nucleotide frequencies) ####

contexts = c(1:6)
contexts = c(rev(contexts)*-1, 0, contexts)

ctx_freq = my_data %>%
  mutate(context = list(contexts)) %>%
  unnest(context) %>%
  mutate(value = substr(motif, 11 + context, 11 + context))

ctx_freq$value %<>% factor(levels = c('A', 'T', 'G', 'C'))

plot_data <- ctx_freq %>%
  group_by(BPS, context, value, group) %>%
  summarise(count = n(), .groups = "drop") %>%
  complete(BPS, context, value, group, fill = list(n = 0)) # Fill missing combinations

plot_data %<>%
  mutate(focus = case_when(grepl('AT>', BPS) ~ 'A', TRUE ~ 'G')) %>% # adds col for focus base (A or G)
  left_join(genome_context %>% rename(gen.perc = perc)) %>% # adds info for the genome_context bias at each position
  mutate(weight.count = (0.25 / gen.perc) * count) # weights by genome structure

plot_data %<>%
  left_join(
    plot_data %>%
      group_by(BPS, context, group) %>%
      summarise(total = sum(weight.count))
  ) %>%
  mutate(weight.perc = weight.count/total,
         perc.sem = sqrt(weight.perc * (1 - weight.perc) / total)) %>%
  ungroup() %>%
  select(-total)


# chi squared test to determine if significant bias
stats = plot_data %>%
  group_by(BPS, context, focus, group) %>%
  summarise(n = sum(count))

p = tibble()
for (i in 1:nrow(stats)) {
  if (stats[[i,'context']] == 0) {p %<>% bind_rows(c(chi = NA, p = NA)) ; next}
  x = semi_join(ctx_freq, stats[i,], by = join_by(BPS, context, group))$value %>% table()
  
  expected = semi_join(genome_context, stats[i,], by = join_by(context, focus)) %>%
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


# stops context 0 from skewing results with maximum bias
plot_data$bias[plot_data$context == 0] = 
  max(plot_data$bias[plot_data$context != 0])
plot_data$weight.perc[plot_data$context == 0 & !is.na(plot_data$count)] = 1

plot_data %<>% mutate(nt.bias = abs(0.25 - weight.perc))

a = ggplot(plot_data, aes(x = context, y = weight.perc, fill = value, alpha = bias)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = pal_bases, name = 'Nucleotide') +
  scale_alpha_continuous(
    range = c(0, 1),
    name = 'Bias',
    limits = c(0, 0.5),
    breaks = c(0, 0.25, 0.5),
    labels = c('0', '0.25', '> 0.5')) +
  scale_y_continuous(breaks = c(0.25, 0.5, 0.75), labels = c('25%', '50%', '75%'),
                     sec.axis = dup_axis(name = 'Weighted frequency')) +
  scale_x_continuous(breaks = seq(-6, 6, 2), labels = c('-6', '-4', '-2', '0', '+2', '+4', '+6')) +
  geom_hline(yintercept = c(0.5), linewidth = 0.5) +
  geom_hline(yintercept = c(0.25, 0.75), linewidth = 0.25) +
  # facet_wrap(~BPS, labeller = as_labeller(labels)) +
  facet_grid(rows = vars(BPS), cols = vars(group), switch = 'y',
             labeller = as_labeller(function(x) sub(">", "→", x))
  ) +
  coord_cartesian(xlim = c(-6, 6)) +
  labs(y = NULL, x = "Context Position", fill = "Value") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y.left = element_text(family = 'Aptos'), 
        axis.text.y.left = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(size = 8)) + 
  guides(fill = guide_legend(position = 'top', order = 1),
         alpha = guide_legend(position = 'top', order = 2,
                              override.aes = list(color = 'black',
                                                  linewidth = 0.1,
                                                  fill = '#B6C999')))

fits = plot_data %>%
  filter(p.adj < 0.05) %>%
  group_by(group, BPS) %>%
  summarise(min = min(bias), .groups = 'drop') %>%
  left_join(plot_data %>%
              filter(p.adj >= 0.05) %>%
              group_by(group, BPS) %>%
              summarise(max = max(bias), .groups = 'drop')) %>%
  mutate(value = (min + max) / 2) %>%
  mutate(value = ifelse(is.na(value), min, value))




## plot for bias at each position, first estimates min significant bias
counts = my_data %>%
  group_by(group, BPS) %>%
  summarise(count = n(), .groups = 'drop') %>%
  arrange(group, BPS)


b = ggplot(filter(plot_data, context != 0), aes(x = context, y = bias, color = BPS)) +
  geom_hline(data = fits, aes(yintercept = value), linetype = 'dashed', alpha = 0.5) +
  geom_vline(xintercept = 0, color = 'grey') +
  geom_point(aes(alpha = p.adj < 0.05), size = 1.2) +
  geom_label(data = counts, aes(label = paste0('n=', count)),
             x = -6.4 , y = 1, color = 'black', fill = 'white', size = 3,
             linewidth = 0, hjust = 0) +
  # geom_label(data = counts, aes(label = loss),
  # x = 6.4 , y = 1, color = 'black', fill = 'white', size = 3,
  # label.size = 0, hjust = 1) +
  facet_grid(rows = vars(BPS), cols = vars(group), switch = 'y',
             labeller = as_labeller(function(x) sub(">", "→", x))) +
  scale_x_continuous(breaks = seq(-6, 6, 2), labels = c('-6', '-4', '-2', '0', '+2', '+4', '+6')) +
  scale_y_continuous(sec.axis = dup_axis(name = 'Deviation from expected frequencies'),
                     breaks = c(0, 0.5, 1)) +
  scale_color_manual(values = pal_spec, guide = NULL) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.1), guide = NULL) + 
  coord_cartesian(ylim = c(0, 1.2)) +
  labs(y = NULL, x = "Context Position", color = "BPS") +
  theme_bw() +
  theme(axis.text.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        strip.text.y.left = element_text(family = 'Aptos'),
        strip.background = element_rect(fill = NA, color = NA))


cowplot::plot_grid(plotlist = list(a, b), ncol = 1, labels = c('B', 'C'), align = 'v', rel_heights = c(1, 0.9))
ggsave('fig1bc.png', width = 7, height = 9)


#### Figure S3 (Dam and Dcm hotspots) ####
plot_data = my_data %>%
  mutate(type = case_when(
    substr(motif, 10, 13) == 'GATC' ~ 'G**A**TC',
    substr(motif, 8, 12) %in% c('CCTGG', 'CCAGG') ~ 'C**C**WGG',
    TRUE ~ 'other'
  ))


plot_data %<>%
  group_by(group, BPS, type) %>%
  summarise(count = n()) %>%
  left_join(plot_data %>%
              group_by(group, BPS) %>%
              summarise(total = n())) %>%
  mutate(prop = count/total,
         type = factor(type, levels = c('G**A**TC', 'C**C**WGG', 'other')))


expected = tibble(
  group = rep('Whole\ngenome', 3),
  type = c('G**A**TC', 'C**C**WGG', 'other'),
  gen.count = c(38240, 24090, gen.length - 24090 - 38240)
) %>% mutate(prop = gen.count / gen.length,
             type = factor(type, levels = c('G**A**TC', 'C**C**WGG', 'other')))


plot_data$BPS %<>% factor(levels = c('AT>GC', 'AT>CG', 'AT>TA', 'GC>AT', 'GC>TA', 'GC>CG'))


p1 = plot_data %>%
  mutate(type = factor(type, levels = c('G**A**TC', 'C**C**WGG', 'other'))) %>%
  ggplot(aes(x = BPS, y = prop, fill = type)) +
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
            aes(x = BPS, y = rate, fill = type)) +
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

ggsave('figS3.png', width = 8.5, height = 6)



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
  group_by(group, BPS, run.length) %>%
  summarise(count = n(), .groups = 'drop')


plot_data %<>%
  left_join(plot_data %>%
              group_by(group, BPS) %>%
              summarise(total = sum(count))) %>%
  mutate(perc = count / total,
         focus = substr(BPS, 1, 1))



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
      geom_point(data = plot_data %>% filter(BPS == mut),
                 aes(x = run.length, y = perc), color = 'black') +
      geom_text(data = plot_data %>% filter(BPS == mut),
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
           run.length >= 3 & BPS %in% c('AT>CG', 'GC>CG') ~ 'C3+',
           TRUE ~ NA_character_
         )) %>%
  filter(is.na(hotspot))


contexts = c(1:6)
contexts = c(rev(contexts)*-1, 0, contexts)

ctx_freq = filtered_data %>%
  mutate(context = list(contexts)) %>%
  unnest(context) %>%
  mutate(value = substr(motif, 11 + context, 11 + context))

ctx_freq$value %<>% factor(levels = c('A', 'T', 'G', 'C'))

plot_data <- ctx_freq %>%
  group_by(BPS, context, value, group) %>%
  summarise(count = n(), .groups = "drop") %>%
  complete(BPS, context, value, group, fill = list(n = 0)) # Fill missing combinations

plot_data %<>%
  mutate(focus = case_when(grepl('AT>', BPS) ~ 'A', TRUE ~ 'G')) %>% # adds col for focus base (A or G)
  left_join(genome_context) %>% # adds info for the genome_context bias at each position
  mutate(weight.count = (0.25 / perc) * count) # weights by genome structure

plot_data %<>%
  left_join(
    plot_data %>%
      group_by(BPS, context, group) %>%
      summarise(total = sum(weight.count))
  ) %>%
  mutate(perc = weight.count/total,
         perc.sem = sqrt(perc * (1 - perc) / total)) %>%
  ungroup() %>%
  select(-total)


# chi squared test to determine if significant bias
stats = plot_data %>%
  group_by(BPS, context, focus, group) %>%
  summarise(n = sum(count))

p = tibble()
for (i in 1:nrow(stats)) {
  if (stats[[i,'context']] == 0) {p %<>% bind_rows(c(chi = NA, p = NA)) ; next}
  x = semi_join(ctx_freq, stats[i,], by = join_by(BPS, context, group))$value %>% table()
  
  expected = semi_join(genome_context, stats[i,], by = join_by(context, focus)) %>%
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


# stops context 0 from skewing results with maximum bias
plot_data$bias[plot_data$context == 0] = 
  max(plot_data$bias[plot_data$context != 0])
plot_data$perc[plot_data$context == 0 & !is.na(plot_data$count)] = 1


a = ggplot(plot_data, aes(x = context, y = perc, fill = value, alpha = bias)) +
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
  # facet_wrap(~BPS, labeller = as_labeller(labels)) +
  facet_grid(rows = vars(BPS), cols = vars(group), switch = 'y',
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
  group_by(group, BPS) %>%
  summarise(og.count = n(), .groups = 'drop')

filt.counts = filtered_data %>%
  group_by(group, BPS) %>%
  summarise(count = n(), .groups = 'drop') %>%
  arrange(group, BPS)

filt.counts %<>%
  left_join(og.counts) %>%
  mutate(difference = og.count - count) %>%
  mutate(loss = round(difference / og.count, 2) * 100) %>%
  mutate(difference = paste0('-', difference)) %>%
  mutate(loss = paste0('-', loss, '%'))

# caluclate minimum pvalue threshold
fits = plot_data %>%
  filter(p.adj < 0.05) %>%
  group_by(group, BPS) %>%
  summarise(min = min(bias), .groups = 'drop') %>%
  left_join(plot_data %>%
              filter(p.adj >= 0.05) %>%
              group_by(group, BPS) %>%
              summarise(max = max(bias), .groups = 'drop')) %>%
  mutate(value = (min + max) / 2) %>%
  mutate(value = ifelse(is.na(value), min, value))


b = ggplot(filter(plot_data, context != 0), aes(x = context, y = bias, color = BPS)) +
  geom_hline(data = fits, aes(yintercept = value), linetype = 'dashed', alpha = 0.5) +
  geom_vline(xintercept = 0, color = 'grey') +
  geom_point(aes(alpha = p <= 0.05)) +
  geom_label(data = filt.counts, aes(label = paste0('n=', count)),
             x = -6.4 , y = 1, color = 'black', fill = 'white', size = 3,
             label.size = 0, hjust = 0) +
  geom_label(data = filt.counts, aes(label = loss),
             x = 6.4 , y = 1, color = 'black', fill = 'white', size = 3,
             label.size = 0, hjust = 1) +
  facet_grid(rows = vars(BPS), cols = vars(group), switch = 'y',
             labeller = as_labeller(function(x) sub(">", "→", x))) +
  scale_x_continuous(breaks = seq(-6, 6, 2), labels = c('-6', '-4', '-2', '0', '+2', '+4', '+6')) +
  scale_y_continuous(sec.axis = dup_axis(name = 'Absolute difference from expected frequencies'),
                     breaks = c(0, 0.5, 1)) +
  scale_color_manual(values = pal_spec, guide = NULL) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.1), guide = NULL) + 
  coord_cartesian(ylim = c(0, 1.2)) +
  labs(y = NULL, x = "Context Position", color = "BPS") +
  theme_bw() +
  theme(axis.text.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        strip.text.y.left = element_text(family = 'Aptos'),
        strip.background = element_rect(fill = NA, color = NA))


cowplot::plot_grid(plotlist = list(a, b), ncol = 1, labels = c('A', 'B'))
ggsave('figS5.png', width = 7, height = 8.5)


#### Figure 2 / Figure S6 (transient misalignment hotspots) ####

BPS.conv = rep(spectrum, each = 2)
names(BPS.conv) = c('AG', 'TC', 'GA', 'CT', 'AC', 'TG', 'GT', 'CA', 'AT', 'TA', 'GC', 'CG')


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
  
  # identifies which BPS would result from transient misalignment
  mutate(
    BPS.pur = case_when(
      cat.pur == 'nascent' ~ BPS.conv[paste0(focus, base.3)],
      cat.pur == 'template' ~ BPS.conv[paste0(focus, base.5)],
      TRUE ~ NA_character_
    ),
    BPS.pyr = case_when(
      cat.pyr == 'nascent' ~ BPS.conv[paste0(focus, base.5)],
      cat.pyr == 'template' ~ BPS.conv[paste0(focus, base.3)],
      TRUE ~ NA_character_
    )
  )

# diagnostic tibble
all.runs %>% select(motif, base.5, run.length.5, cat.pyr, BPS.pyr, base.3, run.length.3, cat.pur, BPS.pur)  %>%
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
             BPS.pur, BPS.pyr),
    names_to = c(".value", "template"),
    names_pattern = "(.*)\\.(pur|pyr)"
  )


# summarise counts
all.runs %<>%
  group_by(focus, template, run.length, cat, BPS) %>%
  summarise(gen.count = n(), .groups = 'drop') %>%
  mutate(BPS = factor(BPS, levels = spectrum),
         cat = factor(cat, levels = c('no', 'template', 'nascent')))

# diagnostic plot
ggplot(all.runs, aes(x = run.length, y = gen.count, shape = cat, color = BPS)) +
  # geom_bar(stat = 'identity', position = 'dodge') +
  geom_point() +
  scale_color_manual(values = pal_spec) +
  scale_y_log10() +
  facet_grid(cols = vars(template), rows = vars(BPS)) +
  theme_bw()


# combines runs >6 into 6+ category
max.length = 6
all.runs %<>%
  filter(run.length < max.length) %>%
  bind_rows(filter(all.runs, run.length >= max.length) %>%
              group_by(focus, template, cat, BPS) %>%
              summarise(gen.count = sum(gen.count), .groups = 'drop') %>%
              bind_cols(run.length = max.length))


### does the same for the mutation data

plot_data = my_data %>%
  mutate(focus = substr(BPS, 1, 1)) %>%
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
  
  # identifies which BPS would result from transient misalignment
  mutate(
    BPSTM.pur = case_when(
      cat.pur == 'nascent' ~ BPS.conv[paste0(focus, base.3)],
      cat.pur == 'template' ~ BPS.conv[paste0(focus, base.5)],
      TRUE ~ NA_character_
    ),
    BPSTM.pyr = case_when(
      cat.pyr == 'nascent' ~ BPS.conv[paste0(focus, base.5)],
      cat.pyr == 'template' ~ BPS.conv[paste0(focus, base.3)],
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
             BPSTM.pur, BPSTM.pyr),
    names_to = c(".value", "template"),
    names_pattern = "(.*)\\.(pur|pyr)"
  )


plot_data %<>%
  mutate(cat = case_when(
    cat == 'no' ~ cat,
    BPS ==  BPSTM ~ cat,
    TRUE ~ 'no'
  )) %>%
  group_by(group, BPS, focus, template, run.length, cat) %>%
  summarise(count = n(), .groups = 'drop') %>%
  left_join(rates_summary %>% select(group, generations))


# combines runs >6 into 6+ category
plot_data %<>% filter(run.length < max.length) %>%
  bind_rows(filter(plot_data, run.length >= max.length) %>%
              group_by(group, BPS, focus, template, cat, generations) %>%
              summarise(count = sum(count), .groups = 'drop') %>%
              bind_cols(run.length = max.length))


plot_data %<>%
  left_join(plot_data %>% group_by(group, BPS, template, run.length) %>%
              summarise(total = sum(count), .groups = 'drop')) %>%
  mutate(perc = count / total)

plot_data$cat %<>% factor(levels = c('nascent', 'template', 'no'))
plot_data$group %<>% factor(levels = groups)



plot_grid(
  ggplot(plot_data %>% filter(template == 'pur'), aes(x = run.length, y = perc)) +
    geom_bar(stat = 'identity', aes(fill = cat)) +
    geom_richtext(data = rates_groups, aes(label = mean_rate %>% sci_label(digits = 0, F)),
                  x = 0.5, y = 0, size = 2, hjust = 0, vjust = 0,
                  label.size = 0, label.padding = unit(1, 'pt')) + 
    scale_fill_manual(values = c('#56BCC2', '#E77D72', 'grey'),guide = NULL) +
    scale_x_continuous(breaks = c(1:6), labels = c(1, 2, 3, 4, 5, '6+')) +
    scale_y_continuous(breaks = c(0.25, 0.75), labels = c('25%', '75%')) +
    facet_grid(rows = vars(BPS), cols = vars(group),
               labeller = as_labeller(function(x) sub('proofreading', 'Pr', x) %>% sub('>', '→', .))) +
    labs(x = 'Run length', y = 'Proportion of mutations') +
    theme_bw() +
    theme(strip.background = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          strip.text.x = element_text(family = 'Aptos'),
          strip.text.y = element_text(color = 'white')),
  ggplot(plot_data  %>% filter(template == 'pyr'), aes(x = run.length, y = perc)) +
    geom_bar(stat = 'identity', aes(fill = cat)) +
    scale_fill_manual(values = c('#56BCC2', '#E77D72', 'grey'),guide = NULL) +
    scale_x_continuous(breaks = c(1:6), labels = c(1, 2, 3, 4, 5, '6+')) +
    scale_y_continuous(breaks = c(0.25, 0.75), labels = c('25%', '75%')) +
    facet_grid(rows = vars(BPS), cols = vars(group),
               labeller = as_labeller(function(x) sub('proofreading', 'Pr', x) %>% sub('>', '→', .))) +
    labs(x = 'Run length', y = ' ') +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text = element_text(family = 'Aptos'),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank())
)

ggsave('figS6b.png', width = 9, height = 6)




one_data = plot_data %>%
  filter(run.length <= 1) %>%
  group_by(group, BPS, focus, template, generations) %>% summarise(count = sum(count), .groups = 'drop') %>%
  left_join(all.runs %>%
              filter(run.length <= 1) %>%
              group_by(focus, template) %>% summarise(gen.count = sum(gen.count), .groups = 'drop') %>%
              select(focus, template, gen.count)) %>%
  mutate(one.rate = count / generations / gen.count)


plot_data = bind_rows(
  plot_data %>%
    filter(cat == 'no' & run.length > 1) %>%
    left_join(all.runs %>% filter(cat == 'no') %>% select(-BPS)),
  plot_data %>%
    filter(cat != 'no' & run.length > 1) %>%
    left_join(all.runs %>% filter(cat != 'no'))
) %>%
  mutate(rate = count / generations / gen.count)


plot_data2 = plot_data %>% bind_rows(one_data %>% rename(rate = one.rate) %>% mutate(run.length = 1, cat = 'no'))


plot_data2$cat %<>% factor(levels = c('nascent', 'template', 'no'))
plot_data2$group %<>% factor(levels = groups)

p2 = plot_grid(
  ggplot(plot_data2  %>% filter(template == 'pur'), aes(x = run.length, y = rate)) +
    geom_bar(stat = 'identity', aes(fill = cat)) +
    scale_fill_manual(values = c('#56BCC2', '#E77D72', 'grey'), guide = NULL) +
    scale_x_continuous(breaks = c(1:6), labels = c(1, 2, 3, 4, 5, '6+')) +
    scale_y_continuous(breaks = 1e-10, labels = c('   ')) +
    facet_wrap(~BPS+group, scales = 'free_y', ncol = 4) + 
    labs(x = 'Run length', y = NULL, title = 'Assuming purine templated mispair') +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_blank(),
          plot.background = element_blank(),
          panel.grid = element_blank(),
          # panel.border = element_blank(),
          axis.ticks = element_blank(),
    ),
  ggplot(plot_data2  %>% filter(template == 'pyr'), aes(x = run.length, y = rate)) +
    geom_bar(stat = 'identity', aes(fill = cat)) +
    scale_fill_manual(values = c('#56BCC2', '#E77D72', 'grey'), guide = NULL) +
    scale_x_continuous(breaks = c(1:6), labels = c(1, 2, 3, 4, 5, '6+')) +
    scale_y_continuous(breaks = 1e-10, labels = c('   ')) +
    facet_wrap(~BPS+group, scales = 'free_y', ncol = 4) + 
    labs(x = 'Run length', y = NULL, title = 'Assuming pyrimidine templated mispair') +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_blank(),
          plot.background = element_blank(),
          panel.grid = element_blank(),
          # panel.border = element_blank(),
          axis.ticks = element_blank(),
    )
)



plot_data3 = plot_data2 %>%
  group_by(group, BPS, template) %>%
  summarise(max.rate = max(rate)) %>%
  mutate(max.rate = max.rate* 1e10)

p1 = plot_grid(
  ggplot(plot_data3 %>% filter(template == 'pur'), aes(x = 1, y = max.rate)) +
    geom_bar(stat = 'identity', fill = 'black') +
    scale_y_log10(labels = function(x) sci_label(x * 1e-11), breaks = c(1, 100, 10000)) +
    facet_wrap(~BPS+group, ncol = 4) +
    labs(x = ' ', y = NULL, title = ' ') +
    coord_cartesian(xlim = c(0.5, 15), ylim = c(1, max(plot_data3$max.rate))) +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_blank(),
          axis.text.x = element_text(color = 'white'),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.border = element_blank()),
  ggplot(plot_data3 %>% filter(template == 'pyr'), aes(x = 1, y = max.rate)) +
    geom_bar(stat = 'identity', fill = 'black') +
    scale_y_log10(labels = function(x) sci_label(x * 1e-11), breaks = c(1, 100, 10000)) +
    facet_wrap(~BPS+group, ncol = 4) +
    labs(x = ' ', y = NULL, title = ' ') +
    coord_cartesian(xlim = c(0.5, 15), ylim = c(1, max(plot_data3$max.rate))) +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_blank(),
          axis.text.x = element_text(color = 'white'),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.border = element_blank())
)

axis = plot_grid(
  ggplot(plot_data3 %>% filter(template == 'pur'), aes(x = 1, y = max.rate)) +
    geom_bar(stat = 'identity', fill = 'white') +
    scale_y_log10(labels = function(x) sci_label(x * 1e-11), breaks = c(1, 100, 10000)) +
    facet_wrap(~BPS+group, ncol = 4) +
    labs(x = ' ', y = 'Mutation rate (with log scaling bar)', title = ' ') +
    coord_cartesian(xlim = c(0.5, 15), ylim = c(1, max(plot_data3$max.rate))) +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_blank(),
          axis.text.x = element_text(color = 'white'),
          axis.text.y = element_markdown(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank()),
  ggplot(plot_data3 %>% filter(template == 'pyr'), aes(x = 1, y = max.rate)) +
    geom_bar(stat = 'identity', fill = 'white') +
    scale_y_log10(labels = function(x) sci_label(x * 1e-11), breaks = c(1, 100, 10000)) +
    facet_wrap(~BPS+group, ncol = 4) +
    labs(x = ' ', y = NULL, title = ' ') +
    coord_cartesian(xlim = c(0.5, 15), ylim = c(1, max(plot_data3$max.rate))) +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_blank(),
          axis.text.x = element_text(color = 'white'),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank())
)


ggdraw() +
  draw_plot(axis) +
  draw_plot(p1, x = 0.06, width = 0.94) +
  draw_plot(p2, x = 0.06, width = 0.94)

ggsave('figS6a.png', width = 9, height = 6)



digits_left <- function(x) ifelse(x == 0 | is.na(x), 1, floor(log10(abs(x))) + 1)

digits_left(5.2e-10 * 1e10)

plot_data %<>% left_join(one_data %>% select(group, BPS, template, one.rate)) %>%
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
        facet_grid(rows = vars(BPS), cols = vars(group), switch = 'y',
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

ggsave('fig2bc.png', width = 9, height = 5)



#### Figure S7 (reoccurring invdividual mutations) ####

plot_data = my_data %>%
  group_by(pos, ref, BPS, group, motif) %>%
  summarise(count = n(), .groups = 'drop')

plot_data %>%
  group_by(count) %>%
  summarise(total = n()) %>%
  ggplot() +
  geom_point(aes(x = count, y = total), size = 3) +
  scale_x_continuous(breaks = c(1:9)) +
  scale_y_log10(labels = function(x) sci_label(x)) +
  labs(x = 'Times mutation was observed', y = 'Number of different mutations') +
  theme_bw() +
  theme(axis.text.y = element_markdown())



plot_data$pos %<>% as_factor()
plot_data$pos <- factor(plot_data$pos, levels = rev(levels(plot_data$pos)))

plot_data %<>%
  mutate(motif = case_when(
    substr(motif, 11, 11) == 'A' ~ Vectorize(rev.comp)(motif),
    TRUE ~ motif
  ))


plot_data %<>%
  filter(count > 3) %>%
  bind_cols(hotspot = c(
    'TC<sub>4</sub>',
    'GC<sub>6</sub>',
    'TC<sub>7</sub>',
    'unknown',
    'unknown',
    'TC<sub>5</sub>',
    'TC<sub>5</sub>',
    'GC<sub>7</sub>',
    'TC<sub>7</sub>',
    'TC<sub>8</sub>',
    'TC<sub>8</sub>',
    'TC<sub>8</sub>'
  ))


plot_data2 = plot_data %>%
  separate_longer_position(motif, width = 1) %>%
  group_by(pos) %>%
  mutate(motif_pos = row_number()) %>%
  ungroup() %>%
  mutate(motif = factor(motif, levels = names(pal_bases)))


a = ggplot(plot_data2,
           aes(x = motif_pos, y = pos, fill = motif)) +
  geom_tile(height = 0.9) +
  geom_rect(aes(xmin = 10.5, xmax = 11.5, ymin = 0.5, ymax = 12.5), fill = 'transparent', color = 'black', size = 1) +
  scale_x_continuous(breaks = seq(1, 21, 2), labels = seq(-10, 10, 2)) +
  scale_fill_manual(values = pal_bases) +
  labs(x = "Context position", y = "Genomic location", fill = "Nucleotide") +
  theme_minimal() +
  coord_cartesian(xlim = c(0.9, 21.1)) +
  theme(panel.grid = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(0, -10, 0, 0)) +
  guides(fill = guide_legend(position = 'top'))
a

b = ggplot(plot_data %>% filter(count > 3),
           aes(y = pos)) +
  geom_text(aes(label = BPS %>% sub('>','→', .), x = 'a',
                color = BPS),
            family = 'Aptos') +
  geom_richtext(aes(label = hotspot, x = 'b'),
                label.size = 0) +
  geom_text(aes(label = count, x = 'c')) +
  scale_x_discrete(breaks = c('a', 'b', 'c'), labels = c('BPS', 'Hotspot', 'Count'),
                   position = 'top') +
  scale_color_manual(values = pal_spec, guide = NULL) +
  coord_cartesian(xlim = c(1.2, 2.5)) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10, face = 'bold'),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(0, 10, 0, 20))

b

plot_grid(a, b, align = 'h', rel_widths = c(3.7, 1))


ggsave('figS7.svg', width = 9, height = 5)

#### Figure 3 / Figure S8 (GC>CG and AT>CG hotspots) ####

### run all this experiment with either chosen.BPS = 'GC>CG' or 'AT>CG ###
chosen.BPS = 'GC>CG'

c.runs.genome = genome.motifs %>%
  mutate(run.length = attr(regexpr("^C*", substring(motif, 12, 21), perl = TRUE), "match.length"),
         five.base = substr(motif, 10, 10)) %>%
  group_by(focus, five.base, run.length) %>%
  summarise(gen.count = n(), .groups = 'drop') %>%
  mutate(five.base = factor(five.base, levels = c('A', 'T', 'G', 'C'))) %>%
  filter(focus == substr(chosen.BPS, 1, 1))


# fills dataframe with missing combinations, gives them NA value
c.runs.genome %<>% complete(focus, five.base, run.length)

trunc_signif <- function(x) {
  pow <- floor(log10(abs(x)))
  digits = ifelse(pow > 2, 2, 1)
  scale <- 10^(pow - digits + 1)
  floor(x / scale) * scale
}

c = ggplot(c.runs.genome, aes(x = run.length, y = gen.count)) +
  geom_bar(aes(fill = five.base), stat = "identity", position = "dodge") +
  scale_x_continuous(breaks = function(x) (x[1]+x[2])/2) +
  scale_y_continuous(breaks = function(y) trunc_signif(y),
                     labels = function(y) prettyNum(y, big.mark = ','),
                     expand = expansion(mult = 0)) +
  scale_fill_manual(values = pal_bases, guide = NULL) +
  facet_wrap(~-run.length, scales = "free", ncol = 1) +
  labs(x = 'No. 3\' Cs', y = 'Count of Gs in genome') +
  theme_bw() +
  theme(strip.text = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y.left = element_text(family = 'Aptos Mono'),
        axis.text.x = element_text(size = 7))


plot_data = my_data %>%
  filter(BPS == chosen.BPS) %>%
  mutate(run.length = attr(regexpr("^C*", substring(motif, 12, 21), perl = TRUE), "match.length"),
         five.base = substr(motif, 10, 10)) %>%
  group_by(group, BPS, five.base, run.length) %>%
  summarise(count = n(), .groups = 'drop')


plot_data %<>%
  left_join(rates_summary %>% select(group, generations)) %>%
  mutate(focus = substr(BPS, 1, 1)) %>%
  left_join(c.runs.genome) %>%
  mutate(rate = count / generations / gen.count,
         group = factor(group, levels = groups),
         five.base = factor(five.base, levels = c('A', 'T', 'G', 'C')))


zero_data = plot_data %>%
  filter(run.length == 0)

plot_data %<>%
  filter(run.length > 0) %>%
  left_join(zero_data %>%
              select(group, BPS, five.base, rate) %>%
              rename(zero.rate = rate)) %>%
  mutate(rel.rate = rate / zero.rate)


# calculates n= labels
labels = plot_data %>%
  group_by(group, BPS, run.length) %>%
  summarise(count = sum(count), height = max(rel.rate) * 8, .groups = 'drop')


library(ggtext)
b = ggplot(plot_data, aes(x = run.length, y = rel.rate)) +
  geom_bar(aes(fill = five.base), stat = "identity",
           position = position_dodge2(preserve = "single")) +
  geom_richtext(data = labels, aes(x = run.length, y = height,
                                   label = count),
                inherit.aes = F, size = 3, label.size = 0, fill = NA,
                family = 'Aptos Mono') +
  facet_wrap(~group, ncol = 1) +
  scale_fill_manual(values = pal_bases) +
  scale_x_continuous(breaks = c(1:10)) +
  scale_y_log10(breaks = c(1, 1e2, 1e4),
                labels = c(1, paste0(rep(10, 2), '<sup>', c(2,4), '</sup>')),
                expand = expansion(mult = 0)) +
  coord_cartesian(xlim = c(0.8, max(plot_data$run.length))) +
  geom_hline(yintercept = 1) +
  labs(x = paste0('Number of consecutive Cs 3\' of focal ',
                  substr(chosen.BPS, 1, 1), ' nucleotide'),
       y = 'Fold increase in mutation rate') +
  theme_bw() +
  guides(fill = guide_legend(position = 'top', title = '5\' nucleotide')) +
  theme(axis.text.y = element_markdown(),
        # axis.title.x.bottom = element_text(hjust = 0),
        legend.text = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.ticks.x = element_blank(),
        # strip.text = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        legend.justification = 'center')


zero_data %<>%
  mutate(group = gsub('proofreading', '', group)) %>%
  mutate(group = gsub('MMR', '', group)) %>%
  mutate(group = factor(group, levels = c('(+) (+)', '(+) (-)', '(-) (+)', '(-) (-)')))


limits <- if (chosen.BPS == "GC>CG") {
  c(7e-11, 3.9e-10, 7e-9, 7e-9)
} else {
  c(8e-11, 3e-10, 4e-8, 12e-9)
}

labels <- zero_data %>%
  group_by(group, BPS, run.length) %>%
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
                                   label = paste0('n=', count)),
                size = 3, label.size = 0, label.padding = unit(0, 'pt'),
                vjust = 0.9, family = 'Aptos Mono') +
  facet_wrap(~group,
             scales = 'free_y',
             ncol = 1) +
  coord_cartesian(clip = "on", expand = FALSE) +
  # facet_wrap(~BPS + group, ncol = 4, scales = 'free_y') +
  scale_fill_manual(values = pal_bases) +
  scale_x_continuous(breaks = c(0:6), labels = c(paste(0:5), '6+')) +
  scale_y_continuous(labels = function(x) sci_label(x, digits = 0, drop_one = F)) +
  labs(x = '', y = 'A:T→C:G mutation rate') +
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
        axis.title.y = element_text(family = 'Aptos'),
        axis.text.y.left = element_markdown())


plot_grid(a, b, c, ncol = 3, rel_widths = c(0.18, 0.6, 0.18),
          align = 'h', axis = 'bt')



# main text for GC>CG
ggsave('fig3.png', height = 6, width = 7)

# supplementary for AT>CG
ggsave('figS8.png', height = 6, width = 7)




#### Figure 4 (strand bias) ####

# Coordinates for OriC and Ter in E. coli K-12
OriC <- 3925696
Ter <- 1640202

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
  group_by(group, BPS, strand) %>%
  summarise(count = n()) %>%
  left_join(genome.strands) %>%
  left_join(rates_summary %>% select(group, generations)) %>%
  mutate(rate = count / genome.count / generations / 2) %>% # divides by 2 to make equivalent with non-split graph
  mutate(group = factor(group, levels = groups)) %>%
  mutate(strand = factor(strand, levels = c('LDST', 'LGST')))

stats = plot_data %>%
  select(group, BPS, strand, rate) %>%
  pivot_wider(names_from = strand, values_from = rate) %>%
  mutate(
    total = LDST + LGST,
    LDST_pct = round(LDST / total * 100),
    LGST_pct = 100 - LDST_pct,
    ratio = paste0(LDST_pct, ":", LGST_pct),
    ymax = max(LDST, LGST)
  ) %>%
  select(group, BPS, ratio, ymax)

stats %<>%
  left_join(stats %>%
              group_by(group) %>%
              summarise(max = max(ymax)))

limits = tibble(
  group = groups,
  lim = c(6e-11, 2e-8, 6e-7, 6e-7)
)

ggplot(plot_data) +
  geom_bar(aes(x = BPS, y = rate, fill = strand, color = BPS),
           stat = 'identity', position = 'dodge') +
  geom_label(data = stats, aes(x = BPS, y = ymax + (0.2 * max), label = ratio), size = 2,
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
      geom_bar(aes(x = BPS, y = rate, fill = strand, color = BPS),
               stat = 'identity', position = 'dodge') +
      geom_label(data = stats2, aes(x = BPS, label = ratio), y = stats2$ymax + 0.1 * max(breaks[[i]]), size = 2,
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


legend = ggplot(plot_data, aes(x = BPS, y = rate, fill = strand)) +
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
  mutate(context = list(contexts)) %>%
  unnest(context) %>%
  mutate(value = substr(motif, 11 + context, 11 + context))

plot_data$value %<>% factor(levels = c('A', 'T', 'G', 'C'))

plot_data %<>%
  group_by(group, BPS, strand, context, value) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(!is.na(count)) %>%
  complete(group, BPS, strand, context, value, fill = list(count = 0)) # Fill missing combinations


plot_data %<>%
  mutate(focus = substr(BPS, 1, 1)) %>% # adds col for focus base (A or G)
  left_join(genome_context %>% rename(gen.perc = perc)) %>% # adds info for the genome_context bias at each position
  mutate(weight.count = (0.25 / gen.perc) * count) # weights by genome structure


plot_data %<>% left_join(
  plot_data %>%
    group_by(group, BPS, strand, context) %>%
    summarise(weight.total = sum(weight.count), .groups = 'drop')
) %>%
  mutate(weight.perc = weight.count / weight.total)

plot_data %<>% mutate(weight.perc = ifelse(context == 0 & value == 'G', 1, weight.perc))


b = ggplot(filter(plot_data, BPS == 'GC>CG'), aes(x = context, y = weight.perc, fill = value)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = pal_bases, name = 'Base', guide = NULL) +
  scale_y_continuous(breaks = c(0.25, 0.5, 0.75), sec.axis = dup_axis(name = 'Weighted frequency'),
                     labels = c('25%', '50%', '75%')) +
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
  filter(context != 0) %>%
  select(group, BPS, context, value, strand, weight.perc) %>%
  pivot_wider(names_from = strand, values_from = weight.perc)


regression_data <- plot_data %>%
  group_by(BPS, group) %>%
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
  left_join(regression_data, by = c("BPS", "group"))

# plot_data %<>% filter(BPS == 'GC>CG')

# Plot with manual annotations
c = ggplot(plot_data %>% filter(BPS == 'GC>CG'),
           aes(x = LDST, y = LGST, fill = value)) +
  geom_hline(yintercept = 0.25) +
  geom_vline(xintercept = 0.25) +
  geom_abline(slope = 1, intercept = 0, color = 'darkgrey', alpha = 0.5) +
  geom_abline(aes(intercept = intercept, slope = slope), color = "black", linetype = "dashed", size = 0.75) +
  geom_point(color = 'black', shape = 21, stroke = 0.1, size = 2.5) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75), labels = c('0%', '25%', '50%', '75%')) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75), labels = c('0%', '25%', '50%', '75%')) +
  coord_cartesian(xlim = c(0, max(c(plot_data$LDST, plot_data$LGST))), ylim = c(0, max(c(plot_data$LDST, plot_data$LGST)))) +
  scale_fill_manual(values = pal_bases, guide = NULL) +
  facet_grid(rows = vars(BPS), cols = vars(group),
             labeller = as_labeller(function(x) sub(">", "→", x))) +
  labs(x = 'Weighted frequency across G:C→C:G sites where G templates the... leading strand',
       y = '... lagging strand') +
  theme_minimal() +
  theme(
    strip.text.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(family = 'Aptos'),
    axis.title.x = element_text(family = 'Aptos'))


### Fig 4D

d = ggplot(regression_data) +
  geom_bar(aes(x = BPS, y = slope, fill = BPS), stat = 'identity') +
  geom_hline(yintercept = 1, size = 0.5, color = 'black') +
  geom_text(aes(x = BPS, y = slope - 0.12,
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
ggsave('fig4.png', width = 8*scale, height = 10*scale)



#### Figure S9 (Strand-comparison dot plots for all BPSs) ####

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
  facet_grid(rows = vars(BPS), cols = vars(group),
             labeller = as_labeller(function(x) sub(">", "→", x))) +
  labs(x = 'Nucleotide frequency across sites where purine templates the leading strand', y = 'Nucleotide frequency across sites where purine templates the lagging strand') +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    strip.text.y = element_text(family = 'Aptos'),
  ) +
  guides(fill = guide_legend(position = 'left', title = 'Nucleotide'))

ggsave('figS9.png', width = 7*1.1, height = 8.6*1.1)


### Supplementary Table 2

regression_data %>%
  select(group, BPS, slope, slope_test_p, r_squared) %>%
  arrange(group, BPS) %>%
  mutate(BPS = sub('>', '→', BPS)) %>%
  write_csv('Table S2.csv')




#### Figure S10 (Strand bias of A:T→C:G and G:C→C:G hotspots) ####

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
         focus = substr(BPS, 1, 1)) %>%
  group_by(group, BPS, focus, strand, run.length) %>%
  summarise(count = n(), .groups = 'drop')


plot_data %<>%
  left_join(rates_summary %>% select(group, generations)) %>%
  left_join(c.runs.genome) %>%
  mutate(rate = count / generations / gen.count,
         group = factor(group, levels = groups))


zero_data = plot_data %>%
  filter(run.length == 0)

plot_data %<>%
  left_join(zero_data %>% rename(zero.rate = rate) %>% select(group, BPS, strand, zero.rate)) %>%
  mutate(rel.rate = rate / zero.rate)

zero_data %<>%
  mutate(group = gsub('proofreading', '', group)) %>%
  mutate(group = gsub('MMR', '', group)) %>%
  mutate(group = factor(group, levels = c('(+) (+)', '(+) (-)', '(-) (+)', '(-) (-)')))



plot_list = list() ; for (chosen.BPS in c('GC>CG', 'AT>CG')) {
  
  plot_data2 = plot_data %>% filter(BPS == chosen.BPS)
  zero_data2 = zero_data %>% filter(BPS == chosen.BPS)
  
  labels <- zero_data2 %>%
    group_by(group, BPS, run.length) %>%
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
    # facet_wrap(~BPS + group, ncol = 4, scales = 'free_y') +
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
    group_by(group, BPS, run.length) %>%
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
                    substr(chosen.BPS, 1, 1)),
         y = 'Fold increase in mutation rate') +
    theme_bw() +
    guides(fill = guide_legend(position = 'top',
                               title = paste(sub('>', '→', chosen.BPS),
                                             'sites where',
                                             substr(chosen.BPS, 1, 1),
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

ggsave('figS10.png', width = 9, height = 8)




#### Figure 5 / Figure S12 (wider GC content) ####

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
        group_by(group, BPS, context) %>%
        summarise(GC_mean = mean(GC), GC_sd = sd(GC), n = n(), .groups = "drop")
    )
  if (i %% 10 == 0) {
    print(i)
  }
}

# the above takes ages to run, so save it so you don't have to run it again
window1000 = plot_data
write_csv(window1000, 'sliding_window.csv')

plot_data = read_csv('sliding_window.csv')

plot_data %<>%
  mutate(ts.tv = case_when(
    BPS %in% c('AT>GC', 'GC>AT') ~ 1,
    BPS %in% c('AT>CG', 'GC>TA') ~ 2,
    BPS %in% c('AT>TA', 'GC>CG') ~ 2
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
                aes(x = context, y = GC_mean, group = BPS, color = BPS)) +
      geom_hline(yintercept = length(which(mg1655 %in% c('G', 'C'))) / length(mg1655),
                 color = 'black', linetype = 'dashed') +
      geom_smooth(method = 'loess', span = 1, se = F) +
      geom_ribbon(aes(ymin = GC_mean - GC_sd / sqrt(n), ymax = GC_mean + GC_sd / sqrt(n), fill = BPS),
                  linetype = 0, alpha = 0.3) +
      scale_fill_manual(values = pal_spec,
                        labels = c('A:T→G:C', 'G:C→A:T', 'A:T→C:G', 'G:C→T:A', 'A:T→T:A', 'G:C→C:G'),
                        guide = NULL) +
      scale_color_manual(values = pal_spec, guide = NULL) +
      scale_y_continuous(breaks = c(0.45, 0.5, 0.55), labels = c('45%', '50%', '55%')) +
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
                aes(x = context, y = GC_mean, group = BPS, color = BPS)) +
      geom_hline(yintercept = length(which(mg1655 %in% c('G', 'C'))) / length(mg1655),
                 color = 'black', linetype = 'dashed') +
      geom_smooth(method = 'loess', span = 1, se = F) +
      geom_ribbon(aes(ymin = GC_mean - GC_sd / sqrt(n), ymax = GC_mean + GC_sd / sqrt(n), fill = BPS),
                  linetype = 0, alpha = 0.3) +
      scale_fill_manual(values = pal_spec,
                        labels = c('A:T→G:C', 'G:C→A:T', 'A:T→C:G', 'G:C→T:A', 'A:T→T:A', 'G:C→C:G'),
                        guide = NULL) +
      scale_color_manual(values = pal_spec, guide = NULL) +
      scale_y_continuous(position = 'right', breaks = c(0.45, 0.5, 0.55), labels = c('45%', '50%', '55%')) +
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
  draw_text('Starting position of 20 bp sliding context window', x = 0.5, y = 0.025, size = 11)


a = ggdraw() +
  draw_plot(plot, x = 0.025, y = 0, width = 0.975, height = 0.975) + 
  draw_label('Mean GC% of context window', x = 0.015, y = 0.5, angle = 90, size = 11)



### Figure 5B

window = 100

wrap.genome.wide = as.character(mg1655) %>% paste0(collapse = '')

wrap.genome.wide = paste0(
  substr(wrap.genome.wide, gen.length - window + 1, gen.length),
  wrap.genome.wide,
  substr(wrap.genome.wide, 1, window)
)

genome.GC = tibble(
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


ggplot(genome.GC) +
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
  genome.GC %>% mutate(BPS = ifelse(focus == 'A', 'AT>GC', 'GC>AT')),
  genome.GC %>% mutate(BPS = ifelse(focus == 'A', 'AT>CG', 'GC>TA')),
  genome.GC %>% mutate(BPS = ifelse(focus == 'A', 'AT>TA', 'GC>CG'))
) %>%
  mutate(BPS = factor(BPS, levels = spectrum))

plot_data %<>% pivot_longer(cols = c(five.GC, three.GC), names_to = 'context', values_to = 'GC') 
base_plot %<>% pivot_longer(cols = c(five.GC, three.GC), names_to = 'context', values_to = 'GC') 


stats = tibble()
for (i in 1:4) {
  for (j in 1:6) {
    for (k in 1:2) {
      data1 = filter(plot_data, group == groups[i] & BPS == spectrum[j] & context == c('five.GC', 'three.GC')[k])
      data2 = filter(base_plot, BPS == spectrum[j] & context == c('five.GC', 'three.GC')[k])
      
      ks.result = ks.test(data1$GC, data2$GC)
      
      stats %<>% bind_rows(tibble(
        group = groups[i], BPS = spectrum[j], context = c('five.GC', 'three.GC')[k],
        p = ks.result$p.value, D = ks.result$statistic
      ))
    }
  } ; print(paste(i, '/', 4))
}

stats %<>% mutate(p.adj = p.adjust(p, method = 'fdr'))

stats %<>% left_join(my_data %>%
                       group_by(group, BPS) %>%
                       summarise(count = paste0('n =  ', n())) %>%
                       mutate(context = 'five.GC'))



### Figure S12 (all BPS / snp combos)

# run only for supplementary figure
stats %<>% mutate(label = case_when(
  p.adj < 1e-10 ~ paste0('*p* < 10<sup>-10</sup><br>D = ', round(D, 2)),
  p.adj >= 0.05 ~ paste0('N.S.<br>D = ', round(D, 2)),
  TRUE ~ paste0('*p* = ', sci_label(p.adj, digits = 1), '<br>D = ', round(D, 2))
))


plot_data %<>% mutate(BPS = factor(BPS, levels = spectrum))
base_plot %<>% mutate(BPS = factor(BPS, levels = spectrum))

binwidth = 0.03
ggplot() +
  geom_histogram(data = base_plot %>% filter(context == 'five.GC'),
                 aes(y = GC, x = -after_stat(density)),
                 binwidth = binwidth, fill = 'white', color = 'black', linewidth = 0.5) +
  geom_histogram(data = base_plot %>% filter(context == 'three.GC'),
                 aes(y = GC, x = after_stat(density)),
                 binwidth = binwidth, fill = 'white', color = 'black', linewidth = 0.5) +
  geom_histogram(data = plot_data %>% filter(context == 'five.GC'),
                 aes(y = GC, x = -after_stat(density), fill = BPS),
                 binwidth = binwidth, alpha = 0.7) +
  geom_histogram(data = plot_data %>% filter(context == 'three.GC'),
                 aes(y = GC, x = after_stat(density), fill = BPS),
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
  facet_grid(cols = vars(group), rows = vars(BPS),
             labeller = as_labeller(function(x) sub('>', '→', x))) +
  labs(x = 'Probability Density (left is 5\' context, right is 3\' context)',
       y = 'GC-content of ±100 bp context',
       fill = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(family = 'Aptos'),
        legend.text = element_text(family = 'Aptos')) +
  guides(fill = guide_legend(position = 'top', nrow = 1))


ggsave('figS12.png', width = 7*1.1, height = 9*1.1)


# run only for main text figure
stats %<>%
  mutate(label = case_when(
    p.adj >= 0.05 ~ paste0('N.S.<br>', 'D = ', round(D, 2)),
    p.adj >= 0.01 ~ paste0('<i>p</i> < 0.05<br>D = ', round(D, 2)),
    p.adj >= 0.001 ~ paste0('<i>p</i> < 0.01<br>D = ', round(D, 2)),
    TRUE ~ paste0('<i>p</i> < 0.001<br>D = ', round(D, 2))
  ))


plot_data %<>% filter(BPS %in% spectrum[c(1, 4, 6)] & group %in% groups[c(1:2)]) %>%
  droplevels()
base_plot %<>% filter(BPS %in% spectrum[c(1, 4, 6)]) %>% 
  droplevels()
stats %<>% filter(BPS %in% spectrum[c(1, 4, 6)] & group %in% groups[c(1:2)]) %>%
  droplevels()


binwidth = 0.03
b = ggplot() +
  geom_histogram(data = base_plot %>% filter(context == 'five.GC'),
                 aes(y = GC, x = -after_stat(density)),
                 binwidth = binwidth, fill = 'white', color = 'black', linewidth = 0.5) +
  geom_histogram(data = base_plot %>% filter(context == 'three.GC'),
                 aes(y = GC, x = after_stat(density)),
                 binwidth = binwidth, fill = 'white', color = 'black', linewidth = 0.5) +
  geom_histogram(data = plot_data %>% filter(context == 'five.GC'),
                 aes(y = GC, x = -after_stat(density), fill = BPS),
                 binwidth = binwidth, alpha = 0.7) +
  geom_histogram(data = plot_data %>% filter(context == 'three.GC'),
                 aes(y = GC, x = after_stat(density), fill = BPS),
                 binwidth = binwidth, alpha = 0.7) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  geom_hline(yintercept = length(which(mg1655 %in% c('G', 'C'))) / length(mg1655),
             linetype = 'dashed', linewidth = 0.4) +
  geom_richtext(data = stats,
                aes(x = ifelse(context == 'five.GC', -4.5, 4.5), label = label),
                y = 0.16, label.padding = unit(0, 'in'), label.size = 0, hjust = 0.5, size = 2.6) +
  geom_richtext(data = stats,
                aes(x = -0.3, y = 0.82, label = count),
                label.size = 0, hjust = 0.5, size = 3, label.padding = unit(3, 'points')) +
  scale_y_continuous(breaks = c(0.25, 0.5, 0.75), labels = c('25%', '50%', '75%'), sec.axis = dup_axis()) +
  scale_x_continuous(breaks = c(-8, -4, 0, 4, 8), labels = c(8, 4, 0, 4, 8)) +
  scale_fill_manual(values = pal_spec, guide = NULL) +
  facet_wrap(~group+BPS, ncol = 6, drop = T) +
  # facet_grid(cols = vars(BPS), rows = vars(group)) +
  coord_cartesian(xlim = c(-6.5, 6.5)) +
  labs(x = 'Probability Density (left is 5\' context, right is 3\' context)',
       y = 'GC% of ±100 bp',
       fill = NULL) +
  theme_minimal() +
  theme(strip.text.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.right = element_blank())


plot_grid(a, b, ncol = 1, rel_heights = c(0.8, 0.4))

ggsave('fig5.png', width = 9.5, height = 5)




#### Figure S11 (Bias towards GC-rich regions for G:C→C:G is not exclusive to GC3+ mutations) ####

plot_data = my_data %>%
  filter(BPS == 'GC>CG') %>%
  mutate(run.length = attr(regexpr("^C*", substring(motif, 12, 21), perl = TRUE), "match.length"),
         hotspot = ifelse(run.length > 2, TRUE, FALSE)) %>%
  mutate(context = list(contexts)) %>%
  unnest(context) %>%
  rowwise() %>%
  mutate(GC = GC_context(mg1655, pos, context, ref)) %>%
  ungroup() %>%
  group_by(group, BPS, hotspot, context) %>%
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


ggsave('figS11.png', width = 10, height = 5)

