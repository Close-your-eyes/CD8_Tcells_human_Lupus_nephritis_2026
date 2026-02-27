
# download files from figshare: 10.6084/m9.figshare.31417085

# install required package
# devtools::install_github("close-your-eyes/scexpr")
# devtools::install_github("close-your-eyes/colrr")
# install.packages("ggplot2")
# install.packages("ggtext")
# install.packages("ggpubr")
# install.packages("dplyr")
# install.packages("patchwork")
# install.packages("cowplot")
# install.packages("abdiv")
# install.packages("purrr")
# install.packages("rlang")


library(scexpr)
library(ggplot2) 

# read seurat objects in getwd()
SO <- readRDS("SO_urine_blood_hg38_CR7_non_contam_rep1_wo_Introns_strict_QC_patient_integration_SCT_harmony_1_400_15_230515-113449.rds")
SO_CD8_eff <- readRDS("SO_urine_blood_hg38_CR7_non_contam_rep1_wo_Introns_strict_QC_patient_integration_CD8_eff_SCT_harmony_1_300_12_231016-163715.rds")


## ---- Fig1A -------
Fig1A <- scexpr::feature_plot2(SO,
                               features = "clusters",
                               pt_size = 0.1,
                               reduction = "tsne",
                               name_anno_pos = NULL,
                               col_pal_d_args = list(name = SO@misc$clusters_color),
                               axes_lim_set = list(x = c(-60,60))) +
  ggtext::geom_richtext(data = SO@misc$label_df, aes(label = label), label.colour = NA, size = 3,
                        color = "black", fill = "grey95",
                        label.padding = unit(c(rep(0.05,4)), "lines"),
                        alpha = 0.8) +
  theme(strip.text.x = ggtext::element_markdown(size = 10),
        strip.background = element_rect(color = "white", fill = "grey95"),
        legend.position = "none") +
  facet_wrap(vars("CD4<sup>+</sup> or CD8<sup>+</sup> T cells")) +
  annotate("text", x = -40, y = 51, label = paste0("n = ", length(Seurat::Cells(SO))), size = 3)


## ---- Fig1B ------
Fig1B <- scexpr::feature_plot2(SO,
                               features = "body_fluid",
                               pt_size = 0.1,
                               reduction = "tsne",
                               split_feature = "body_fluid",
                               col_pal_d_args = list(name = colrr::col_pal("custom")[c(5,2)]),
                               contour_feature = "clusters",
                               name_anno_pos = NULL,
                               contour_args = list(breaks = 0.15, linewidth = 0.2, color = "black"),
                               contour_same_across_split = T,
                               axes_lim_set = list(c(-52,52), c(-55, 52))) +
  theme(strip.text.x = ggtext::element_markdown(size = 10),
        strip.background = element_rect(color = "white", fill = "grey95"),
        legend.position = "none")


## ---- Fig1C -------
SO@meta.data$body_fluid2 <- ifelse(SO@meta.data$body_fluid == "blood",
                                   paste0("blood*\nn = ", table(SO@meta.data$body_fluid)["blood"]),
                                   paste0("urine\nn = ", table(SO@meta.data$body_fluid)["urine"]))
SO@meta.data$body_fluid2 <- as.factor(SO@meta.data$body_fluid2)
Fig1C <-
  scexpr::composition_barplot(SO,
                              y = "abs",
                              x_cat = "clusters",
                              fill_cat = "body_fluid2",
                              flip = T,
                              label_only_largest = T,
                              label_rel_decimals = 0,
                              label_size_abs = 4,
                              label_rel_nudge = list(CD8_CCR7 = c(0, 0),
                                                     CD4_IL7R = c(0, 600),
                                                     CD8_CCL5 = c(0, 800),
                                                     CD4_FOXP3 = c(0, 800),
                                                     CD8_MKI67 = c(0, 1000),
                                                     CD8_HSPA1A = c(0, 1300)),
                              plot_rel_labels = T,
                              col_pal = c(colrr::col_pal("custom")[c(5,2)]))[["plot"]] +
  colrr::theme_material(white = T) +
  theme(axis.title = element_blank(), legend.position.inside = c(0.6,0.90),
        axis.text.y = ggtext::element_markdown(color = SO@misc$clusters_color, face = "bold", size = 10)) +
  guides(fill = guide_legend(nrow = 1, position = "inside")) +
  scale_y_continuous(breaks = c(0,4000,8000), labels = c("n = 0", 4000, 8000)) +
  annotate("text", y = 2900, x = 2.5, label = "*CD4:CD8 ratio adjusted", hjust = 0)


## ---- Fig1D ------
heatmap_features <- c("CD3E", "CD8A", "CD8B", "CD4",
                      "GZMK", "GZMB", "PRF1", "NKG7", "CCL5", "HLA-DRB1", "CD38",
                      "MKI67", "TOP2A", "CENPF",
                      "HSPA1B", "DNAJB1", "HSPA1A",
                      "MAL", "LEF1", "CCR7", "SELL", "RPS13", "RPS3A", "RPS12",
                      "AQP3", "IL7R", "LTB", "S100A4",
                      "TNFRSF18", "TNFRSF4", "CTLA4", "FOXP3", "IKZF2")
Fig1D <- scexpr::heatmap_pseudobulk2(SO,
                                     color = "white",
                                     features = heatmap_features,
                                     meta_col = "clusters",
                                     feature_order = "none",
                                     group_order = "none",
                                     theme = colrr::theme_material(white = T),
                                     levels_plot = levels(SO@meta.data$clusters),
                                     theme_args = list(axis.title = ggplot2::element_blank(),
                                                       axis.text.x = ggplot2::element_text(),
                                                       axis.text.y = ggplot2::element_text(size = 10),
                                                       legend.position = "bottom",
                                                       legend.direction = "horizontal"),
                                     legend_fill_args = list(label.theme = ggplot2::element_text(size = 10),
                                                             title.theme = ggplot2::element_text(size = 10),
                                                             title.position = "right",
                                                             title = "transcription\nlevel [z-score]",
                                                             title.hjust = 0.5,
                                                             barwidth = 0.8,
                                                             barheight = 14,
                                                             order = 1),
                                     legend_size_args = list(label.theme = ggplot2::element_text(size = 10),
                                                             title.theme = ggplot2::element_text(size = 10),
                                                             title.position = "right",
                                                             title = "transcription\nfrequency [%]",
                                                             title.hjust = 0.5,
                                                             label.position = "bottom",
                                                             order = 2,
                                                             ncol = NULL,
                                                             nrow = NULL,
                                                             override.aes = list(color = "black")),
                                     colorsteps = c(-1,0,1),
                                     dotplot = T,
                                     axes_flip = T)[["plot"]] +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.text.y = ggtext::element_markdown(color = SO@misc$clusters_color, face = "bold"),
    legend.box.margin = margin(-0.4,0,0,0, unit = "cm"),
    legend.position = "bottom",
    axis.title.x = element_blank(), axis.title.y = element_blank())


## ---- Fig1E ------
names(SO_CD8_eff@misc$label_df) <- tolower(names(SO_CD8_eff@misc$label_df))
Fig1E <- scexpr::feature_plot2(SO_CD8_eff,
                               features = "clusters",
                               reduction = "tsne",
                               pt_size = 0.1,
                               col_pal_d_args = list(name = SO_CD8_eff@misc$clusters_color),
                               name_anno_pos = NULL,
                               axes_lim_set = list(x = c(-70,60))) +
  ggtext::geom_richtext(data = SO_CD8_eff@misc$label_df, aes(label = label), label.colour = NA, size = 3,
                        color = "black", fill = "grey95",
                        label.padding = unit(c(rep(0.05,4)), "lines"),
                        alpha = 0.8) +
  theme(strip.text.x = ggtext::element_markdown(size = 10),
        strip.background = element_rect(color = "white", fill = "grey95"),
        legend.position = "none",
        axis.title.y = element_blank()) +
  facet_wrap(vars("CD8_CCL5 transcriptomes")) +
  annotate("text", x = -47, y = 53, label = paste0("n = ", length(Seurat::Cells(SO_CD8_eff))), size = 3)


## ---- Fig1F --------
Fig1F <- scexpr::feature_plot2(SO_CD8_eff,
                               features = "body_fluid",
                               pt_size = 0.1,
                               reduction = "tsne",
                               split_feature = "body_fluid",
                               col_pal_d_args = list(name = colrr::col_pal("custom")[c(5,2)]),
                               contour_feature = "clusters",
                               name_anno_pos = NULL,
                               contour_args = list(breaks = 0.3, linewidth = 0.2, color = "black"),
                               contour_same_across_split = T,
                               axes_lim_set = list(x = c(-60,63), y = c(-66, 50))) +
  theme(strip.text.x = ggtext::element_markdown(size = 10),
        strip.background = element_rect(color = "white", fill = "grey95"),
        legend.position = "none")

## ----- Fig1G -------
SO_CD8_eff@meta.data$body_fluid2 <- ifelse(SO_CD8_eff@meta.data$body_fluid == "blood",
                                           paste0("blood\nn = ", table(SO_CD8_eff@meta.data$body_fluid)["blood"]),
                                           paste0("urine\nn = ", table(SO_CD8_eff@meta.data$body_fluid)["urine"]))
SO_CD8_eff@meta.data$body_fluid2 <- as.factor(SO_CD8_eff@meta.data$body_fluid2)
Fig1G <-
  scexpr::composition_barplot(SO_CD8_eff,
                              y = "abs",
                              x_cat = "clusters",
                              fill_cat = "body_fluid2",
                              flip = T,
                              label_only_largest = T,
                              label_rel_decimals = 0,
                              label_rel_nudge = list(c(0,0),
                                                     c(0,0),
                                                     c(0,0),
                                                     c(0,0),
                                                     c(0,0)),
                              plot_rel_labels = T,
                              col_pal = c(colrr::col_pal("custom")[c(5,2)]))[["plot"]] +
  colrr::theme_material(white = T) + 
  theme(axis.title = element_blank(), legend.position.inside = c(0.6,0.89),# legend.justification = c(1.2,1.1),
        axis.text.y = ggtext::element_markdown(color = SO_CD8_eff@misc$clusters_color, face = "bold", size = 10),
        legend.key.size = unit(1, "lines")) +
  guides(fill = guide_legend(nrow = 1, position = "inside")) +
  scale_y_continuous(breaks = c(0,750,1500), labels = c("n = 0", 750, 1500)) #name = "n transcriptome"



## ---- Fig1H --------
heatmap_features_CD8eff <- c("GZMK", "CD74", "CXCR3", "CXCR4", "RGS1", "HLA-DRB1", "CD28", "CD38",
                             "CD27", "CXCR6", "GZMB", "GNLY", "PRF1", "LGALS1",  "HAVCR2", "KLRC1",
                             "IL7R", "RPS3A", "TXNIP",
                             "GZMH", "CX3CR1", "FGFBP2", "KLRG1", "KLRD1", "S100A4", "NKG7",
                             "IFNG", "TNF", "CCL3", "CCL4", "CCL3L1", "CCL4L2", "NR4A1", "NR4A2", "NR4A3")
Fig1H <- scexpr::heatmap_pseudobulk2(SO_CD8_eff,
                                     meta_col = "clusters",
                                     features = heatmap_features_CD8eff,
                                     theme = colrr::theme_material(white = T),
                                     theme_args = list(axis.title = ggplot2::element_blank(),
                                                       axis.text.x = ggplot2::element_text(),
                                                       axis.text.y = ggplot2::element_text(size = 10, face = "italic"),
                                                       legend.position = "bottom",
                                                       legend.direction = "horizontal"),
                                     legend_fill_args = list(label.theme = ggplot2::element_text(size = 10),
                                                             title.theme = ggplot2::element_text(size = 10),
                                                             title.position = "right",
                                                             title = "transcription\nlevel [z-score]",
                                                             title.hjust = 0.5,
                                                             title.vjust = 1.3,
                                                             barwidth = 0.8,
                                                             barheight = 14,
                                                             order = 1),
                                     legend_size_args = list(label.theme = ggplot2::element_text(size = 10),
                                                             title.theme = ggplot2::element_text(size = 10),
                                                             title.position = "right",
                                                             title = "transcription\nfrequency [%]",
                                                             title.hjust = 0.5,
                                                             label.position = "bottom",
                                                             order = 2,
                                                             ncol = NULL,
                                                             nrow = NULL,
                                                             override.aes = list(color = "black")),
                                     levels_plot = levels(SO_CD8_eff@meta.data$clusters),
                                     dotplot = T,
                                     colorsteps = c(-1,0,1),
                                     axes_flip = T)[["plot"]] +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        axis.text.y = ggtext::element_markdown(color = SO_CD8_eff@misc$clusters_color, face = "bold"),
        legend.box.margin = margin(-0.4,0,0,0, unit = "cm"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())



## ---- Fig2A -----
SO@meta.data$expanded2 <- ifelse(SO@meta.data$expanded, "expanded", "unique")
names(SO@misc$label_df) <- tolower(names(SO@misc$label_df))
Fig2A <- scexpr::feature_plot2(SO,
                               features = "expanded2",
                               get_data_args = list(na_rm = T),
                               pt_size = 0.1,
                               reduction = "tsne",
                               col_pal_d_args = list(name = c("tomato2", "grey70")),
                               contour_feature = "clusters",
                               contour_args = list(color = "black", breaks = 0.2, linewidth = 0.2),
                               name_anno_pos = NULL,
                               axes_lim_set = list(x=c(-50,50),y=c(-54,54))) +
  ggtext::geom_richtext(data = SO@misc$label_df, aes(label = label), label.colour = NA, size = 3,
                        color = "black", fill = "grey95",
                        label.padding = unit(c(rep(0.05,4)), "lines"),
                        alpha = 0.8) +
  facet_wrap(vars("CD4<sup>+</sup> or CD8<sup>+</sup> T cells")) +
  theme(strip.text.x = ggtext::element_markdown(size = 10), strip.background = element_rect(color = "white", fill = "grey95"),
        legend.position = "inside", legend.position.inside = c(0.1,1), legend.justification = c(0,1))




## ---- Fig2C ------
group_name1 <- "CD8<br>clusters"
group_name2 <- "*CCR7*<sup> -</sup> CD8<br>clusters"
freq_stats <- purrr::map_dfr(setNames(c(group_name1, group_name2),
                                      c(group_name1, group_name2)), function(x) {
                                        
                                        cluster_select <- if (x == group_name2) {
                                          c("CD8_CCL5", "CD8_MKI67", "CD8_HSPA1A")
                                        } else {
                                          c("CD8_CCL5", "CD8_MKI67", "CD8_HSPA1A", "CD8_CCR7")
                                        }
                                        
                                        SO@meta.data %>%
                                          dplyr::filter(clusters %in% cluster_select) %>%
                                          dplyr::filter(!is.na(cl_name))  %>%
                                          dplyr::select(body_fluid, patient, cl_name, expanded) %>%
                                          dplyr::group_by(body_fluid, patient) %>%
                                          dplyr::summarise(freq_expanded = sum(expanded)/dplyr::n(), .groups = "drop")
                                      }, .id = "group") %>%
  dplyr::mutate(group = factor(group, levels = c(group_name1, group_name2)))

Fig2C <- ggplot(freq_stats, aes(x = body_fluid, y = freq_expanded*100)) +
  geom_line(aes(group = patient)) +
  geom_point(size = 3, aes(color = patient)) +
  colrr::theme_material(white = T) + 
  theme(strip.text.x = ggtext::element_markdown(size = 12), legend.text = element_text(size = 12)) +
  ggpubr::stat_compare_means(comparisons = list(c(1,2)), paired = T, method = "t.test", label = "p.signif") + 
  scale_color_manual(values = SO@misc$patient_color) +
  labs(x = NULL, y = "frequency of expanded clones [%]") +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  facet_wrap(vars(group), axes = "all", axis.labels = "margins")


## ----- Fig2D ------
temp_col_pal <- stats::setNames(c("grey40", "hotpink1", "green3", "grey85"), c("exp", "exp+ovlp", "ovlp", "unique"))
p_exp_ovlp_all <- scexpr::feature_plot2(SO,
                                        features = "exp_olvp",
                                        get_data_args = list(na_rm = T),
                                        pt_size = 0.35,
                                        split_feature = "body_fluid",
                                        col_pal_d_args = list(name = temp_col_pal),
                                        reduction = "tsne",
                                        name_anno_pos = NULL,
                                        contour_args = list(color = "black", breaks = 0.2, linewidth = 0.2),
                                        contour_feature = "clusters",
                                        facet_grid_row_var = "CD4<sup>+</sup> or CD8<sup>+</sup> T cells",
                                        axes_lim_set = list(x = c(-55,NA))) +
  ggtext::geom_richtext(data = SO@misc$label_df, aes(label = label), size = 2.5, label.padding = unit(rep(0.01,4),"lines"), label.color = NA) +
  theme(strip.text.y = ggtext::element_markdown(size = 10),
        strip.background = element_rect(fill = "grey95", color = "white"),
        strip.background.y = element_rect(fill = "grey95", color = "white"))

names(SO_CD8_eff@misc$label_df) <- tolower(names(SO_CD8_eff@misc$label_df))
p_exp_ovlp_CD8_eff <- scexpr::feature_plot2(SO_CD8_eff,
                                            features = "exp_olvp",
                                            get_data_args = list(na_rm = T),
                                            pt_size = 0.35,
                                            split_feature = "body_fluid",
                                            col_pal_d_args = list(name = temp_col_pal),
                                            reduction = "tsne",
                                            name_anno_pos = NULL,
                                            contour_args = list(color = "black", breaks = 0.3, linewidth = 0.2),
                                            contour_feature = "clusters",
                                            facet_grid_row_var = "CD8_CCL5 transcriptomes",
                                            strip.background.y = element_rect(fill = "grey95", color = "white"),
                                            axes_lim_set = list(x = c(-59,62), y = c(-67, 40))) +
  ggtext::geom_richtext(data = SO_CD8_eff@misc$label_df, aes(label = label), size = 2.5,
                        label.padding = unit(rep(0.01,4),"lines"), label.color = NA) +
  theme(strip.text.x = element_blank(),
        strip.background = element_rect(fill = "grey95", color = "white"))

p_exp_ovlp <- cowplot::plot_grid(patchwork::wrap_plots(p_exp_ovlp_all,
                                                       patchwork::plot_spacer(),
                                                       p_exp_ovlp_CD8_eff, guides = "collect",
                                                       ncol = 1, heights = c(0.5,-0.026,0.5)))
Fig2D <- p_exp_ovlp


## ------ Fig2E --------
## cluster comparison for exp+ovlp
cl_data <- scexpr::get_data(SO_CD8_eff, c("exp_olvp", "patient", "body_fluid", "clusters", "cl_name", "Epitope.species_observed",
                                          "expanded", "overlapping"), reduction = "tsne", try_df = T) %>%
  dplyr::filter(patient != "Pat1") %>%
  dplyr::mutate(patient = as.factor(as.character(patient))) %>%
  tidyr::drop_na(cl_name)

cl_data_summary <-
  cl_data %>%
  dplyr::group_by(exp_olvp, clusters, patient, body_fluid, .drop = F) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::left_join(cl_data %>%
                     dplyr::group_by(clusters, patient, body_fluid, .drop = F) %>%
                     dplyr::summarise(n_cluster = dplyr::n())) %>%
  dplyr::mutate(freq = n/n_cluster)

cl_data_summary_blood <-
  cl_data_summary %>%
  dplyr::filter(exp_olvp == "exp+ovlp") %>%
  dplyr::filter(body_fluid == "blood")
p_cluster_comp_expovlp_blood_abs <- ggplot(cl_data_summary_blood, aes(x = clusters, y = n)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 3, aes(color = patient)) +
  colrr::theme_material(white = T) +
  ggpubr::stat_compare_means(method = "t.test", paired = T, comparisons = list(c(1,4), c(2,4), c(3,4), c(5,4)), label = "..p.format..",
                             #symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "n.s.")),
                             size = 3.5) +
  scale_color_manual(values = SO_CD8_eff@misc$patient_color) +
  labs(y = "transcriptomes belonging<br>to exp + ovlp clonotypes", x = NULL) +
  theme(axis.text.x = ggtext::element_markdown(face = "bold", color = SO_CD8_eff@misc$clusters_color[which(names(SO_CD8_eff@misc$clusters_color) %in% cl_data_summary_blood$clusters)], size = 10, angle = 20, hjust = 1, vjust = 1),
        axis.title.y = ggtext::element_markdown(), strip.background = element_rect(fill = "grey95", color = "white"),
        strip.text.x = element_text(vjust = 0, size = 12), legend.text = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  facet_wrap(vars(body_fluid))

cl_data_summary_urine <-
  cl_data_summary %>%
  dplyr::filter(exp_olvp == "exp+ovlp") %>%
  dplyr::filter(body_fluid == "urine")
p_cluster_comp_expovlp_urine_abs <- ggplot(cl_data_summary_urine, aes(x = clusters, y = n)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 3, aes(color = patient)) +
  colrr::theme_material(white = T) +
  ggpubr::stat_compare_means(method = "t.test", paired = T, comparisons = list(c(1,2), c(1,3), c(1,4), c(1,5)), label = "..p.format..",
                             #symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "n.s.")),
                             size = 3.5) +
  scale_color_manual(values = SO_CD8_eff@misc$patient_color) +
  labs(y = NULL, x = NULL) +
  theme(axis.text.x = ggtext::element_markdown(face = "bold", color = SO_CD8_eff@misc$clusters_color, size = 10, angle = 20, hjust = 1, vjust = 1),
        strip.background = element_rect(fill = "grey95", color = "white"),
        strip.text.x = element_text(vjust = 0, size = 12), legend.text = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  facet_wrap(vars(body_fluid))

# do not put relative values into main figure
#p_cluster_comp_expovlp_blood, p_cluster_comp_expovlp_urine,
p_cluster_comp_expovlp <- cowplot::plot_grid(patchwork::wrap_plots(p_cluster_comp_expovlp_blood_abs,
                                                                   p_cluster_comp_expovlp_urine_abs,
                                                                   nrow = 1,
                                                                   guides = "collect"))
Fig2E <- p_cluster_comp_expovlp


## ------ Fig2F ---------
n_cluster_per_cl <-
  cl_data %>%
  dplyr::group_by(cl_name) %>%
  dplyr::summarise(n_cluster = dplyr::n_distinct(clusters))

cl_data <-
  cl_data %>%
  dplyr::left_join(n_cluster_per_cl)

cl_data_summary_shared_between_cluster2 <-
  cl_data %>%
  dplyr::filter(body_fluid == "urine") %>%
  dplyr::distinct(patient, cl_name, clusters, n_cluster) %>%  # here decides if clones or clonotypes are plotted; # n_cluster
  dplyr::group_by(patient) %>%
  dplyr::summarise(total = dplyr::n(), shared = sum(n_cluster > 1), .groups = "drop") %>%
  dplyr::mutate(rel_shared = shared/total*100)
cl_data_summary_shared_between_cluster2

cl_data_summary_shared_between_cluster3 <-
  cl_data %>%
  dplyr::filter(expanded) %>%
  dplyr::filter(body_fluid == "urine") %>%
  dplyr::distinct(patient, cl_name, clusters, n_cluster) %>% # here decides if clones or clonotypes are plotted
  dplyr::group_by(patient) %>%
  dplyr::summarise(total = dplyr::n(), shared = sum(n_cluster > 1), .groups = "drop") %>%
  dplyr::mutate(rel_shared = shared/total*100)
cl_data_summary_shared_between_cluster3

## jaccard overlap matrix plot ##
cl_data_sub <-
  cl_data %>%
  dplyr::filter(expanded) %>%
  dplyr::distinct(patient, cl_name, clusters, body_fluid) %>%
  dplyr::filter(body_fluid == "urine")
cl_data_sub_split <- split(cl_data_sub, cl_data_sub$clusters)

ovlp_index_df <- purrr::map_dfr(setNames(levels(cl_data_sub$patient), levels(cl_data_sub$patient)), function(z) {
  purrr::map_dfr(setNames(combn(names(cl_data_sub_split), 2, simplify = F), combn(names(cl_data_sub_split), 2, simplify = F)), function(x) {
    sub1 <- cl_data_sub_split[[x[1]]]
    sub1 <- dplyr::filter(sub1, patient == z)
    sub2 <- cl_data_sub_split[[x[2]]]
    sub2 <- dplyr::filter(sub2, patient == z)
    
    xx <- stack(table(sub1$cl_name))
    yy <- stack(table(sub2$cl_name))
    
    names(xx)[1] <- x[1]
    names(yy)[1] <- x[2]
    
    zz <- dplyr::full_join(xx, yy, by = "ind") %>%
      dplyr::mutate(!!rlang::sym(x[1]) := ifelse(is.na(!!rlang::sym(x[1])), 0, !!rlang::sym(x[1]))) %>%
      dplyr::mutate(!!rlang::sym(x[2]) := ifelse(is.na(!!rlang::sym(x[2])), 0, !!rlang::sym(x[2])))
    
    #y <- length(intersect(sub1$cl_name, sub2$cl_name))/mean(c(length(sub1$cl_name), length(sub2$cl_name)))
    data.frame(ovlp_index = (1 - abdiv::jaccard(zz[,x[1],drop=T], zz[,x[2],drop=T]))*100, c1 = x[1], c2 = x[2])
  })
}, .id = "patient")
# order of compared clusters is irrelevant
ovlp_index_df$c1_c2 <- unlist(sapply(purrr::map2(ovlp_index_df$c1, ovlp_index_df$c2, function(x,y) sort(c(x,y))), paste, collapse = " vs. ", simplify = F))

for (i in seq_along(ovlp_index_df$c1_c2)) {
  for (j in names(SO_CD8_eff@misc$clusters_color)) {
    ovlp_index_df$c1_c2[i] <- gsub(j, paste0("<span style = 'color:", SO_CD8_eff@misc$clusters_color[j], ";'>", j, "</span>"), ovlp_index_df$c1_c2[i])
  }
}

p_ovlp_boxplots <- ggplot(ovlp_index_df, aes(x = reorder(c1_c2, ovlp_index, FUN = median), y = ovlp_index)) +
  geom_boxplot() +
  geom_point(size = 3, aes(color = patient)) +
  colrr::theme_material(white = T) +
  theme(axis.text.x = ggtext::element_markdown(angle = 30, hjust = 1, vjust = 1, size = 10), legend.position = "left",
        legend.text = element_text(size = 12)) +
  labs(y = "Jaccard index * 100\nbetween urinary clusters", x = NULL) +
  scale_color_manual(values = SO_CD8_eff@misc$patient_color) +
  ggpubr::stat_compare_means(paired = T, ref.group = 9, method = "t.test", label = "p.signif",
                             symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "n.s.")),
                             size = 4, vjust = 1) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  ylim(NA, 38)

Fig2F <- p_ovlp_boxplots


## ---- Fig2G -----
## jaccard index of CD8_CX3CR1 blood cluster to urinary clusters
cl_data_sub2 <-
  cl_data %>%
  dplyr::filter(expanded) %>%
  dplyr::distinct(patient, cl_name, clusters, body_fluid) %>%
  dplyr::filter(body_fluid == "urine" | (body_fluid == "blood" & clusters == "CD8_CX3CR1")) %>%
  dplyr::mutate(body_fluid_cluster = paste0(clusters, "_", body_fluid))
cl_data_sub2_split <- split(cl_data_sub2, cl_data_sub2$body_fluid_cluster)

ovlp_index_df2 <- purrr::map_dfr(setNames(levels(cl_data_sub2$patient), levels(cl_data_sub2$patient)), function(z) {
  purrr::map_dfr(setNames(combn(names(cl_data_sub2_split), 2, simplify = F), combn(names(cl_data_sub2_split), 2, simplify = F)), function(x) {
    sub1 <- cl_data_sub2_split[[x[1]]]
    sub1 <- dplyr::filter(sub1, patient == z)
    sub2 <- cl_data_sub2_split[[x[2]]]
    sub2 <- dplyr::filter(sub2, patient == z)
    
    xx <- stack(table(sub1$cl_name))
    yy <- stack(table(sub2$cl_name))
    
    names(xx)[1] <- x[1]
    names(yy)[1] <- x[2]
    
    zz <- dplyr::full_join(xx, yy, by = "ind") %>%
      dplyr::mutate(!!rlang::sym(x[1]) := ifelse(is.na(!!rlang::sym(x[1])), 0, !!rlang::sym(x[1]))) %>%
      dplyr::mutate(!!rlang::sym(x[2]) := ifelse(is.na(!!rlang::sym(x[2])), 0, !!rlang::sym(x[2])))
    
    #y <- length(intersect(sub1$cl_name, sub2$cl_name))/mean(c(length(sub1$cl_name), length(sub2$cl_name)))
    data.frame(ovlp_index = (1 - abdiv::jaccard(zz[,x[1],drop=T], zz[,x[2],drop=T]))*100, c1 = x[1], c2 = x[2])
  })
}, .id = "patient") %>%
  dplyr::filter(c1 == "CD8_CX3CR1_blood")

ovlp_index_df2_order_summary <-
  ovlp_index_df2 %>%
  dplyr::group_by(c2) %>%
  dplyr::summarise(median_ovlp_index = median(ovlp_index)) %>%
  dplyr::arrange(median_ovlp_index)
ovlp_index_df2$c2 <- factor(ovlp_index_df2$c2, ovlp_index_df2_order_summary$c2)
cluster_color_temp <- SO_CD8_eff@misc$clusters_color
names(cluster_color_temp) <- paste0(names(cluster_color_temp), "_urine")

p_ovlp_boxplots2 <- ggplot(ovlp_index_df2, aes(x = c2, y = ovlp_index)) +
  geom_boxplot() +
  geom_point(size = 3, aes(color = patient)) +
  colrr::theme_material(white = T) +
  theme(axis.text.x = ggtext::element_markdown(angle = 30, hjust = 1, vjust = 1, size = 10, color = cluster_color_temp[levels(ovlp_index_df2$c2)]),
        legend.position = "none") +
  labs(y = "Jaccard index * 100\nwith CD8_CX3CR1_blood", x = NULL) +
  scale_color_manual(values = SO_CD8_eff@misc$patient_color) +
  ggpubr::stat_compare_means(paired = T, ref.group = 5, method = "t.test", label = "p.signif",
                             symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "n.s.")),
                             size = 4, vjust = 1) +
  ylim(NA, 35)
Fig2G <- p_ovlp_boxplots2


## ----- Fig3E ---------
epitope_colors <- setNames(c(
  "#104dff",
  "#ff4c4c",
  "#8a4b00",
  "#ceff00",
  "#ffbbff",
  "#00b235",
  "#fbef00",
  "grey85",
  "grey85",
  "grey50",
  "grey85"
), nm = c(
  "CMV",
  "EBV",
  "HSV-2",
  "HomoSapiens",
  "multi",
  "other",
  "CEFX",
  "none",
  "no reactivity found",
  "*in vitro* tested",
  "not tested"
))

SO@meta.data <-
  SO@meta.data %>%
  dplyr::mutate(in_vitro_tested = ifelse(is.na(in_vitro_tested), F, in_vitro_tested)) %>%
  dplyr::mutate(Epitope.species_observed = ifelse(is.na(Epitope.species_observed), "none", Epitope.species_observed)) %>%
  dplyr::mutate(in_vitro_tested_Epitope.species_observed = ifelse(Epitope.species_observed != "none", Epitope.species_observed, in_vitro_tested)) %>%
  dplyr::mutate(in_vitro_tested_Epitope.species_observed = ifelse(in_vitro_tested_Epitope.species_observed == "FALSE", "not tested", in_vitro_tested_Epitope.species_observed)) %>%
  dplyr::mutate(in_vitro_tested_Epitope.species_observed = ifelse(in_vitro_tested_Epitope.species_observed == "TRUE", "*in vitro* tested", in_vitro_tested_Epitope.species_observed)) %>%
  dplyr::mutate(in_vitro_tested_Epitope.species_observed = ifelse(is.na(cl_name), "no TCR", in_vitro_tested_Epitope.species_observed)) %>%
  dplyr::mutate(in_vitro_tested_Epitope.species_observed = ifelse(in_vitro_tested_Epitope.species_observed %in% c("not tested", "no TCR"), "not tested", in_vitro_tested_Epitope.species_observed))

SO_CD8_eff@meta.data <-
  SO_CD8_eff@meta.data %>%
  dplyr::mutate(in_vitro_tested = ifelse(is.na(in_vitro_tested), F, in_vitro_tested)) %>%
  dplyr::mutate(Epitope.species_observed = ifelse(is.na(Epitope.species_observed), "none", Epitope.species_observed)) %>%
  dplyr::mutate(in_vitro_tested_Epitope.species_observed = ifelse(Epitope.species_observed != "none", Epitope.species_observed, in_vitro_tested)) %>%
  dplyr::mutate(in_vitro_tested_Epitope.species_observed = ifelse(in_vitro_tested_Epitope.species_observed == "FALSE", "not tested", in_vitro_tested_Epitope.species_observed)) %>%
  dplyr::mutate(in_vitro_tested_Epitope.species_observed = ifelse(in_vitro_tested_Epitope.species_observed == "TRUE", "*in vitro* tested", in_vitro_tested_Epitope.species_observed)) %>%
  dplyr::mutate(in_vitro_tested_Epitope.species_observed = ifelse(is.na(cl_name), "no TCR", in_vitro_tested_Epitope.species_observed)) %>%
  dplyr::mutate(in_vitro_tested_Epitope.species_observed = ifelse(in_vitro_tested_Epitope.species_observed %in% c("not tested", "no TCR"), "not tested", in_vitro_tested_Epitope.species_observed))

names(SO@misc$label_df) <- tolower(names(SO@misc$label_df))
p_react_tcr <- scexpr::feature_plot2(SO,
                                     pt_size = 0.5,
                                     features = "in_vitro_tested_Epitope.species_observed",
                                     reduction = "tsne",
                                     col_pal_d_args = list(name = epitope_colors),
                                     split_feature = "body_fluid",
                                     name_anno_pos = NULL,
                                     contour_args = list(color = "black", breaks = 0.15, linewidth = 0.2),
                                     contour_feature = "clusters",
                                     facet_grid_row_var = "CD4<sup>+</sup> or CD8<sup>+</sup> T cells",
                                     axes_lim_set = list(x = c(-53,50), y = c(-55, 55))) +
  ggtext::geom_richtext(data = SO@misc$label_df, aes(label = label), size = 3, label.padding = unit(rep(0.01,4),"lines"), label.color = NA) +
  theme(strip.background = element_rect(fill = "grey95", color = "white"),
        legend.position = "right", legend.text = ggtext::element_markdown(),
        strip.text.y = ggtext::element_markdown())


names(SO_CD8_eff@misc$label_df) <- tolower(names(SO_CD8_eff@misc$label_df))
p_react_tcr_CD8_eff <- scexpr::feature_plot2(SO_CD8_eff,
                                             pt_size = 0.5,
                                             features = "in_vitro_tested_Epitope.species_observed",
                                             reduction = "tsne",
                                             col_pal_d_args = list(name = epitope_colors),
                                             split_feature = "body_fluid",
                                             name_anno_pos = NULL,
                                             contour_args = list(color = "black", breaks = 0.2, linewidth = 0.2),
                                             facet_grid_row_var = "CD8_CCL5 transcriptomes",
                                             contour_feature = "clusters",
                                             axes_lim_set = list(x = c(-59,62), y = c(-67, 40))) +
  ggtext::geom_richtext(data = SO_CD8_eff@misc$label_df, aes(label = label), size = 3, label.padding = unit(rep(0.01,4),"lines"), label.color = NA) +
  theme(strip.background = element_rect(fill = "grey95", color = "white"),
        legend.position = "right", legend.text = ggtext::element_markdown(),
        strip.text.y = ggtext::element_markdown(),
        strip.text.x = element_blank())
p_react_tcr_comb <- cowplot::plot_grid(patchwork::wrap_plots(p_react_tcr, p_react_tcr_CD8_eff, guides = "collect", ncol = 1)) #+ patchwork::plot_layout(guides = "collect")) #& theme(legend.box.spacing = margin(0.1)))
Fig3E <- p_react_tcr_comb





