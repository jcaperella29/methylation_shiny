library(shiny)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(reshape2)
library(plotly)
library(pwr)
library(uwot)
library(randomForest)
library(caret)
library(pROC)
library(shiny)


ui <- fluidPage(
  titlePanel("DNA Methylation Analysis"),
  shiny::tags$head(
    shiny::tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
  ),
  
  sidebarLayout(
    sidebarPanel(
  fileInput("file", "Upload Methylation Data (CSV)", accept = ".csv"),
  
  selectInput("region", "Select Region:", choices = c("All", "Island", "Shore", "Shelf")),
  uiOutput("condition_selector"),

  h4("Differential Methylation"),
  actionButton("run_de", "Run DE Analysis"),
  actionButton("run_volcano", "Run Volcano Plot"),

  h4("Visualizations"),
  actionButton("run_boxplot", "Draw Boxplot"),
  actionButton("run_heatmap", "Draw Heatmap"),
  actionButton("run_pca", "Run PCA"),
  actionButton("run_umap", "Run UMAP"),

  h4("Pathway Analysis Options"),
  radioButtons("de_direction", "Select Genes for Pathway Analysis:",
               choices = c("All" = "all", "Hypermethylated (Up)" = "up", "Hypomethylated (Down)" = "down"),
               selected = "all"),
  selectInput("enrichr_db", "Select Enrichment Database:",
              choices = c("KEGG_2021_Human", "Reactome_2022", "GO_Biological_Process_2023"),
              selected = "KEGG_2021_Human"),
  actionButton("run_pathway", "Run Pathway Analysis"),

  h4("Power Analysis"),
  numericInput("power_effect", "Effect size (Cohen's d):", value = 0.5, min = 0.1, step = 0.1),
  numericInput("power_alpha", "Significance level (α):", value = 0.05, min = 0.001, step = 0.005),
  numericInput("power_power", "Desired power (1 - β):", value = 0.8, min = 0.1, max = 0.99, step = 0.05),
  uiOutput("group_selector"),
  actionButton("run_power", "Run Power Analysis"),

  h4("Random Forest"),
  actionButton("run_rf", "Run Random Forest"),

  h4("Download Results"),
  downloadButton("download_de", "Download DE Results"),
  downloadButton("download_pathway", "Download Pathway Results"),
  downloadButton("download_power_table", "Download Power Table")
),

    
    mainPanel(
  tabsetPanel(
    tabPanel("README", uiOutput("readme_text")),
    tabPanel("DE_Results", DT::dataTableOutput("de_table")),
    tabPanel("Volcano Plot", plotlyOutput("volcano_plot")),
    
    tabPanel("Boxplot", plotlyOutput("boxplot")),
    tabPanel("Heatmap", plotlyOutput("heatmap")),
    tabPanel("PCA Plot", plotlyOutput("pca_plot")),
    tabPanel("UMAP Plot", plotlyOutput("umap_plot")),
    
    tabPanel("Pathway Analysis",
             tabsetPanel(
               tabPanel("Tabular Results", DT::dataTableOutput("enrichr_table")),
               tabPanel("Barplot", plotlyOutput("enrichr_barplot"))
             )),
    
    tabPanel("Power Analysis",
             tabsetPanel(
               tabPanel("Tabular Results", DT::dataTableOutput("power_table")),
               tabPanel("Power Plot", plotlyOutput("power_plot"))
             )),
    
    tabPanel("Random Forest",
             tabsetPanel(
               tabPanel("Prediction", DT::dataTableOutput("rf_pred_table"),
                        downloadButton("download_rf_pred", "Download Prediction")),
               tabPanel("Metrics", DT::dataTableOutput("rf_metrics_table"),
                        downloadButton("download_rf_metrics", "Download Metrics")),
               tabPanel("Feature Importance", DT::dataTableOutput("rf_importance_table"),
                        downloadButton("download_rf_importance", "Download Importance"))
             ))
  )
)))

  
  



server <- function(input, output) {
  #reactives
  enrichr_results <- reactiveVal(NULL)
  power_table_result <- reactiveVal(NULL)
  volcano_data <- reactiveVal(NULL)
  umap_data <- reactiveVal(NULL)
  pca_data <- reactiveVal(NULL)
  rf_pred <- reactiveVal(NULL)
  rf_metrics <- reactiveVal(NULL)
  rf_importance <- reactiveVal(NULL)


  # ---- Helper: assign group/condition from sample name by prefix ----
assign_group <- function(sample_names, groups) {
  # sample_names: character vector like c("Tumor1","Normal2")
  # groups: character vector like c("Tumor","Normal")
  sapply(as.character(sample_names), function(s) {
    matched <- groups[sapply(groups, function(g) startsWith(s, g))]
    if (length(matched) > 0) matched[1] else NA_character_
  })
}

  # ---- Enrichr helpers: unify + barplot ----
  fmt_enrich_for_table <- function(res) {
    if (is.null(res) || nrow(res) == 0) return(NULL)
    
    # Make sure these columns exist
    needed <- c("DB", "Term", "Overlap", "P.value", "Adjusted.P.value", "Genes")
    missing <- setdiff(needed, colnames(res))
    for (nm in missing) {
      res[[nm]] <- NA
    }
    
    df <- res[, needed, drop = FALSE]
    df <- df[order(df$Adjusted.P.value, df$P.value), , drop = FALSE]
    rownames(df) <- NULL
    df
  }

  render_enrich_barplot <- function(df, db = NULL, top_n = 10) {
    req(df, nrow(df) > 0)
    
    if (!is.null(db) && "DB" %in% colnames(df)) {
      df <- df[df$DB == db, , drop = FALSE]
    }
    if (nrow(df) == 0) return(NULL)
    
    df <- df[order(df$Adjusted.P.value), , drop = FALSE]
    df <- head(df, top_n)
    df$Term <- factor(df$Term, levels = rev(df$Term))
    
    plot_ly(
      df,
      x = ~-log10(Adjusted.P.value),
      y = ~Term,
      type = "bar",
      orientation = "h",
      text = ~paste0("FDR: ", signif(Adjusted.P.value, 3),
                     "<br>Genes: ", Genes)
    ) %>%
      layout(
        title = paste("Top Enriched Pathways:", if (!is.null(db)) db else ""),
        xaxis = list(title = "-log10(FDR)"),
        yaxis = list(title = "")
      )
  }

  
  data <- reactive({
    req(input$file)
    df <- read.csv(input$file$datapath)
    df
  })
  
  filtered_data <- reactive({
  req(input$condition)   # IMPORTANT: since we use selected groups
  df <- data()

  if (input$region != "All") {
    df <- df[df$Region == input$region, ]
  }

  df_melt <- melt(df,
                  id.vars = c("Gene", "Region"),
                  variable.name = "Sample",
                  value.name = "Methylation")

  # Use the SAME rule as everywhere else
  df_melt$Condition <- assign_group(df_melt$Sample, input$condition)

  # Keep only selected groups
  df_melt <- df_melt[!is.na(df_melt$Condition) & df_melt$Condition %in% input$condition, ]
  df_melt
})

  output$boxplot <- renderPlotly({
    req(input$run_boxplot)
    req(input$condition)
    showNotification("Boxplot generation beginning!", type = "message")
    
    isolate({
      df <- data()
      
      df_melt <- melt(df, id.vars = c("Gene", "Region"), variable.name = "Sample", value.name = "Methylation")
      
      df_melt$Condition <- sapply(as.character(df_melt$Sample), function(s) {
        matched <- input$condition[sapply(input$condition, function(g) startsWith(as.character(s), g))]
        if (length(matched) > 0) matched[1] else NA
      })
      
      df_melt <- df_melt[df_melt$Condition %in% input$condition, ]
      
      p <- ggplot(df_melt, aes(x = Region, y = Methylation, fill = Condition)) +
        geom_boxplot() +
        theme_minimal() +
        ggtitle("Methylation by Region and Condition")
      
      showNotification("Boxplot generation complete!", type = "message")
      ggplotly(p)
    })
  })
  
  
  
  output$heatmap <- renderPlotly({
    req(input$run_heatmap)
    req(input$condition)
    showNotification("Heatmap generation beginning!", type = "message")
    
    isolate({
      df <- data()
      
      # Find relevant columns
      sample_cols <- colnames(df)[!(colnames(df) %in% c("Gene", "Region"))]
      selected_cols <- sample_cols[sapply(sample_cols, function(s) any(startsWith(s, input$condition)))]
      
      if (length(selected_cols) == 0) return(NULL)
      
      meth_matrix <- df[, selected_cols, drop = FALSE]
      rownames(meth_matrix) <- make.unique(df$Gene)
      meth_matrix <- apply(meth_matrix, 2, as.numeric)
      
      # Remove problematic rows
      meth_matrix <- meth_matrix[complete.cases(meth_matrix) &
                                   apply(meth_matrix, 1, function(x) all(is.finite(x))), , drop = FALSE]
      if (nrow(meth_matrix) == 0 || ncol(meth_matrix) == 0) return(NULL)
      
      # Cluster
      row_clust <- hclust(dist(meth_matrix))
      col_clust <- hclust(dist(t(meth_matrix)))
      meth_matrix <- meth_matrix[row_clust$order, col_clust$order]
      
      p <- plot_ly(
        z = meth_matrix,
        x = colnames(meth_matrix),
        y = rownames(meth_matrix),
        type = "heatmap",
        colorscale = "viridis"
      ) %>%
        layout(title = "Clustered Methylation Heatmap")
      
      showNotification("Heatmap generation completed!", type = "message")
      return(p)
    })
  })
  
  
  de_results <- reactiveVal(NULL)
  
  observeEvent(input$run_de, {
    df <- data()
    req(input$condition)
    if (length(input$condition) != 2) {
      showNotification("Please select exactly two groups for DE analysis.", type = "error")
      return()
    }
    
    showNotification("Differential methylation analysis beginning!", type = "message")
    
    selected_groups <- input$condition
    
    # Match columns for selected groups
    all_cols <- colnames(df)[!(colnames(df) %in% c("Gene", "Region"))]
    sample_labels <- sapply(all_cols, function(x) {
      matched <- selected_groups[which.max(sapply(selected_groups, function(g) startsWith(x, g)))]
      if (is.na(matched)) return(NA) else return(matched)
    })
    
    valid_cols <- !is.na(sample_labels)
    meth_matrix <- df[, all_cols[valid_cols], drop = FALSE]
    sample_labels <- sample_labels[valid_cols]
    rownames(meth_matrix) <- df$Gene
    
    # Clean rows
    meth_matrix <- meth_matrix[complete.cases(meth_matrix), , drop = FALSE]
    
    # T-test per probe
    pvals <- apply(meth_matrix, 1, function(x) {
      tryCatch(t.test(x ~ sample_labels)$p.value, error = function(e) NA)
    })
    
    logFC <- rowMeans(meth_matrix[, sample_labels == selected_groups[1], drop = FALSE], na.rm = TRUE) -
      rowMeans(meth_matrix[, sample_labels == selected_groups[2], drop = FALSE], na.rm = TRUE)
    
    padj <- p.adjust(pvals, method = "fdr")
    
    result <- data.frame(
      Gene = rownames(meth_matrix),
      logFC = logFC,
      pvalue = pvals,
      padj = padj,
      row.names = NULL
    ) %>% arrange(padj)
    
    de_results(result)
    showNotification("Differential methylation analysis complete!", type = "message")
  })
  
 
  observeEvent(input$run_pathway, {
    req(de_results())
    showNotification("Pathway analysis starting...", type = "message")
    
    # Filter DE genes based on FDR and direction
    sig <- de_results() %>% dplyr::filter(padj <= 0.05)
    
    direction <- input$de_direction
    if (direction == "up") {
      sig <- sig %>% dplyr::filter(logFC > 0)
    } else if (direction == "down") {
      sig <- sig %>% dplyr::filter(logFC < 0)
    }
    
    if (nrow(sig) == 0) {
      showNotification("No significant genes in selected direction", type = "warning")
      enrichr_results(NULL)
      return()
    }
    
    # Map CpG probes to gene symbols
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    probe_annot <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    sig_genes <- probe_annot[sig$Gene, "UCSC_RefGene_Name"]
    sig_genes <- unique(unlist(strsplit(sig_genes, ";")))
    sig_genes <- na.omit(sig_genes)
    sig_genes <- sig_genes[sig_genes != ""]
    
    if (length(sig_genes) < 1) {
      showNotification("No valid gene symbols mapped from probes", type = "warning")
      enrichr_results(NULL)
      return()
    }
    
    # Run Enrichr on all three DBs, then unify
    library(enrichR)
    databases <- c("KEGG_2021_Human", "Reactome_2022", "GO_Biological_Process_2023")
    raw_list <- enrichr(sig_genes, databases)
    
    combined <- do.call(rbind, lapply(names(raw_list), function(db) {
      df <- raw_list[[db]]
      if (is.null(df) || nrow(df) == 0) return(NULL)
      df$DB <- db
      df
    }))
    
    combined <- fmt_enrich_for_table(combined)
    enrichr_results(combined)
    
    showNotification("Pathway analysis completed!", type = "message")
  })

    output$enrichr_table <- DT::renderDataTable({
    req(enrichr_results())
    df <- enrichr_results()
    db <- input$enrichr_db
    
    if ("DB" %in% colnames(df)) {
      df <- df[df$DB == db, , drop = FALSE]
    }
    if (nrow(df) == 0) {
      showNotification("Selected database has no results", type = "warning")
      return(NULL)
    }
    
    DT::datatable(df,
                  options = list(pageLength = 10),
                  rownames = FALSE)
  })

  output$de_table <- DT::renderDataTable({
    req(de_results())
    DT::datatable(de_results(), 
                  options = list(pageLength = 10),
                  rownames = FALSE)
  })
  output$download_de <- downloadHandler(
    filename = function() {
      paste0("differential_methylation_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(de_results())
      write.csv(de_results(), file, row.names = FALSE)
    }
  )
  output$download_pathway <- downloadHandler(
    filename = function() {
      paste0("pathway_analysis_", input$enrichr_db, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(enrichr_results())
      df <- enrichr_results()
      db <- input$enrichr_db
      
      if ("DB" %in% colnames(df)) {
        df <- df[df$DB == db, , drop = FALSE]
      }
      if (nrow(df) == 0) {
        showNotification("Selected pathway database has no results.", type = "error")
        write.csv(data.frame(), file, row.names = FALSE)
        return()
      }
      
      write.csv(df, file, row.names = FALSE)
    }
  )

  output$enrichr_barplot <- renderPlotly({
    req(enrichr_results())
    df <- enrichr_results()
    render_enrich_barplot(df, db = input$enrichr_db, top_n = 10)
  })

  
  group_choices <- reactive({
    df <- data()
    sample_cols <- colnames(df)[!(colnames(df) %in% c("Gene", "Region"))]
    groups <- unique(gsub("[0-9]+$", "", sample_cols))  # remove trailing numbers (e.g. Tumor1 → Tumor)
    sort(groups)
  })

  observeEvent(input$run_power, {
    req(data())
    req(input$condition)
    
    selected <- input$condition
    if (length(selected) != 2) {
      showNotification("Please select exactly two groups for power analysis.", type = "error")
      return()
    }
    
    df <- data()
    sample_cols <- colnames(df)[!(colnames(df) %in% c("Gene", "Region"))]
    
    # Match columns to selected groups
    sample_labels <- sapply(as.character(sample_cols), function(name) {
      matched <- selected[sapply(selected, function(g) startsWith(name, g))]
      if (length(matched) > 0) matched[1] else NA
    })
    
    valid_idx <- !is.na(sample_labels)
    sample_labels <- sample_labels[valid_idx]
    sample_cols <- sample_cols[valid_idx]
    
    if (length(unique(sample_labels)) < 2) {
      showNotification("Unable to find samples for both selected groups.", type = "error")
      return()
    }
    
    n_per_group <- min(table(sample_labels))
    
    # Power analysis calculation
    pwr_result <- pwr.t.test(n = n_per_group,
                             d = input$power_effect,
                             sig.level = input$power_alpha,
                             type = "two.sample",
                             alternative = "two.sided")
    
    # Output table
    power_table <- data.frame(
      Parameter = c("Effect size (d)", "Significance level (α)", "Sample size per group", "Estimated power"),
      Value = c(input$power_effect, input$power_alpha, n_per_group, round(pwr_result$power, 3))
    )
    
    output$power_table <- DT::renderDataTable({
      DT::datatable(power_table, rownames = FALSE, options = list(dom = 't'))
    
      
    })
    power_table_result(power_table)  # Store for download
    
    
    # Output power curve
    sample_sizes <- seq(5, 100, by = 1)
    powers <- sapply(sample_sizes, function(n) {
      pwr.t.test(n = n, d = input$power_effect,
                 sig.level = input$power_alpha,
                 type = "two.sample")$power
    })
    
    output$power_plot <- renderPlotly({
      plot_ly() %>%
        add_trace(
          x = sample_sizes,
          y = powers,
          type = "scatter",
          mode = "lines",
          name = "Power Curve",
          line = list(color = "darkblue", width = 3)
        ) %>%
        add_trace(
          x = c(n_per_group),
          y = c(pwr_result$power),
          type = "scatter",
          mode = "markers+text",
          text = paste0("N = ", n_per_group),
          textposition = "top right",
          marker = list(color = "green", size = 10),
          name = "Your Sample Size"
        ) %>%
        add_trace(
          x = sample_sizes,
          y = rep(input$power_power, length(sample_sizes)),
          type = "scatter",
          mode = "lines",
          line = list(color = "red", dash = "dash"),
          name = "Target Power"
        ) %>%
        layout(
          title = paste("Power Curve (", selected[1], " vs ", selected[2], ")"),
          xaxis = list(title = "Sample Size per Group"),
          yaxis = list(title = "Power"),
          showlegend = TRUE
        )
    })
    
    
    
    showNotification("Power analysis completed!", type = "message")
  })
  
  output$condition_selector <- renderUI({
    req(data())  # ensure file is uploaded
    df <- data()
    
    # Get only sample columns
    sample_cols <- colnames(df)[!(colnames(df) %in% c("Gene", "Region"))]
    
    # Strip trailing digits: Tumor1, Tumor2 → Tumor
    group_names <- unique(gsub("[0-9]+$", "", sample_cols))
    
    checkboxGroupInput("condition", "Select Groups for Analysis:",
                       choices = group_names,
                       selected = head(group_names, 2))
  })
  output$download_power_table <- downloadHandler(
    filename = function() {
      paste0("power_analysis_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(power_table_result())
      write.csv(power_table_result(), file, row.names = FALSE)
    }
  )
  observeEvent(input$run_volcano, {
    req(de_results())
    df <- de_results()
    
    df$neglogP <- -log10(pmax(df$pvalue, 1e-300))
    df$Significant <- ifelse(df$padj <= 0.05, "FDR ≤ 0.05", "Not Sig")
    
    volcano_data(df)
    showNotification("Volcano plot generated!", type = "message")
  })
  
  
  
  output$volcano_plot <- renderPlotly({
    req(volcano_data())
    df <- volcano_data()
    
    plot_ly(
      data = df,
      x = ~logFC,
      y = ~neglogP,
      text = ~paste("Gene:", Gene,
                    "<br>logFC:", round(logFC, 3),
                    "<br>p-adj:", signif(padj, 3)),
      color = ~Significant,
      colors = c("gray", "red"),
      type = "scatter",
      mode = "markers"
    ) %>%
      layout(
        title = "Volcano Plot: logFC vs -log10(p-value)",
        xaxis = list(title = "log2 Fold Change"),
        yaxis = list(title = "-log10(p-value)"),
        legend = list(x = 0.01, y = 0.99)
      )
  })
  
  
  

  
  observeEvent(input$run_umap, {
    req(data())
    req(input$condition)
    showNotification("generating UMAP",type= "message")
    df <- data()
    sample_cols <- colnames(df)[!(colnames(df) %in% c("Gene", "Region"))]
    selected_cols <- sample_cols[sapply(sample_cols, function(s) any(startsWith(s, input$condition)))]
    
    if (length(selected_cols) < 2) {
      showNotification("Not enough samples for UMAP.", type = "error")
      return()
    }
    
    meth_matrix <- df[, selected_cols, drop = FALSE]
    meth_matrix <- apply(meth_matrix, 2, as.numeric)
    rownames(meth_matrix) <- df$Gene
    meth_matrix <- meth_matrix[complete.cases(meth_matrix), , drop = FALSE]
    
    # Transpose: rows = samples, cols = probes
    umap_result <- uwot::umap(t(meth_matrix), n_neighbors = 15, min_dist = 0.1, metric = "euclidean")
    
    umap_df <- as.data.frame(umap_result)
    colnames(umap_df) <- c("UMAP1", "UMAP2")
    umap_df$Sample <- rownames(umap_df)
    umap_df$Group <- sapply(umap_df$Sample, function(s) {
      matched <- input$condition[sapply(input$condition, function(g) startsWith(s, g))]
      if (length(matched) > 0) matched[1] else NA
    })
    
    umap_data(umap_df)
    showNotification("UMAP completed!", type = "message")
  })
  output$umap_plot <- renderPlotly({
    req(umap_data())
    df <- umap_data()
    
    plot_ly(
      data = df,
      x = ~UMAP1,
      y = ~UMAP2,
      type = "scatter",
      mode = "markers",
      color = ~Group,
      text = ~Sample,
      marker = list(size = 10)
    ) %>%
      layout(
        title = "UMAP of Samples",
        xaxis = list(title = "UMAP1"),
        yaxis = list(title = "UMAP2")
      )
  })
  
  
  observeEvent(input$run_pca, {
    req(data())
    req(input$condition)
    showNotification("Generating PCA plot",type = "message")
    df <- data()
    
    # Extract columns corresponding to selected groups
    sample_cols <- colnames(df)[!(colnames(df) %in% c("Gene", "Region"))]
    selected_cols <- sample_cols[sapply(sample_cols, function(s) any(startsWith(s, input$condition)))]
    
    if (length(selected_cols) < 2) {
      showNotification("Not enough samples for PCA.", type = "error")
      return()
    }
    
    meth_matrix <- df[, selected_cols, drop = FALSE]
    meth_matrix <- apply(meth_matrix, 2, as.numeric)
    rownames(meth_matrix) <- df$Gene
    meth_matrix <- meth_matrix[complete.cases(meth_matrix), , drop = FALSE]
    
    # Run PCA on transposed matrix (samples = rows)
    pca_result <- prcomp(t(meth_matrix), scale. = TRUE)
    
    pca_df <- as.data.frame(pca_result$x[, 1:2])
    colnames(pca_df) <- c("PC1", "PC2")
    pca_df$Sample <- rownames(pca_df)
    pca_df$Group <- sapply(pca_df$Sample, function(s) {
      matched <- input$condition[sapply(input$condition, function(g) startsWith(s, g))]
      if (length(matched) > 0) matched[1] else NA
    })
    
    pca_data(pca_df)
    showNotification("PCA completed!", type = "message")
  })
  output$pca_plot <- renderPlotly({
    req(pca_data())
    df <- pca_data()
    
    plot_ly(
      data = df,
      x = ~PC1,
      y = ~PC2,
      type = "scatter",
      mode = "markers",
      color = ~Group,
      text = ~Sample,
      marker = list(size = 10)
    ) %>%
      layout(
        title = "PCA Plot of Samples",
        xaxis = list(title = "PC1"),
        yaxis = list(title = "PC2")
      )
  })
  observeEvent(input$run_rf, {
    req(data())
    req(de_results())
    req(input$condition)
    if (length(input$condition) != 2) {
      showNotification("Select exactly two groups for classification.", type = "error")
      return()
    }
    
    
    showNotification("Running Random Forest classfication", type = "message")
    # Step 1: get significant probes
    sig_probes <- de_results() %>% filter(padj <= 0.05) %>% pull(Gene)
    if (length(sig_probes) < 2) {
      showNotification("No significant probes available for Random Forest.", type = "error")
      return()
    }
    
    df <- data()
    sample_cols <- colnames(df)[!(colnames(df) %in% c("Gene", "Region"))]
    
    # Select only matching samples for selected groups
    sample_labels <- sapply(sample_cols, function(s) {
      matched <- input$condition[sapply(input$condition, function(g) startsWith(s, g))]
      if (length(matched) > 0) matched[1] else NA
    })
    
    valid_idx <- !is.na(sample_labels)
    sample_cols <- sample_cols[valid_idx]
    sample_labels <- sample_labels[valid_idx]
    
    # Build data matrix
    meth_matrix <- df[df$Gene %in% sig_probes, sample_cols, drop = FALSE]
    meth_matrix <- apply(meth_matrix, 2, as.numeric)
    meth_matrix <- t(meth_matrix)  # samples in rows
    colnames(meth_matrix) <- make.unique(df[df$Gene %in% sig_probes, "Gene"])
    
    # Convert to data.frame for training
    rf_df <- data.frame(meth_matrix)
    rf_df$Group <- as.factor(sample_labels)
    
    # Train/Test split
    set.seed(123)
    idx <- caret::createDataPartition(rf_df$Group, p = 0.7, list = FALSE)
    train <- rf_df[idx, ]
    test <- rf_df[-idx, ]
    
    # Train Random Forest
    rf_model <- randomForest(Group ~ ., data = train, importance = TRUE)
    
    # Predict
    pred <- predict(rf_model, newdata = test, type = "response")
    probs <- predict(rf_model, newdata = test, type = "prob")[, 2]  # for AUC
    
    # Metrics
    actual <- test$Group
    pred_factor <- factor(pred, levels = levels(actual))
    cm <- caret::confusionMatrix(pred_factor, actual)
    roc <- pROC::roc(actual, probs, levels = rev(levels(actual)))
    
    # Store prediction table
    rf_pred(data.frame(Sample = rownames(test), Actual = actual, Predicted = pred))
    
    # Store metrics
    rf_metrics(data.frame(
      Accuracy = round(cm$overall["Accuracy"], 3),
      Sensitivity = round(cm$byClass["Sensitivity"], 3),
      Specificity = round(cm$byClass["Specificity"], 3),
      AUC = round(auc(roc), 3)
    ))
    
    # Store importance
    imp <- randomForest::importance(rf_model)
    imp_df <- data.frame(Feature = rownames(imp), MeanDecreaseGini = imp[, "MeanDecreaseGini"])
    rf_importance(imp_df[order(-imp_df$MeanDecreaseGini), ])
    
    showNotification("Random Forest analysis complete!", type = "message")
  })
  output$rf_pred_table <- DT::renderDataTable({
    req(rf_pred())
    DT::datatable(rf_pred(), rownames = FALSE)
  })
  
  output$rf_metrics_table <- DT::renderDataTable({
    req(rf_metrics())
    DT::datatable(rf_metrics(), rownames = FALSE)
  })
  
  output$rf_importance_table <- DT::renderDataTable({
    req(rf_importance())
    DT::datatable(rf_importance(), rownames = FALSE)
  })
  
  output$download_rf_pred <- downloadHandler(
    filename = function() paste0("rf_prediction_", Sys.Date(), ".csv"),
    content = function(file) write.csv(rf_pred(), file, row.names = FALSE)
  )
  
  output$download_rf_metrics <- downloadHandler(
    filename = function() paste0("rf_metrics_", Sys.Date(), ".csv"),
    content = function(file) write.csv(rf_metrics(), file, row.names = FALSE)
  )
  
  output$download_rf_importance <- downloadHandler(
    filename = function() paste0("rf_importance_", Sys.Date(), ".csv"),
    content = function(file) write.csv(rf_importance(), file, row.names = FALSE)
  )
  
  
  output$readme_text <- renderText({
    readLines("README.txt")
  })
  
  
}

shinyApp(ui = ui, server = server)





