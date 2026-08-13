#### The App

# ============================================================
# Single Cell Explorer
# ============================================================

library(shiny)
library(Seurat)
library(ggplot2)
library(dplyr)
library(bslib)
library(pals)
library(Polychrome)

# ============================================================
# 1. Dataset configuration
# ============================================================

dataset_files <- c(
  "Fig 2 scRNAseq of CNS Macrophage Niches" = "web_data/fig2.rds",
  "Fig 5 scRNAseq of sdCAMs" = "web_data/fig5_cams.rds",
  "Fig 5 scRNAseq of Microglia" = "web_data/fig5_mg.rds",
  "Fig 5 scRNAseq of subdural immune cells" = "web_data/fig5_full.rds",
  "Fig 6 scATACseq of CNS Macrophages and Blood Monocytes (gene activity matrix only)" = "web_data/fig6.rds"
)

fig2_cell_subtype_colors <- readRDS(
  "colors/fig2_cell_subtype_colors.rds"
)

# ============================================================
# Color helper
# ============================================================

get_group_colors <- function(obj, metadata, dataset) {
  
  # Get the groups in the same order used for plotting
  groups <- levels(
    droplevels(
      factor(obj@meta.data[[metadata]])
    )
  )
  
  # ----------------------------------------------------------
  # Fig 2 cell subtype: use custom palette
  # ----------------------------------------------------------
  
  if (
    dataset == "Fig 2 scRNAseq of CNS Macrophage Niches" &&
    metadata == "cell_subtype_final"
  ) {
    
    # Check that every group has a color
    missing_colors <- setdiff(
      groups,
      names(fig2_cell_subtype_colors)
    )
    
    if (length(missing_colors) > 0) {
      
      stop(
        paste(
          "Missing colors in Fig 2 palette for:",
          paste(missing_colors, collapse = ", ")
        )
      )
      
    }
    
    return(
      fig2_cell_subtype_colors[groups]
    )
  }
  
  
  # ----------------------------------------------------------
  # All other metadata: generate a consistent palette
  # ----------------------------------------------------------
  
  seed_colors <- c(
    "#1034A6",
    "#C30F16",
    "#CCBB44", 
    "#66CCEE",
    "#AA3377",
    "#BBBBBB",
    "#008000"
  )
  
  colors <- createPalette(
    length(groups),
    seedcolors = seed_colors
  )
  
  names(colors) <- groups
  
  colors
  
}

# ============================================================
# Condition color helper
# ============================================================

get_condition_colors <- function(obj) {
  
  conditions <- levels(
    droplevels(
      factor(obj@meta.data[["condition"]])
    )
  )
  
  colors <- createPalette(
    length(conditions),
    seedcolors = c(
      "#1034A6",
      "#C30F16",
      "#CCBB44",
      "#66CCEE",
      "#AA3377",
      "#BBBBBB",
      "#008000"
    )
  )
  
  names(colors) <- conditions
  
  colors
}


# ============================================================
# 2. UI
# ============================================================

ui <- page_sidebar(
  
  title = "Single Cell Explorer for Fliegauf et al. Nature Immunology 2026",
  
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly"
  ),
  
  sidebar = sidebar(
    
    width = 320,
    
    h4("Dataset"),
    
    selectInput(
      inputId = "dataset",
      label = NULL,
      choices = names(dataset_files),
      selected = names(dataset_files)[2]
    ),
    
    hr(),
    
    h4("UMAP coloring"),
    
    selectInput(
      inputId = "metadata",
      label = "Color cells by:",
      choices = NULL
    ),
    
    hr(),
    
    h4("Gene"),
    
    selectizeInput(
      inputId = "gene",
      label = "Search for a gene:",
      choices = NULL,
      selected = NULL,
      
      options = list(
        placeholder = "Type a gene name..."
      )
    ),
    
    hr(),
    
    checkboxInput(
      inputId = "raster",
      label = "Rasterize cells",
      value = TRUE
    ),
    
    hr(),
    
    strong("Current dataset"),
    
    textOutput("cell_count"),
    
    textOutput("gene_count"),
    
    br(),
    
    strong("Selected gene"),
    
    textOutput("percent_expressing"),
    
    textOutput("mean_expression"),
    
    textOutput("median_expression")
  ),
  
  
  # ==========================================================
  # Main content
  # ==========================================================
  
  navset_tab(
    
    # --------------------------------------------------------
    # UMAP TAB
    # --------------------------------------------------------
    
    nav_panel(
      
      "UMAP",
      
      br(),
      
      layout_columns(
        
        card(
          card_header("Cells colored by metadata"),
          
          plotOutput(
            outputId = "umap",
            height = "650px"
          )
        ),
        
        card(
          card_header("Gene expression"),
          
          plotOutput(
            outputId = "feature",
            height = "650px"
          )
        ),
        
        col_widths = c(6, 6)
      )
    ),
    
    
    # --------------------------------------------------------
    # EXPRESSION TAB
    # --------------------------------------------------------
    
    nav_panel(
      
      "Expression",
      
      br(),
      
      card(
        card_header("Expression by metadata"),
        
        plotOutput(
          outputId = "violin",
          height = "700px"
        )
      ),
      
      br(),
      
      card(
        card_header("Mean expression"),
        
        plotOutput(
          outputId = "barplot",
          height = "700px"
        )
      )
    ),   
    
    # --------------------------------------------------------
    # DATASET TAB
    # --------------------------------------------------------
    
    nav_panel(
      
      "Dataset information",
      
      br(),
      
      layout_columns(
        
        card(
          card_header("Dataset information"),
          
          verbatimTextOutput(
            outputId = "dataset_info"
          )
        ),
        
        card(
          card_header("Selected metadata"),
          
          tableOutput(
            outputId = "metadata_summary"
          )
        ),
        
        col_widths = c(6, 6)
      )
    )
  )
)


# ============================================================
# 3. Server
# ============================================================

server <- function(input, output, session) {
  
  
  # ==========================================================
  # Load selected dataset
  # ==========================================================
  
  current_object <- reactive({
    
    req(input$dataset)
    
    path <- dataset_files[[input$dataset]]
    
    validate(
      need(
        file.exists(path),
        paste("Could not find:", path)
      )
    )
    
    obj <- readRDS(path)
    
    obj
    
  })
  
  
  # ==========================================================
  # Update metadata choices
  # ==========================================================
  
  observeEvent(
    
    current_object(),
    
    {
      
      obj <- current_object()
      
      metadata <- colnames(obj@meta.data)
      
      updateSelectInput(
        session = session,
        inputId = "metadata",
        choices = metadata,
        selected = metadata[1]
      )
      
    },
    
    ignoreInit = FALSE
  )
  
  
  # ==========================================================
  # Update gene choices
  # ==========================================================
  
  observeEvent(
    
    current_object(),
    
    {
      
      obj <- current_object()
      
      genes <- if ("RNA" %in% Assays(obj)) {
        
        rownames(obj[["RNA"]])
        
      } else {
        
        rownames(obj[["RNA_called_peaks"]])
        
      }
      
      default_gene <- if (
        "CD3D" %in% genes
      ) {
        "CD3D"
      } else {
        genes[1]
      }
      
      updateSelectizeInput(
        
        session = session,
        
        inputId = "gene",
        
        choices = genes,
        
        selected = default_gene,
        
        server = TRUE
      )
      
    },
    
    ignoreInit = FALSE
  )
  
  
  # ==========================================================
  # UMAP
  # ==========================================================
  
  output$umap <- renderPlot({
    
    obj <- current_object()
    
    req(input$metadata)
    
    
    # --------------------------------------------------------
    # Get colors for current metadata
    # --------------------------------------------------------
    
    group_colors <- get_group_colors(
      obj = obj,
      metadata = input$metadata,
      dataset = input$dataset
    )
    
    
    # --------------------------------------------------------
    # UMAP
    # --------------------------------------------------------
    
    p <- DimPlot(
      
      object = obj,
      
      reduction = "umap",
      
      group.by = input$metadata,
      
      cols = group_colors,
      
      pt.size = 2,
      
      alpha = 0.5,
      
      raster = input$raster,
      
      shuffle = TRUE,
      
      seed = 123
      
    ) +
      
      labs(
        title = input$dataset,
        subtitle = paste(
          "Colored by:",
          input$metadata
        )
      ) +
      
      theme_classic(base_size = 14) +
      
      NoAxes() +
      
      theme(
        legend.position = "bottom",
        legend.text = element_text(size = 12),
        legend.key.size = grid::unit(0.8, "cm")
      ) +
      
      guides(
        color = guide_legend(
          nrow = if (
            input$dataset == "Fig 2 scRNAseq of CNS Macrophage Niches" &&
            input$metadata == "cell_subtype_final"
          ) {
            8
          } else {
            3
          },
          byrow = TRUE,
          override.aes = list(
            size = 6,
            alpha = 1
          )
        )
      )
    
    
    # --------------------------------------------------------
    # Make UMAP cells semi-transparent
    # --------------------------------------------------------
    
    p$layers[[1]]$aes_params$alpha <- 0.65
    
    
    p
    
  })
  
  
  # ==========================================================
  # Feature Plot
  # ==========================================================
  
  output$feature <- renderPlot({
    
    obj <- current_object()
    
    req(input$gene)
    
    validate(
      
      need(
        
        input$gene %in% rownames(obj[["RNA"]]),
        
        paste(
          "Gene",
          input$gene,
          "was not found in this dataset."
        )
        
      )
      
    )
    
    FeaturePlot(
      
      object = obj,
      
      features = input$gene,
      
      reduction = "umap",
      
      slot = "data",
      
      pt.size = 2,
      
      raster = input$raster,
      
      order = TRUE
      
    ) +
      
      labs(
        title = input$gene,
        subtitle = "RNA expression"
      ) +
      
      theme_classic(base_size = 14) +
      NoAxes()
    
  })
  
  
  # ==========================================================
  # Violin plot
  # ==========================================================
  
  output$violin <- renderPlot({
    
    obj <- current_object()
    
    req(input$gene)
    req(input$metadata)
    
    validate(
      need(
        input$gene %in% rownames(obj[["RNA"]]),
        "Gene not found."
      )
    )
    
    # --------------------------------------------------------
    # Special case: condition
    # Group by cell_type and split by condition
    # --------------------------------------------------------
    
    if (input$metadata == "condition") {
      
      condition_colors <- get_condition_colors(obj)
      
      VlnPlot(
        object = obj,
        features = input$gene,
        layer = "data",
        group.by = "cell_type",
        split.by = "condition",
        cols = condition_colors,
        pt.size = 0
      ) +
        labs(
          title = input$gene,
          x = "cell_type",
          y = "Normalized expression"
        ) +
        theme_classic(base_size = 13) +
        theme(
          axis.text.x = element_text(
            angle = 45,
            hjust = 1
          ),
          legend.position = "bottom"
        )
      
    } else {
      
      # ------------------------------------------------------
      # Default behavior for all other metadata
      # ------------------------------------------------------
      
      group_colors <- get_group_colors(
        obj = obj,
        metadata = input$metadata,
        dataset = input$dataset
      )
      
      VlnPlot(
        object = obj,
        features = input$gene,
        layer = "data",
        group.by = input$metadata,
        cols = group_colors,
        pt.size = 0
      ) +
        labs(
          title = input$gene,
          x = input$metadata,
          y = "Normalized expression"
        ) +
        theme_classic(base_size = 13) +
        theme(
          axis.text.x = element_text(
            angle = 45,
            hjust = 1
          ),
          legend.position = "none"
        )
    }
  })
  
  
  # ==========================================================
  # Mean expression bar plot
  # ==========================================================
  
  output$barplot <- renderPlot({
    
    obj <- current_object()
    
    req(input$gene)
    req(input$metadata)
    
    validate(
      need(
        input$gene %in% rownames(obj[["RNA"]]),
        "Gene not found."
      )
    )
    
    # --------------------------------------------------------
    # Special case: condition
    # Group by cell_type and split by condition
    # --------------------------------------------------------
    
    if (input$metadata == "condition") {
      
      condition_colors <- get_condition_colors(obj)
      
      expression_data <- FetchData(
        object = obj,
        vars = c(
          input$gene,
          "cell_type",
          "condition"
        ),
        layer = "data"
      )
      
      colnames(expression_data) <- c(
        "expression",
        "cell_type",
        "condition"
      )
      
      summary_data <- expression_data %>%
        group_by(
          cell_type,
          condition
        ) %>%
        summarise(
          mean_expression = mean(
            expression,
            na.rm = TRUE
          ),
          .groups = "drop"
        )
      
      ggplot(
        summary_data,
        aes(
          x = cell_type,
          y = mean_expression,
          fill = condition
        )
      ) +
        geom_col(
          position = position_dodge()
        ) +
        scale_fill_manual(
          values = condition_colors,
          drop = FALSE
        ) +
        labs(
          title = paste(
            "Mean",
            input$gene,
            "expression"
          ),
          x = "cell_type",
          y = "Mean normalized expression",
          fill = "condition"
        ) +
        theme_classic(base_size = 13) +
        theme(
          axis.text.x = element_text(
            angle = 45,
            hjust = 1
          ),
          legend.position = "bottom"
        )
      
    } else {
      
      # ------------------------------------------------------
      # Default behavior for all other metadata
      # ------------------------------------------------------
      
      group_colors <- get_group_colors(
        obj = obj,
        metadata = input$metadata,
        dataset = input$dataset
      )
      
      expression_data <- FetchData(
        object = obj,
        vars = c(
          input$gene,
          input$metadata
        ),
        layer = "data"
      )
      
      colnames(expression_data) <- c(
        "expression",
        "group"
      )
      
      summary_data <- expression_data %>%
        group_by(group) %>%
        summarise(
          mean_expression = mean(
            expression,
            na.rm = TRUE
          ),
          median_expression = median(
            expression,
            na.rm = TRUE
          ),
          cells = n(),
          percent_expressing =
            mean(expression > 0) * 100,
          .groups = "drop"
        )
      
      summary_data$group <- factor(
        summary_data$group,
        levels = names(group_colors)
      )
      
      ggplot(
        summary_data,
        aes(
          x = group,
          y = mean_expression,
          fill = group
        )
      ) +
        geom_col() +
        scale_fill_manual(
          values = group_colors,
          drop = FALSE
        ) +
        labs(
          title = paste(
            "Mean",
            input$gene,
            "expression"
          ),
          x = input$metadata,
          y = "Mean normalized expression"
        ) +
        theme_classic(base_size = 13) +
        theme(
          axis.text.x = element_text(
            angle = 45,
            hjust = 1
          ),
          legend.position = "none"
        )
    }
  })
  
  
  # ==========================================================
  # Cell count
  # ==========================================================
  
  output$cell_count <- renderText({
    
    paste(
      "Cells:",
      format(
        ncol(current_object()),
        big.mark = ","
      )
    )
    
  })
  
  
  # ==========================================================
  # Gene count
  # ==========================================================
  
  output$gene_count <- renderText({
    
    paste(
      "Genes:",
      format(
        nrow(current_object()),
        big.mark = ","
      )
    )
    
  })
  
  
  # ==========================================================
  # Gene statistics
  # ==========================================================
  
  gene_statistics <- reactive({
    
    obj <- current_object()
    
    req(input$gene)
    
    validate(
      
      need(
        
        input$gene %in% rownames(obj[["RNA"]]),
        
        "Gene not found."
        
      )
      
    )
    
    values <- FetchData(
      
      object = obj,
      
      vars = input$gene,
      
      layer = "data"
      
    )[[input$gene]]
    
    
    list(
      
      percent_expressing =
        mean(values > 0) * 100,
      
      mean =
        mean(values, na.rm = TRUE),
      
      median =
        median(values, na.rm = TRUE)
      
    )
    
  })
  
  
  output$percent_expressing <- renderText({
    
    stats <- gene_statistics()
    
    paste0(
      "% expressing: ",
      round(
        stats$percent_expressing,
        1
      ),
      "%"
    )
    
  })
  
  
  output$mean_expression <- renderText({
    
    stats <- gene_statistics()
    
    paste0(
      "Mean: ",
      round(
        stats$mean,
        3
      )
    )
    
  })
  
  
  output$median_expression <- renderText({
    
    stats <- gene_statistics()
    
    paste0(
      "Median: ",
      round(
        stats$median,
        3
      )
    )
    
  })
  
  
  # ==========================================================
  # Dataset information
  # ==========================================================
  
  output$dataset_info <- renderPrint({
    
    obj <- current_object()
    
    cat(
      "Dataset:",
      input$dataset,
      "\n\n"
    )
    
    cat(
      "Cells:",
      format(
        ncol(obj),
        big.mark = ","
      ),
      "\n"
    )
    
    cat(
      "Genes:",
      format(
        nrow(obj),
        big.mark = ","
      ),
      "\n\n"
    )
    
    cat(
      "Assay:",
      DefaultAssay(obj),
      "\n"
    )
    
    cat(
      "RNA layer: data\n"
    )
    
    cat(
      "Reduction: UMAP\n"
    )
    
  })
  
  
  # ==========================================================
  # Metadata summary
  # ==========================================================
  
  output$metadata_summary <- renderTable({
    
    obj <- current_object()
    
    req(input$metadata)
    
    counts <- table(
      
      obj@meta.data[[input$metadata]],
      
      useNA = "ifany"
      
    )
    
    data.frame(
      
      Group = names(counts),
      
      Cells = as.integer(counts),
      
      row.names = NULL
      
    )
    
  })
  
}


# ============================================================
# 4. Run application
# ============================================================

shinyApp(
  ui = ui,
  server = server
)
