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
    metadata == "cell_type"
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
# Tissue color helper
# ============================================================

get_tissue_colors <- function(obj) {
  
  tissues <- levels(
    droplevels(
      factor(obj@meta.data[["tissue"]])
    )
  )
  
  colors <- createPalette(
    length(tissues),
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
  
  names(colors) <- tissues
  
  colors
}


# ============================================================
# 2. UI
# ============================================================

ui <- page_sidebar(
  
  title = "Single Cell Explorer for Fliegauf et al. Nature Immunology 2026 - Prinz Laboratory",
  
  fillable = FALSE,
  
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly"
  ) |>
    bs_add_rules("
      .navbar {
        padding-top: 0.35rem;
        padding-bottom: 0.35rem;
        min-height: 48px;
      }
      
      .navbar-brand {
        font-size: 1.1rem;
      }
    "),
  
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
    
    h4("Assess gene expression by"),
    
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
      multiple = FALSE,
      options = list(
        placeholder = "Type a gene name...",
        maxOptions = 1000
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
          card_header("UMAP"),
          
          plotOutput(
            outputId = "umap",
            height = "650px"
          )
        ),
        
        card(
          card_header("Gene expression"),
          
          plotOutput(
            outputId = "feature",
            height = "550px"
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
        
        card_header("Expression by tissue, cell type and condition"),
        
        conditionalPanel(
          
          condition = paste0(
            "input.dataset == ",
            jsonlite::toJSON(
              "Fig 2 scRNAseq of CNS Macrophage Niches",
              auto_unbox = TRUE
            ),
            " && input.metadata == 'tissue'"
          ),
          
          div(
            class = "alert alert-warning",
            style = "margin-bottom: 1rem;",
            
            strong("⚠️ Interpretation warning: "),
            
            "Some cell types are represented by very few cells in certain tissues. ",
            "These observations may reflect contamination or imperfect cell-type assignment ",
            "rather than a biologically meaningful population, e.g. duralCAM in Cortex. ",
            
            "The violin plots and mean-expression bar plots should therefore be ",
            "interpreted with the cell counts in the Dataset information tab in mind."
          )
        ),
        
        layout_columns(
          
          plotOutput(
            outputId = "violin",
            height = "500px"
          ),
          
          plotOutput(
            outputId = "expression_umap",
            height = "300px"
          ),
          
          col_widths = c(8, 4)
        )
      ),
      
      br(),
      
      card(
        card_header("Mean expression"),
        
        plotOutput(
          outputId = "barplot",
          height = "500px"
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
      
      default_metadata <- if ("cell_type" %in% metadata) {
        "cell_type"
      } else {
        metadata[1]
      }
      
      updateSelectInput(
        session = session,
        inputId = "metadata",
        choices = metadata,
        selected = default_metadata
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
      
      genes <- rownames(obj[["RNA"]])
      
      default_gene <- if (
        "Lyve1" %in% genes
      ) {
        "Lyve1"
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
    # Determine legend layout
    # --------------------------------------------------------
    
    is_fig2_subtype <- (
      input$dataset == "Fig 2 scRNAseq of CNS Macrophage Niches" &&
        input$metadata == "cell_type"
    )
    
    if (is_fig2_subtype) {
      
      legend_nrow <- 7
      
    } else {
      
      legend_nrow <- 2
      
    }
    
    legend_ncol <- ceiling(
      length(group_colors) / legend_nrow
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
        legend.direction = "horizontal",
        legend.text = element_text(size = 11),
        legend.key.size = grid::unit(0.7, "cm"),
        legend.box = "vertical"
      )
    
    
    # --------------------------------------------------------
    # Explicit legend layout
    # --------------------------------------------------------
    
    p <- p +
      guides(
        color = guide_legend(
          nrow = legend_nrow,
          ncol = legend_ncol,
          byrow = TRUE,
          override.aes = list(
            size = 5,
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
  # Small cell-type UMAP for Expression tab
  # ==========================================================
  
  output$expression_umap <- renderPlot({
    
    obj <- current_object()
    
    validate(
      need(
        "cell_type" %in% colnames(obj@meta.data),
        "cell_type metadata not available."
      )
    )
    
    cell_type_colors <- get_group_colors(
      obj = obj,
      metadata = "cell_type",
      dataset = input$dataset
    )
    
    DimPlot(
      
      object = obj,
      
      reduction = "umap",
      
      group.by = "cell_type",
      
      cols = cell_type_colors,
      
      label = TRUE,
      
      repel = TRUE,
      
      pt.size = 1,
      
      raster = input$raster,
      
      shuffle = TRUE,
      
      seed = 123
      
    ) +
      
      theme_classic(base_size = 11) +
      
      NoAxes() +
      
      theme(
        legend.position = "none",
        plot.margin = margin(5, 5, 5, 5)
      )
    
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
      
      
      # --------------------------------------------------------
      # Special case: Fig 2 tissue
      #
      # Three-dimensional stratification:
      #
      #   ROWS       = tissue
      #   X-axis     = cell type
      #   VIOLINS    = condition
      #
      # Importantly, each tissue gets ONLY the cell types that
      # actually occur in that tissue.
      # --------------------------------------------------------
      
    } else if (
      input$dataset == "Fig 2 scRNAseq of CNS Macrophage Niches" &&
      input$metadata == "tissue"
    ) {
      
      validate(
        need(
          all(
            c("cell_type", "tissue", "condition") %in%
              colnames(obj@meta.data)
          ),
          "Fig 2 requires cell_type, tissue and condition metadata."
        )
      )
      
      
      # ------------------------------------------------------
      # Extract expression + all three grouping variables
      # ------------------------------------------------------
      
      expression_data <- FetchData(
        object = obj,
        vars = c(
          input$gene,
          "cell_type",
          "tissue",
          "condition"
        ),
        layer = "data"
      )
      
      colnames(expression_data) <- c(
        "expression",
        "cell_type",
        "tissue",
        "condition"
      )
      
      
      # ------------------------------------------------------
      # Remove incomplete metadata
      # ------------------------------------------------------
      
      expression_data <- expression_data %>%
        filter(
          !is.na(expression),
          !is.na(cell_type),
          !is.na(tissue),
          !is.na(condition)
        )
      
      
      # ------------------------------------------------------
      # Get the actual order of tissues in the dataset
      # ------------------------------------------------------
      
      tissue_levels <- unique(
        as.character(expression_data$tissue)
      )
      
      
      # ------------------------------------------------------
      # Get condition order
      # ------------------------------------------------------
      
      condition_levels <- unique(
        as.character(expression_data$condition)
      )
      
      
      # ------------------------------------------------------
      # CRITICAL PART
      #
      # Do NOT globally factor cell_type before faceting.
      #
      # Instead, construct a separate cell_type factor for
      # each tissue containing ONLY the cell types present
      # in that tissue.
      #
      # This prevents empty cell-type positions appearing in
      # tissues where that cell type does not exist.
      # ------------------------------------------------------
      
      expression_data <- expression_data %>%
        mutate(
          tissue = factor(
            tissue,
            levels = tissue_levels
          ),
          condition = factor(
            condition,
            levels = condition_levels
          )
        ) %>%
        group_by(tissue) %>%
        mutate(
          cell_type_plot = factor(
            cell_type,
            levels = unique(as.character(cell_type))
          )
        ) %>%
        ungroup()
      
      
      # ------------------------------------------------------
      # Condition colors
      # ------------------------------------------------------
      
      condition_colors <- get_condition_colors(obj)
      
      
      # Make sure the colors are named according to the
      # actual condition levels in this dataset.
      condition_colors <- condition_colors[
        condition_levels
      ]
      
      
      # ------------------------------------------------------
      # Plot
      # ------------------------------------------------------
      
      ggplot(
        expression_data,
        aes(
          x = cell_type_plot,
          y = expression,
          fill = condition
        )
      ) +
        
        geom_violin(
          position = position_dodge(
            width = 0.8
          ),
          scale = "width",
          trim = TRUE,
          na.rm = TRUE
        ) +
        
        facet_grid(
          rows = vars(tissue),
          scales = "free_x",
          space = "free_x",
          drop = TRUE
        ) +
        
        scale_fill_manual(
          values = condition_colors,
          drop = FALSE
        ) +
        
        labs(
          title = input$gene,
          subtitle = "Expression by tissue, cell type and condition",
          x = "Cell type",
          y = "Normalized expression",
          fill = "Condition"
        ) +
        
        theme_classic(
          base_size = 13
        ) +
        
        theme(
          axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            vjust = 1
          ),
          
          strip.text.y = element_text(
            angle = 0,
            face = "bold"
          ),
          
          strip.background = element_blank(),
          
          legend.position = "bottom",
          
          panel.spacing.y = unit(
            0.8,
            "lines"
          )
        )
      
      
      # --------------------------------------------------------
      # Default behavior for all other metadata
      # --------------------------------------------------------
      
    } else {
      
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
      
      
      # --------------------------------------------------------
      # Special case: Fig 2 combined tissue × cell type × condition
      #
      # X-axis  = cell type
      # Fill    = condition
      # Rows    = tissue
      #
      # Only combinations that actually exist in the data
      # are shown. Missing cell types in a tissue therefore
      # remain absent rather than producing empty bars.
      # --------------------------------------------------------
      
    } else if (
      input$dataset == "Fig 2 scRNAseq of CNS Macrophage Niches" &&
      input$metadata == "tissue"
    ) {
      
      condition_colors <- get_condition_colors(obj)
      
      expression_data <- FetchData(
        object = obj,
        vars = c(
          input$gene,
          "cell_type",
          "tissue",
          "condition"
        ),
        layer = "data"
      )
      
      colnames(expression_data) <- c(
        "expression",
        "cell_type",
        "tissue",
        "condition"
      )
      
      # Remove cells with missing metadata
      expression_data <- expression_data %>%
        filter(
          !is.na(expression),
          !is.na(cell_type),
          !is.na(tissue),
          !is.na(condition)
        )
      
      # ------------------------------------------------------
      # Preserve the actual ordering in the dataset
      # ------------------------------------------------------
      
      tissue_levels <- levels(
        droplevels(
          factor(expression_data$tissue)
        )
      )
      
      cell_type_levels <- levels(
        droplevels(
          factor(expression_data$cell_type)
        )
      )
      
      condition_levels <- levels(
        droplevels(
          factor(expression_data$condition)
        )
      )
      
      expression_data <- expression_data %>%
        mutate(
          tissue = factor(
            tissue,
            levels = tissue_levels
          ),
          cell_type = factor(
            cell_type,
            levels = cell_type_levels
          ),
          condition = factor(
            condition,
            levels = condition_levels
          )
        )
      
      # ------------------------------------------------------
      # Calculate mean expression for each actual combination
      #
      # Because we group only observed combinations, cell types
      # absent from a tissue do NOT generate zero/empty bars.
      # ------------------------------------------------------
      
      summary_data <- expression_data %>%
        group_by(
          tissue,
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
      
      # ------------------------------------------------------
      # Plot
      # ------------------------------------------------------
      
      ggplot(
        summary_data,
        aes(
          x = cell_type,
          y = mean_expression,
          fill = condition
        )
      ) +
        
        geom_col(
          position = position_dodge(
            width = 0.8
          ),
          width = 0.7
        ) +
        
        facet_grid(
          tissue ~ .,
          scales = "free_x",
          space = "free_x",
          drop = TRUE
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
          subtitle = "Expression by tissue, cell type and condition",
          x = "Cell type",
          y = "Mean normalized expression",
          fill = "Condition"
        ) +
        
        theme_classic(
          base_size = 13
        ) +
        
        theme(
          axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            vjust = 1
          ),
          strip.text.y = element_text(
            angle = 0,
            face = "bold"
          ),
          strip.background = element_blank(),
          legend.position = "bottom",
          panel.spacing = unit(
            0.8,
            "lines"
          )
        )
      
      
      # --------------------------------------------------------
      # Default behavior for all other metadata
      # --------------------------------------------------------
      
    } else {
      
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
    
    # --------------------------------------------------------
    # General dataset information
    # --------------------------------------------------------
    
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
    
    
    # --------------------------------------------------------
    # Fig 2 + tissue
    # Show cell type × tissue
    # --------------------------------------------------------
    
    if (
      input$dataset == "Fig 2 scRNAseq of CNS Macrophage Niches" &&
      input$metadata == "tissue" &&
      "tissue" %in% colnames(obj@meta.data) &&
      "cell_type" %in% colnames(obj@meta.data)
    ) {
      
      cat(
        "\n\n"
      )
      
      cat(
        "Cells by cell type and tissue:\n\n"
      )
      
      counts <- table(
        obj@meta.data$cell_type,
        obj@meta.data$tissue,
        useNA = "ifany"
      )
      
      print(counts)
      
    } else if (
      
      # ------------------------------------------------------
      # Condition
      # Show cell type × condition
      # ------------------------------------------------------
      
      input$metadata == "condition" &&
      "condition" %in% colnames(obj@meta.data) &&
      "cell_type" %in% colnames(obj@meta.data)
      
    ) {
      
      cat(
        "\n\n"
      )
      
      cat(
        "Cells by cell type and condition:\n\n"
      )
      
      counts <- table(
        obj@meta.data$cell_type,
        obj@meta.data$condition,
        useNA = "ifany"
      )
      
      print(counts)
    }
    
  })
  
  
  # ==========================================================
  # Selected metadata
  # ==========================================================
  
  output$metadata_summary <- renderTable({
    
    obj <- current_object()
    
    req(input$metadata)
    
    
    # --------------------------------------------------------
    # Fig 2 + tissue
    # Show tissue × cell type
    # --------------------------------------------------------
    
    if (
      input$dataset == "Fig 2 scRNAseq of CNS Macrophage Niches" &&
      input$metadata == "tissue" &&
      "tissue" %in% colnames(obj@meta.data) &&
      "cell_type" %in% colnames(obj@meta.data)
    ) {
      
      counts <- table(
        obj@meta.data$cell_type,
        obj@meta.data$tissue,
        useNA = "ifany"
      )
      
      result <- as.data.frame.matrix(counts)
      
      result <- cbind(
        cell_type = rownames(result),
        result,
        row.names = NULL
      )
      
      return(result)
    }
    
    
    # --------------------------------------------------------
    # Condition
    # Show condition × cell type
    # --------------------------------------------------------
    
    if (
      input$metadata == "condition" &&
      "condition" %in% colnames(obj@meta.data) &&
      "cell_type" %in% colnames(obj@meta.data)
    ) {
      
      counts <- table(
        obj@meta.data$cell_type,
        obj@meta.data$condition,
        useNA = "ifany"
      )
      
      result <- as.data.frame.matrix(counts)
      
      result <- cbind(
        cell_type = rownames(result),
        result,
        row.names = NULL
      )
      
      return(result)
    }
    
    
    # --------------------------------------------------------
    # Default behavior
    # --------------------------------------------------------
    
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
