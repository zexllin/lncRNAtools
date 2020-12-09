#####################################
######    plot tow y axis    ########
ggplot2.two_y_axis <- function(g1, g2) {
  g1 <- ggplotGrob(g1)
  g2 <- ggplotGrob(g2)

  # Get the location of the plot panel in g1.
  # These are used later when transformed elements of g2 are put back into g1
  pp <- c(subset(g1$layout, name == 'panel', se = t:r))

  # Overlap panel for second plot on that of the first plot
  g1 <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == 'panel')]], pp$t, pp$l, pp$b, pp$l)

  # Then proceed as before:

  # ggplot contains many labels that are themselves complex grob;
  # usually a text grob surrounded by margins.
  # When moving the grobs from, say, the left to the right of a plot,
  # Make sure the margins and the justifications are swapped around.
  # The function below does the swapping.
  # Taken from the cowplot package:
  # https://github.com/wilkelab/cowplot/blob/master/R/switch_axis.R

  hinvert_title_grob <- function(grob){

    # Swap the widths
    widths <- grob$widths
    grob$widths[1] <- widths[3]
    grob$widths[3] <- widths[1]
    grob$vp[[1]]$layout$widths[1] <- widths[3]
    grob$vp[[1]]$layout$widths[3] <- widths[1]

    # Fix the justification
    grob$children[[1]]$hjust <- 1 - grob$children[[1]]$hjust
    grob$children[[1]]$vjust <- 1 - grob$children[[1]]$vjust
    grob$children[[1]]$x <- unit(1, 'npc') - grob$children[[1]]$x
    grob
  }

  # Get the y axis title from g2
  index <- which(g2$layout$name == 'ylab-l') # Which grob contains the y axis title?
  ylab <- g2$grobs[[index]]        # Extract that grob
  ylab <- hinvert_title_grob(ylab)     # Swap margins and fix justifications

  # Put the transformed label on the right side of g1
  g1 <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pp$r)
  g1 <- gtable_add_grob(g1, ylab, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'ylab-r')

  # Get the y axis from g2 (axis line, tick marks, and tick mark labels)
  index <- which(g2$layout$name == 'axis-l')  # Which grob
  yaxis <- g2$grobs[[index]]          # Extract the grob

  # yaxis is a complex of grobs containing the axis line, the tick marks, and the tick mark labels.
  # The relevant grobs are contained in axis$children:
  #   axis$children[[1]] contains the axis line;
  #   axis$children[[2]] contains the tick marks and tick mark labels.

  # First, move the axis line to the left
  yaxis$children[[1]]$x <- unit.c(unit(0, 'npc'), unit(0, 'npc'))

  # Second, swap tick marks and tick mark labels
  ticks <- yaxis$children[[2]]
  ticks$widths <- rev(ticks$widths)
  ticks$grobs <- rev(ticks$grobs)

  # Third, move the tick marks
  ticks$grobs[[1]]$x <- ticks$grobs[[1]]$x - unit(1, 'npc') + unit(3, 'pt')

  # Fourth, swap margins and fix justifications for the tick mark labels
  ticks$grobs[[2]] <- hinvert_title_grob(ticks$grobs[[2]])

  # Fifth, put ticks back into yaxis
  yaxis$children[[2]] <- ticks

  # Put the transformed yaxis on the right side of g1
  g1 <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pp$r)
  g1 <- gtable_add_grob(g1, yaxis, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'axis-r')
  grid.newpage()
  grid.draw(g1)
}


#####################################
######     ID conversion       ######

ensembl2symbolFun <- function(ensemblID, info='symbol') {
  geneInfo <- biotype[match(ensemblID, biotype$ensemblID),]
  geneSymbol <- geneInfo$geneSymbol

  if (info=='symbol') {
    return (geneSymbol)
  } else if (info=='all') {
    return (geneInfo)
  }
}


symbol2ensemblFun <- function(symbol, info='ensemblID') {
  geneInfo <- biotype[match(symbol, biotype$geneSymbol),]
  ensemblID <- geneInfo$ensemblID

  if (info=='ensemblID') {
    return (ensemblID)
  } else if (info=='all') {
    return (geneInfo)
  }
}


organizeEnrichFun <- function(go) {

  Terms <- paste(go$ID, go$Description, sep='~')
  Counts <- go$Count

  GeneRatio <- go$GeneRatio
  BgRatio <- go$BgRatio

  pValue <- go$pvalue
  FDR <- go$p.adjust

  listTotal <- vapply(go$GeneRatio, function(v)
    convertRatioFun(v, type='bg'), numeric(1))
  popHits <- vapply(go$BgRatio, function(v)
    convertRatioFun(v, type='hit'), numeric(1))
  popTotal <- vapply(go$BgRatio, function(v)
    convertRatioFun(v, type='bg'), numeric(1))

  foldEnrichment <- as.vector(Counts/listTotal*popTotal/popHits)

  geneID <- go$geneID
  geneSymbol <- unlist(lapply(strsplit(geneID, '/', fixed=TRUE),
                              function(v) paste(ensembl2symbolFun(v), collapse = '/')))

  goOutput <- data.frame(Terms, Counts, GeneRatio, BgRatio, pValue, FDR,
                         foldEnrichment, geneID, geneSymbol)

  return (goOutput)
}


###
convertRatioFun <- function(v, type='bg') {
  ratio <- strsplit(v, '/', fixed=TRUE)

  if (type=='bg') {
    num <- as.numeric(as.character(ratio[[1]][2]))
  } else if (type=='hit') {
    num <- as.numeric(as.character(ratio[[1]][1]))
  }

  return (num)
}
