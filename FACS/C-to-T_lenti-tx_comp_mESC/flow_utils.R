inv.log.trans <- function(x) {
    ts <- flowWorkspace::flowjo_log_trans()
    inv <- ts[["inverse"]](x)
    return(inv)
}
log.trans <- function(x) {
    ts <- flowWorkspace::flowjo_log_trans()
    trs <- ts[["transform"]](x)
    return(trs)
}

compute_dynamic_gate <- function(gs, inv.bins) {
    # Dynamic GFP gate will be applied for P1/P2 gate

    df.stat <- data.frame()
    for (min.cutoff in inv.bins) {
        # Create dynamic FL1-A gate
        rg_gfp <- flowCore::rectangleGate(
            "FL1-A" = c(min.cutoff, Inf), "FSC-A" = c(-Inf, Inf),
            filterId = "GFP_dynamic"
        )
        gates <- flowWorkspace::gs_get_pop_paths(gs)
        if ("/P1/P2/GFP_dynamic" %in% gates) {
            flowWorkspace::gs_pop_remove("GFP_dynamic", gs = gs)
            flowWorkspace::gs_pop_add(gs, rg_gfp, parent = "P2")
            flowWorkspace::recompute(gs)
        } else {
            flowWorkspace::gs_pop_add(gs, rg_gfp, parent = "P2")
            flowWorkspace::recompute(gs)
        }
        stat <- flowWorkspace::gs_pop_get_stats(gs, "GFP_dynamic", type = "percent") %>%
            dplyr::mutate(percent, cutoff = inv.log.trans(min.cutoff))
        df.stat <- rbind(df.stat, stat)
    }
    return(df.stat)
}

update_flowFrame_keywords <- function(flowFrame, exprs.m, desc = NULL, data.range = "data") {
    params <- flowCore::parameters(flowFrame)
    pdata <- flowCore::pData(params)

    if (is.null(desc)) {
        desc <- colnames(exprs.m)
    }

    for (i in 1:ncol(flowFrame)) {
        s <- paste("$P", i, "S", sep = "")
        n <- paste("$P", i, "N", sep = "")
        r <- paste("$P", i, "R", sep = "")
        b <- paste("$P", i, "B", sep = "")
        e <- paste("$P", i, "E", sep = "")

        keyval <- list()
        if (!is.na(desc[i])) {
            keyval[[s]] <- desc[i]
        }
        keyval[[n]] <- colnames(exprs.m)[i]

        if (data.range == "data") {
            keyval[[r]] <- ceiling(max(exprs.m[, i], na.rm = TRUE))
        } else if (is.numeric(data.range)) {
            keyval[[r]] <- data.range
        } else {
            stop("Invalid data.range parameter")
        }
        keyval[[b]] <- 32
        keyval[[e]] <- "0,0"
        flowCore::keyword(flowFrame) <- keyval

        pdata[i, "minRange"] <- min(exprs.m[, i], na.rm = TRUE)
        pdata[i, "maxRange"] <- max(exprs.m[, i], na.rm = TRUE)
    }
    flowCore::pData(params) <- pdata
    flowCore::parameters(flowFrame) <- params

    return(flowFrame)
}

copy_keywords <- function(source.frame, target.frame, kw.list) {
    source.keywords <- flowCore::keyword(source.frame)
    for (kw in kw.list) {
        if (!is.null(source.keywords[[kw]])) {
            flowCore::keyword(target.frame) <- source.keywords[kw]
        }
    }
    return(target.frame)
}

as_flowFrame <- function(exprs.m, source.frame = NULL) {
    flow.frame <- flowCore::flowFrame(exprs.m)
    flow.frame <- update_flowFrame_keywords(flow.frame, exprs.m)

    if (!is.null(source.frame)) {
        num.cols <- ncol(flow.frame)
        kw.list <- paste("$P", 1:num.cols, "S", sep = "")
        kw.list <- c(kw.list, paste("$P", 1:num.cols, "N", sep = ""))
        kw.list <- c(kw.list, "$CYT", "$CYTSN", "$DATE", "$FIL", "$BTIM", "$ETIM")
        flow.frame <- copy_keywords(source.frame, flow.frame, kw.list)
        marker.names <- as.character(flowCore::parameters(source.frame)$desc)
        names(marker.names) <- as.character(flowCore::parameters(source.frame)$name)

        # Use the channel name for channels where the description is missing
        w <- is.na(marker.names)
        marker.names[w] <- names(marker.names)[w]
        flowCore::markernames(flow.frame) <- marker.names
    }
    return(flow.frame)
}

concatenate_fcs_files <- function(files.list, output.file = NULL) {
    # premessa::concatenate_fcs_files
    # modified based on the following code:
    # https://rdrr.io/github/ParkerICI/premessa/src/R/fcs_io.R

    m <- lapply(
        files.list,
        function(x) {
            flowCore::read.FCS(filename = x, truncate_max_range = FALSE)
        }
    )

    # Use the first flowFrame as reference
    flow.frame <- m[[1]]
    m <- lapply(m, function(x) {
        flowCore::exprs(x)
    })

    m <- do.call(rbind, m)
    ret <- as_flowFrame(m, flow.frame)
    if (!is.null(output.file)) {
        write_flowFrame(ret, output.file)
    } else {
        return(ret)
    }
}

write_flowFrame <- function(flowFrame, path) {
    f.name <- basename(path)
    flowCore::keyword(flowFrame)[["$FIL"]] <- f.name
    flowCore::write.FCS(flowFrame, path)
    return(invisible(NULL))
}