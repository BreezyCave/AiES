# Copyright 2025 National Center of Neurology and Psychiatry
# BSD 3-Clause License (see LICENSE file)
#
# Internal helpers shared by axDistmap(), axSvm() and axQnt().
# Nothing in this file is exported.
#
# Rationale: the image-processing pipeline used to be duplicated in axDistmap()
# and axQnt(). Keeping a single implementation here consolidates the segmentation process.


# ---------------------------------------------------------------------------
# Reference imaging scale
# ---------------------------------------------------------------------------
# The classifiers built with this package were trained on images acquired at
# 0.51190453 um/px (1392 x 1040 px = 712.57 x 532.38 um, 20x objective), which
# were then resized to 900 px width. The effective scale of the training data
# is therefore:
#
#   0.51190453 * 1392 / 900 = 0.791746 um/px
#
# Every pixel-unit parameter of the pipeline is defined at this scale:
#
#   subBack  = 30 px^2   ->  18.8 um^2  (objects below ~4.9 um equivalent diameter)
#   thresh   = 8, 8      ->  17 x 17 px = 13.5 x 13.5 um  (w/h are half-widths)
#   haralick = 1, 2 px   ->  0.79 / 1.58 um  (fixed in px; not rescaled)
#
# Resampling an image to .aiesRefPixelSize therefore restores all of them at once,
# enabling built-in models to be transferable across cameras. At any other working
# scale, subBack and the threshold window are rescaled by .aiesScaleParams() to keep
# their physical size, but the Haralick offsets stay at 1 and 2 px, so their physical
# span changes.

.aiesRefPixelSize <- 0.791746      # um/px, effective scale of the training data
.aiesRefSubBack   <- 30            # px^2 at .aiesRefPixelSize
.aiesRefThreshWH  <- 8             # half-width, px, at .aiesRefPixelSize
.aiesRefHaralick  <- c(1L, 2L)     # GLCM offsets, px, at .aiesRefPixelSize


# ---------------------------------------------------------------------------
# Feature name mapping
# ---------------------------------------------------------------------------
# Short name (used by the built-in models and in every exported file)
# -> original EBImage column name
#
# .aiesSegment() renames every feature to its short name, so that axDistmap(),
# axSvm() and axQnt() all exchange short names only. "Area" (EBImage's s.area)
# is the area weight used for AII / DI; it is never a predictor.

.aiesNameMap <- c(
    "Area" = "s.area",        "Per"  = "s.perimeter",
    "Ram"  = "s.radius.mean", "Ras"  = "s.radius.sd",
    "Ran"  = "s.radius.min",  "Rax"  = "s.radius.max",
    "ASM"  = "h.asm.s1", "CON" = "h.con.s1", "COR" = "h.cor.s1",
    "VAR"  = "h.var.s1", "IDM" = "h.idm.s1", "SAV" = "h.sav.s1",
    "SVA"  = "h.sva.s1", "SEN" = "h.sen.s1", "ENT" = "h.ent.s1",
    "DVA"  = "h.dva.s1", "DEN" = "h.den.s1", "F12" = "h.f12.s1",
    "F13"  = "h.f13.s1",
    "SM2"  = "h.asm.s2", "ON2" = "h.con.s2", "OR2" = "h.cor.s2",
    "AR2"  = "h.var.s2", "DM2" = "h.idm.s2", "AV2" = "h.sav.s2",
    "SV2"  = "h.sva.s2", "SE2" = "h.sen.s2", "NT2" = "h.ent.s2",
    "VA2"  = "h.dva.s2", "EN2" = "h.den.s2", "F22" = "h.f12.s2",
    "F32"  = "h.f13.s2",
    "mcx"  = "m.cx",     "mcy" = "m.cy",
    "mjx"  = "m.majoraxis", "ecc" = "m.eccentricity",
    "the"  = "m.theta",     "cir" = "Cir"
)

## original column name -> short name
.aiesRevNameMap <- stats::setNames(names(.aiesNameMap), .aiesNameMap)

## Columns that are never used as predictors by default (axSvm, model = "custom"
## without `features`):
##   Group, s.area, Area, PixelSize  grouping variable, area weight, scale record
##   mcx, mcy, the                   centroid coordinates and orientation: they describe
##                                   where and how an object lies in the image, not its state
##   DVA, VA2                        identical to CON, ON2: EBImage computes h.dva with the
##                                   same C code as h.con (src/haralick.c)
.aiesNonFeature <- c("Group", "s.area", "Area", "PixelSize",
                     "mcx", "mcy", "the", "DVA", "VA2")


# ---------------------------------------------------------------------------
# Built-in model presets
# ---------------------------------------------------------------------------
# The two models reported in the accompanying manuscript (Table 1).
# NOTE: the label "g0.1" / "g0.2" refers to the gamma held fixed *during forward
# selection*; the gamma recorded here is the final value after the local grid
# search, which is why the two differ.

.aiesModels <- list(
    "g0.1_APA" = list(
        label           = "g0.1 APA",
        features        = c("Ras", "ENT", "Per", "SVA", "SV2"),
        gamma           = 0.0871,
        cost            = 1.189,
        targetPixelSize = .aiesRefPixelSize,
        file            = "g0.1_APA.svm"
    ),
    "g0.2_ACC" = list(
        label           = "g0.2 ACC",
        features        = c("CON", "ecc", "mjx", "cir", "Ras", "Ran"),
        gamma           = 0.207,
        cost            = 1.110,
        targetPixelSize = .aiesRefPixelSize,
        file            = "g0.2_ACC.svm"
    )
)

## Union of every feature required by the built-in models, as short names,
## plus Area (area weight).
.aiesStandardColumns <- c(
    "Area", "Per", "Ras", "Ran",
    "ecc", "mjx",
    "ENT", "SVA", "SV2", "CON",
    "cir"
)


# ---------------------------------------------------------------------------
# Name conversion helpers
# ---------------------------------------------------------------------------

## Rename original EBImage column names to the short names used by the models.
## Columns without a mapping (e.g. Group) are left untouched, and columns that
## already carry short names are left untouched as well.
.aiesToShort <- function(df) {
    cn <- colnames(df)
    hit <- cn %in% names(.aiesRevNameMap)
    cn[hit] <- unname(.aiesRevNameMap[cn[hit]])
    colnames(df) <- cn
    df
}

## Accept a feature vector given in either convention and return short names.
.aiesAsShort <- function(features) {
    out <- features
    hit <- out %in% names(.aiesRevNameMap)
    out[hit] <- unname(.aiesRevNameMap[out[hit]])
    out
}

## The predictors an already-trained svm object expects, as short names.
.aiesModelFeatures <- function(svm_model) {
    ft <- NULL
    if (!is.null(svm_model$terms)) ft <- attr(svm_model$terms, "term.labels")
    if (is.null(ft) && !is.null(svm_model$SV)) ft <- colnames(svm_model$SV)
    ft
}


# ---------------------------------------------------------------------------
# Physical scale
# ---------------------------------------------------------------------------

## Extract one attribute value from an XML string (first occurrence).
.aiesXmlAttr <- function(xml, attr, numeric = TRUE) {
    pat <- paste0("[[:space:]]", attr, "=\"[^\"]*\"")
    m <- regmatches(xml, regexpr(pat, xml))
    if (length(m) == 0) return(NA)
    v <- sub(paste0("^[[:space:]]", attr, "=\"([^\"]*)\"$"), "\\1", m)
    if (numeric) suppressWarnings(as.numeric(v)) else v
}


## Read the image metadata that the pipeline needs, without reading pixels.
## Returns a list with
##   pixelSize        um/px, or NULL if it cannot be determined reliably
##   pixelSource      where pixelSize came from
##   containerBits    bits per sample of the TIFF container (8, 16, ...)
##   significantBits  bits actually used by the camera (OME SignificantBits), or NULL
##
## Sources, in order of preference:
##   1. OME-XML in the ImageDescription tag (Nikon NIS-Elements, Bio-Formats, ...)
##   2. TIFF resolution tags with a physical unit (inch / cm)
##   3. ImageJ: resolution tag without unit plus 'unit=micron' in the description
##      (ImageJ saves micrometer calibration this way; XResolution = px per um)
## Guessing is avoided: when none applies, pixelSize is NULL and the user is
## asked to supply `pixelSize`.
.aiesTiffInfo <- function(file) {
    out <- list(pixelSize = NULL, pixelSource = NULL,
                containerBits = NULL, significantBits = NULL)
    if (!requireNamespace("tiff", quietly = TRUE)) return(out)

    ## suppressWarnings: libtiff reports vendor-private tags (e.g. ImageJ 50838/50839) as unknown
    info <- tryCatch(suppressWarnings(tiff::readTIFF(file, payload = FALSE)),
                     error = function(e) NULL)
    if (is.null(info)) return(out)

    if (!is.null(info$bits.per.sample)) {
        out$containerBits <- as.integer(info$bits.per.sample)[1]
    }

    desc <- info$description
    desc <- if (is.null(desc) || length(desc) == 0) "" else paste(desc, collapse = "\n")

    ## 1. OME-XML
    if (grepl("<OME", desc, fixed = TRUE)) {
        px <- .aiesXmlAttr(desc, "PhysicalSizeX")
        if (!is.na(px) && px > 0) {
            unit <- .aiesXmlAttr(desc, "PhysicalSizeXUnit", numeric = FALSE)
            ## OME default unit is micrometer; anything not recognized is treated as such
            mult <- if (is.na(unit)) 1 else switch(unit, nm = 1e-3, mm = 1e3, 1)
            out$pixelSize   <- px * mult
            out$pixelSource <- "OME-XML"
        }
        sb <- .aiesXmlAttr(desc, "SignificantBits")
        if (!is.na(sb) && sb > 0) out$significantBits <- as.integer(sb)
    }

    ## 2./3. Resolution tags
    if (is.null(out$pixelSize)) {
        xres <- info$x.resolution
        unit <- info$resolution.unit
        if (!is.null(xres) && is.finite(xres) && xres > 0) {
            if (!is.null(unit) && unit %in% c(3, "3", "cm", "centimeter")) {
                out$pixelSize <- 10000 / xres; out$pixelSource <- "TIFF resolution tag"
            } else if (!is.null(unit) && unit %in% c(2, "2", "inch")) {
                out$pixelSize <- 25400 / xres; out$pixelSource <- "TIFF resolution tag"
            } else if (grepl("unit=[[:space:]]*(micron|um|\u00b5m)", desc)) {
                out$pixelSize <- 1 / xres;     out$pixelSource <- "ImageJ description"
            }
        }
    }
    out
}

## Bring data that use fewer bits than their container (e.g. 12-bit camera data
## saved as 16-bit TIFF) to the full [0, 1] range.
##
## Why this matters: the first step of the pipeline is a sigmoid with a fixed
## midpoint of 0.5. The built-in models were trained on 8-bit images, which
## EBImage maps onto the full [0, 1] range. 12-bit data in a 16-bit container
## are read as values below 0.0625, so the sigmoid operates on its extreme tail
## and the segmentation changes completely, even for an identical image.
.aiesRescaleBits <- function(img, bitDepth, containerBits) {
    if (is.null(bitDepth) || is.null(containerBits)) return(img)
    if (bitDepth >= containerBits) return(img)

    f <- (2^containerBits - 1) / (2^bitDepth - 1)
    img <- img * f
    over <- mean(img > 1)
    if (over > 0) {
        message(sprintf("  Note: %.2f%% of pixels exceeded the %d-bit range and were clipped.",
                        100 * over, bitDepth))
        img[img > 1] <- 1
    }
    message(sprintf("  Rescaled %d-bit data in a %d-bit container (x%.4f).",
                    bitDepth, containerBits, f))
    img
}


## Decide the resize width and the effective pixel size for one image.
## Returns a list; effPixelSize is NA when the physical scale is unknown, in
## which case parameter auto-scaling is disabled (k = 1, legacy behavior).
## `pixelSize` is the input pixel size already resolved by the caller (argument or
## file metadata); `pixelSource` says which, and is used in the message only.
.aiesResolveScale <- function(file, widthPx,
                              pixelSize = NULL,
                              targetPixelSize = .aiesRefPixelSize,
                              resizeW = 900,
                              pixelSource = "argument") {

    ## 1. validate the input pixel size
    src <- pixelSource
    if (!is.null(pixelSize) && (!is.finite(pixelSize) || pixelSize <= 0)) pixelSize <- NULL

    ## 2. decide the resize width
    if (!is.null(targetPixelSize) && !is.null(pixelSize)) {
        resizeW_eff  <- max(1L, as.integer(round(widthPx * pixelSize / targetPixelSize)))
        effPixelSize <- targetPixelSize
        modeType <- "scale-matched"
        mode <- sprintf("scale-matched (%s: %.5f um/px -> %.5f um/px)",
                        src, pixelSize, targetPixelSize)
    } else {
        if (is.null(resizeW)) {
            resizeW_eff  <- NULL
            effPixelSize <- if (is.null(pixelSize)) NA_real_ else pixelSize
            modeType <- "no-resize"
            mode <- "native resolution"
        } else {
            resizeW_eff  <- as.integer(resizeW)
            effPixelSize <- if (is.null(pixelSize)) NA_real_ else pixelSize * widthPx / resizeW
            modeType <- "fixed-width"
            mode <- sprintf("fixed width (%d px)", resizeW_eff)
        }
    }

    ## Report the fact when the physical scale is unknown (processing continues)
    if (!is.null(targetPixelSize) && is.null(pixelSize)) {
        warning(sprintf(paste0(
            "%s: pixel size unknown; no scale matching. Processed as %s; ",
            "pixel-unit parameters are not rescaled."),
            basename(file), mode), call. = FALSE)
    }

    ## 3. parameter scaling factor relative to the reference scale
    k <- if (is.na(effPixelSize)) 1 else .aiesRefPixelSize / effPixelSize

    list(resizeW = resizeW_eff, effPixelSize = effPixelSize, k = k, mode = mode,
         modeType = modeType)
}


## Pixel-unit parameters are rescaled to match the working resolution. `subBack`
## is treated as an area in reference-scale pixels, so it is multiplied by k².
## The local threshold half-width (reference: 8 px) represents a length. Thus,
## it is multiplied by k, rounded to the nearest integer, and cannot be set by
## the user. In contrast, the Haralick GLCM scales are deliberately not rescaled
## and always remain at c(1, 2) px, regardless of the working resolution.
.aiesScaleParams <- function(k, subBack = .aiesRefSubBack, autoParam = TRUE) {
    if (!isTRUE(autoParam) || !is.finite(k) || k == 1) {
        return(list(subBack = subBack,
                    threshWH = .aiesRefThreshWH,
                    haralickScales = .aiesRefHaralick))
    }
    list(subBack        = max(1, round(subBack * k^2)),
         threshWH       = max(1L, as.integer(round(.aiesRefThreshWH * k))),
         haralickScales = .aiesRefHaralick)
}


# ---------------------------------------------------------------------------
# Color handling
# ---------------------------------------------------------------------------
# A TIFF may arrive as a 2D grayscale array, a 3-channel RGB image (e.g. from a
# color camera), or a multi-frame stack. Reduce any of these to one grayscale
# channel, or return NULL if the requested colorMode is inconsistent with the file.
#
#   gray : single channel as is; for multi-channel input, channel 1 is used
#          (correct when the channels duplicate one grayscale image)
#   rgb  : the three channels are averaged

.aiesToSingleChannel <- function(img, mode, label = "") {
    d <- dim(img)
    nch <- if (length(d) == 2L) 1L else d[3]

    if (nch == 1L) {
        if (mode != "gray") {
            message(sprintf("  Note: single-channel image; colorMode = '%s' has no effect.", mode))
        }
        out <- if (length(d) == 2L) img else img[, , 1]
        return(EBImage::Image(out, colormode = "Grayscale"))
    }

    out <- switch(mode,
        gray = {
            message(sprintf("  Note: %d channels/frames present; using the first (colorMode = 'gray').", nch))
            img[, , 1]
        },
        rgb = {
            if (nch != 3L) {
                warning(sprintf("%s: skipped. colorMode = 'rgb' requires 3 channels, but %d found.",
                                label, nch), call. = FALSE)
                return(NULL)
            }
            (img[, , 1] + img[, , 2] + img[, , 3]) / 3
        }
    )
    EBImage::Image(out, colormode = "Grayscale")
}


## Read a file, reduce it to one channel and resample it to the working scale.
## Returns NULL when the file should be skipped.
.aiesReadGray <- function(file, colorMode = "gray",
                          pixelSize = NULL,
                          targetPixelSize = .aiesRefPixelSize,
                          resizeW = 900) {

    tmpImage <- EBImage::readImage(file, all = TRUE)
    dimImg <- dim(tmpImage)
    nChannel <- if (length(dimImg) == 2L) 1L else dimImg[3]
    meta <- .aiesTiffInfo(file)
    message(sprintf("  dim: [%s] | channels/frames: %d | container: %s-bit",
                    paste(dimImg, collapse = " x "), nChannel,
                    if (is.null(meta$containerBits)) "?" else meta$containerBits))

    ## Binary check before any rescaling or interpolation
    if (all(tmpImage %in% c(0, 1))) {
        message(sprintf("  Skipping binary image: %s", file))
        return(NULL)
    }

    ## Bit depth: taken from the OME SignificantBits when present. There is no
    ## user override, so that it cannot conflict with the file metadata. Without
    ## it, the data are assumed to span the full range of the container.
    bd <- meta$significantBits
    bitRescaled <- !is.null(bd) && !is.null(meta$containerBits) && bd < meta$containerBits
    if (!is.null(bd)) {
        tmpImage <- .aiesRescaleBits(tmpImage, bd, meta$containerBits)
    } else if (max(tmpImage) < 4096 / 65536) {
        warning(sprintf(paste0(
            "%s: maximum intensity is only %.4f of the full range. ",
            "If this is 12-bit camera data stored unscaled in a 16-bit TIFF, export it as ",
            "OME-TIFF or scale it to the full range in the acquisition software; ",
            "otherwise the segmentation will not match the built-in models."),
            basename(file), max(tmpImage)), call. = FALSE)
    }

    ## Pixel size: explicit argument > file metadata > resizeW fallback.
    ## Metadata formats differ between vendors, so the argument is preferred.
    ## When both exist, the argument is used and both values are reported.
    argPixelSize <- pixelSize          # as given by the user, kept for the parameter record
    pixSrc <- "argument"
    if (!is.null(meta$pixelSize)) {
        if (is.null(pixelSize)) {
            pixelSize <- meta$pixelSize
            pixSrc <- meta$pixelSource
            message(sprintf("  pixel size from %s: %.5f um/px", meta$pixelSource, pixelSize))
        } else {
            warning(sprintf(paste0(
                "%s: using pixelSize = %s um/px (argument); ",
                "file metadata (%s) gives %.5f um/px. Check that they agree."),
                basename(file), format(pixelSize), meta$pixelSource, meta$pixelSize),
                call. = FALSE)
        }
    }

    scale <- .aiesResolveScale(file, dimImg[1], pixelSize, targetPixelSize, resizeW, pixSrc)

    ## Reduce channels first, so discarded channels are never resampled
    grayImage <- .aiesToSingleChannel(tmpImage, colorMode, basename(file))
    invisible({rm(tmpImage); gc()})
    if (is.null(grayImage)) return(NULL)

    if (!is.null(scale$resizeW)) {
        grayImage <- EBImage::resize(grayImage, w = scale$resizeW)
    }
    message(sprintf("  %s | working size: [%s]",
                    scale$mode, paste(dim(grayImage), collapse = " x ")))

    attr(grayImage, "aiesScale") <- scale
    ## Facts recorded in the *_Parameters.csv file (see .aiesParamRow)
    attr(grayImage, "aiesInfo") <- list(
        widthPx = dimImg[1], heightPx = dimImg[2], channels = nChannel,
        containerBits = meta$containerBits, significantBits = bd,
        bitRescaled = bitRescaled,
        argPixelSize = argPixelSize,
        metaPixelSize = meta$pixelSize, metaSource = meta$pixelSource,
        workingWidthPx = dim(grayImage)[1])
    grayImage
}


# ---------------------------------------------------------------------------
# Parameter record (*_Parameters.csv)
# ---------------------------------------------------------------------------
# One row per processed image. PixelSize is the effective scale (um/px) of the
# image that was actually analyzed; axSvm() reads it (matched through DataFile)
# to store the training scale with the model.

.aiesParamColumns <- c("FileName", "DataFile", "Date", "AiESVersion",
    "WidthPx", "HeightPx", "Channels", "ContainerBits", "SignificantBits",
    "BitRescaled", "ColorMode", "ArgPixelSize", "MetaPixelSize", "MetaSource",
    "PixelSize", "ScaleMode", "ResizeW", "SubBack", "ThreshWH", "FeatureSet")

## File name relative to its input folder (`roots` may hold several candidates).
.aiesRelPath <- function(file, roots) {
    f <- normalizePath(as.character(file), winslash = "/", mustWork = FALSE)
    for (r in roots) {
        rr <- sub("/+$", "", normalizePath(as.character(r), winslash = "/", mustWork = FALSE))
        if (startsWith(f, paste0(rr, "/"))) {
            return(substr(f, nchar(rr) + 2L, nchar(f)))
        }
    }
    basename(f)
}

.aiesParamRow <- function(fileName, dataFile, sdate, info, colorMode,
                          scale, prm, featureSet) {
    na0 <- function(x) if (is.null(x) || length(x) == 0) NA else x
    ver <- tryCatch(as.character(utils::packageVersion("AiES")),
                    error = function(e) NA_character_)
    data.frame(
        FileName        = fileName,
        DataFile        = na0(dataFile),
        Date            = as.character(sdate),
        AiESVersion     = ver,
        WidthPx         = info$widthPx,
        HeightPx        = info$heightPx,
        Channels        = info$channels,
        ContainerBits   = na0(info$containerBits),
        SignificantBits = na0(info$significantBits),
        BitRescaled     = info$bitRescaled,
        ColorMode       = colorMode,
        ArgPixelSize    = na0(info$argPixelSize),
        MetaPixelSize   = na0(info$metaPixelSize),
        MetaSource      = na0(info$metaSource),
        PixelSize       = scale$effPixelSize,
        ScaleMode       = scale$modeType,
        ResizeW         = info$workingWidthPx,
        SubBack         = prm$subBack,
        ThreshWH        = prm$threshWH,
        FeatureSet      = featureSet,
        stringsAsFactors = FALSE)[, .aiesParamColumns]
}


# ---------------------------------------------------------------------------
# Segmentation and feature extraction
# ---------------------------------------------------------------------------
# Single source of truth for the processing pipeline. Returns a list with the
# labelled segmentation, the normalized distance map, and the feature table
# (all EBImage features plus circularity, under their short names).
#
# Haralick scales are currently always c(1, 2) (see .aiesScaleParams). The
# renaming block below only takes effect if other scales are ever passed; it
# keeps the columns named .s1 / .s2 so that downstream name mapping stays valid.

.aiesSegment <- function(grayImage, subBack = 30, threshWH = 8,
                         haralickScales = c(1L, 2L)) {

    contrastST <- (1 / (1 + (0.5 / grayImage)^5))
    medFltr    <- EBImage::medianFilter(contrastST, 1)
    logTrsf    <- (log1p(medFltr) / log1p(max(medFltr)))
    logTrsfLog <- (log1p(logTrsf) / log1p(max(logTrsf)))
    thrWBInv   <- 1 - EBImage::thresh(logTrsfLog, threshWH, threshWH, -0.04)
    dm      <- EBImage::distmap(thrWBInv)
    nDm     <- EBImage::normalize(dm)
    bnrySeg <- EBImage::bwlabel(nDm)

    ## Remove small objects that image processing cannot classify
    sdat  <- as.data.frame(EBImage::computeFeatures.shape(bnrySeg))
    rmNum <- which(sdat$s.area <= subBack)
    bnrySeg <- EBImage::rmObjects(bnrySeg, rmNum)

    sdat <- as.data.frame(EBImage::computeFeatures.shape(bnrySeg))
    mdat <- as.data.frame(EBImage::computeFeatures.moment(bnrySeg))
    hdat <- as.data.frame(EBImage::computeFeatures.haralick(
        bnrySeg, nDm, haralick.scales = haralickScales))

    ## Restore canonical .s1 / .s2 suffixes when non-default offsets were used.
    ## An intermediate ".S" marker avoids collisions (e.g. scales c(2,3)).
    if (!identical(as.integer(haralickScales), c(1L, 2L))) {
        cn <- colnames(hdat)
        for (i in seq_along(haralickScales)) {
            cn <- sub(paste0("\\.s", haralickScales[i], "$"), paste0(".S", i), cn)
        }
        colnames(hdat) <- sub("\\.S", ".s", cn)
    }

    singleData <- data.frame(cbind(sdat, hdat, mdat))
    invisible({rm(list = c("sdat", "hdat", "mdat")); gc(); gc()})

    singleData <- stats::na.omit(singleData)
    singleData$Cir <- singleData$s.area * pi * 4 / singleData$s.perimeter^2

    ## Short names from here on (s.area -> Area, Cir -> cir, ...)
    singleData <- .aiesToShort(singleData)

    list(features = singleData, segmentation = bnrySeg, distmap = nDm)
}
