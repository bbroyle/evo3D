# ------------------------------------------ #
# EXTENSION_UTILS.R
# utilities for incorporating multiple msa and pdbs
# also special logic for multimers so each codon recieves one patch
# Brad Broyles
# ------------------------------------------ #


# .split_pdb_column() ----

#' split pdb-mapped columns by codon
#'
#' collapses pdb-level residue mappings into codon-level summaries.
#' keeps one reference aa, msa id, and codon string per codon_id,
#' and adds concatenated pdb aa and residue ids for each pdb source.
#'
#' @param df data frame with at least \code{codon_id}, \code{ref_aa},
#'   \code{msa}, \code{codon}, \code{pdb}, \code{pdb_aa}, and \code{residue_id}
#'
#' @return data frame with one row per \code{codon_id}, including:
#' \itemize{
#'   \item \code{ref_aa}, \code{msa}, \code{codon}
#'   \item per-pdb columns \code{<pdb>_pdb_aa}, \code{<pdb>_residue_id}
#'   \item \code{resolved} – number of non-gap residues across pdbs
#' }
#' @keywords internal

.split_pdb_column = function(df) {
  # store the intended codon order
  codon_order <- unique(df$codon_id)

  # keep msa and codon along with codon_id
  base = aggregate(
    list(ref_aa = df$ref_aa,
         msa     = df$msa,
         codon   = df$codon),
    by = list(codon_id = df$codon_id),
    FUN = function(x) x[1]   # take the first entry (all identical within codon_id)
  )

  pdb_ids = unique(df$pdb)

  for (p in pdb_ids) {
    sub = df[df$pdb == p, ]

    aa = aggregate(list(val = sub$pdb_aa),
                   by = list(codon_id = sub$codon_id),
                   FUN = function(x) paste(x, collapse = "+"))
    names(aa)[2] = paste0(p, "_pdb_aa")

    res = aggregate(list(val = sub$residue_id),
                    by = list(codon_id = sub$codon_id),
                    FUN = function(x) paste(x, collapse = "+"))
    names(res)[2] = paste0(p, "_residue_id")

    base = merge(base, aa,  by = "codon_id", all.x = TRUE, sort = FALSE)
    base = merge(base, res, by = "codon_id", all.x = TRUE, sort = FALSE)
  }

  # total occurrences per codon (across all pdbs)
  #occ = aggregate(list(count = df$pdb_aa),
  #                by = list(codon_id = df$codon_id),
  #                FUN = length)

  # number of resolved (non-gap) residues per codon
  res = aggregate(list(resolved = df$pdb_aa),
                  by = list(codon_id = df$codon_id),
                  FUN = function(x) sum(x != "-"))

  #base = merge(base, occ, by = "codon_id", all.x = TRUE, sort = FALSE)
  base = merge(base, res, by = "codon_id", all.x = TRUE, sort = FALSE)

  # finally: restore codon order
  base = base[match(codon_order, base$codon_id), ]

  base
}


# .union_distance() ----

#' union codon sets across patches
#'
#' combines multiple patch strings (codon ids separated by '+')
#' into a single union patch. if a codon appears in multiple
#' patches, its count is taken as the maximum observed across patches.
#'
#' @param patches character vector of patch strings
#'
#' @return single patch string with codon ids joined by '+'
#' @keywords internal

.union_distance = function(patches){

  # drop NA patches
  patches = patches[!is.na(patches)]

  if(length(patches) == 0) return(NA)

  sets = strsplit(patches, '\\+')

  #1_msa1+2_msa1+3_msa1+4_msa1 -- turn to counts 1_msa1(1), 2_msa1(1), 3_msa1(2) -- union will keep the highest count for each codon #
  patch_counts = lapply(sets, table)
  all_codons = unique(unlist(sets))

  mat = sapply(patch_counts, function(tab) {
    out = integer(length(all_codons))
    names(out) = all_codons
    out[names(tab)] = tab
    out
  })

  counts = apply(mat, 1, max)

  # unpack back to patch level string #
  parts = do.call(rbind, strsplit(names(counts), "_"))
  codons = as.numeric(parts[,1])
  msas = parts[,2]

  # order by msa then codon numeric
  ord = order(msas, codons)

  # reconstruct, replicating by count
  ordered_ids = rep(names(counts)[ord], times = counts[ord])

  # final patch string
  patch = paste(ordered_ids, collapse = "+")
}

# .union_exposure_distance() ----

#' union codon sets with exposure filter
#'
#' combines multiple patch strings into a single union patch,
#' keeping the maximum observed count per codon across patches.
#' non valid codons can be filtered out by providing a valid set.
#'
#' @param patches character vector of patch strings (codon ids separated by '+')
#' @param valid_codons character vector of codon ids to retain
#'
#' @return single patch string with codon ids joined by '+'
#' @keywords internal

.union_exposure_distance = function(patches, valid_codons){

  # drop NA patches
  patches = patches[!is.na(patches)]
  if(length(patches) == 0) return(NA)

  sets = strsplit(patches, '\\+')

  #1_msa1+2_msa1+3_msa1+4_msa1 -- turn to counts 1_msa1(1), 2_msa1(1), 3_msa1(2) -- union will keep the highest count for each codon #
  patch_counts = lapply(sets, table)
  all_codons = unique(unlist(sets))

  mat = sapply(patch_counts, function(tab) {
    out = integer(length(all_codons))
    names(out) = all_codons
    out[names(tab)] = tab
    out
  })

  counts = apply(mat, 1, max)

  # unpack back to patch level string #
  parts = do.call(rbind, strsplit(names(counts), "_"))
  codons = as.numeric(parts[,1])
  msas = parts[,2]

  # order by msa then codon numeric
  ord = order(msas, codons)

  # reconstruct, replicating by count
  ordered_ids = rep(names(counts)[ord], times = counts[ord])

  # filter out non valid codons

  # final patch string
  patch = paste(ordered_ids, collapse = "+")
}

# .variable_size_merge() ----

#' merge codon-level patches with variable size rules
#'
#' rebuilds codon patches from residue-level assignments using either
#' geometric distance only, or exposure + distance criteria.
#'
#' @param codon_df data frame of codons with columns \code{codon_id}, \code{exposed}, and \code{codon_patch}
#' @param residue_df data frame of residues with columns \code{codon_id} and \code{codon_patch}
#' @param merge_type character. one of \code{"distance"} (default) or \code{"exposure_distance"}
#' @param only_exposed_in_patch logical. if TRUE, all members must be exposed;
#'   if FALSE, only seed codons must be exposed
#'
#' @return updated \code{codon_df} with merged \code{codon_patch} column
#' @keywords internal
#'
#' @details
#' three modes are supported:
#' \itemize{
#'   \item distance – patches merged by geometry only
#'   \item exposure_distance – seeds must be exposed, members may be buried
#'   \item exposure_distance + only_exposed_in_patch = TRUE – all residues must be exposed
#' }

.variable_size_merge <- function(codon_df, residue_df, merge_type = c("distance", "exposure_distance"),
                                 only_exposed_in_patch = FALSE) {
  merge_type <- match.arg(merge_type)

  codon_ids <- unique(codon_df$codon_id)
  codon_ids <- codon_ids[!is.na(codon_ids)]

  if (merge_type == "distance") {
    # Geometry only: only rebuild when >1 residue
    new_p <- lapply(codon_ids, function(x) {
      ro <- which(residue_df$codon_id == x)
      if (length(ro) > 1) {
        .union_distance(residue_df$codon_patch[ro])
      } else if (length(ro) == 1) {
        residue_df$codon_patch[ro]
      } else {
        NA
      }
    })
    hold <- data.frame(codon_id = codon_ids, codon_patch = unlist(new_p))

  } else if (merge_type == "exposure_distance" && !only_exposed_in_patch) {
    # Allow buried members: same as distance, then filter buried seeds
    valid_codons <- codon_df$codon_id[codon_df$exposed]
    valid_codons <- valid_codons[!is.na(valid_codons)]

    new_p <- lapply(codon_ids, function(x) {
      ro <- which(residue_df$codon_id == x)
      if (length(ro) > 1) {
        .union_distance(residue_df$codon_patch[ro])
      } else if (length(ro) == 1) {
        residue_df$codon_patch[ro]
      } else {
        NA
      }
    })
    hold <- data.frame(codon_id = codon_ids, codon_patch = unlist(new_p))
    # Seed-level exposure filter
    hold$codon_patch[!hold$codon_id %in% valid_codons] <- NA

  } else if (merge_type == "exposure_distance" && only_exposed_in_patch) {
    # All residues (seeds + members) must be exposed, always rebuild
    valid_codons <- codon_df$codon_id[codon_df$exposed]
    valid_codons <- valid_codons[!is.na(valid_codons)]

    new_p <- lapply(codon_ids, function(x) {
      ro <- which(residue_df$codon_id == x)
      if (length(ro) > 0) {
        .union_exposure_distance(residue_df$codon_patch[ro], valid_codons)
      } else {
        NA
      }
    })
    hold <- data.frame(codon_id = codon_ids, codon_patch = unlist(new_p))
    # Explicitly drop buried seeds
    hold$codon_patch[!hold$codon_id %in% valid_codons] <- NA
  }

  # Update codon_df with new patches
  codon_df$codon_patch <- NA
  codon_df$codon_patch[match(hold$codon_id, codon_df$codon_id)] <- hold$codon_patch



  return(codon_df)
}

# .fixed_size_merge() ----

#' merge fixed-size 3d windows
#'
#' builds codon- or residue-level fixed size patches across one or more pdbs.
#' handles multimers and multi-pdb sets by unifying residue neighborhoods into
#' a consistent codon-level environment.
#'
#' @param codon_df data frame of codons with columns \code{codon_id}, \code{exposed}, \code{exposed_count}, and \code{codon_patch}
#' @param residue_df data frame of residues with columns \code{codon_id}, \code{residue_id}, \code{codon_patch}, \code{pdb}, and \code{exposed}
#' @param patch_mode character, either \code{"codon"} or \code{"residue"}
#' @param max_patch integer. maximum patch size (number of codons/residues)
#' @param dist_cutoff numeric. maximum distance in angstroms (default NA = no cutoff)
#' @param merge_type character. either \code{"distance"} or \code{"exposure_distance"}
#' @param pdb_info_sets list of pdb_info objects, each with residue-wise distance matrices
#' @param only_exposed_in_patch logical. if TRUE, all residues in patch must be exposed;
#'   if FALSE, only seed codons must be exposed
#'
#' @return updated \code{codon_df} with rebuilt \code{codon_patch} assignments
#' @keywords internal

.fixed_size_merge = function(codon_df, residue_df, patch_mode,
                             max_patch,
                             dist_cutoff, merge_type,
                             pdb_info_sets, only_exposed_in_patch = TRUE){

  # 8 types of merging [patch_mode C/R, merge_type D/E, only_exposed T/F] #
  # CDT - codon deduplicated, closest distance, exposed w/i pdb context
  # CDF* - codon deduplicated, closest distance, any exposure
  # CET - codon deduplicated, closest distance, exposure from merge context
  # CEF* - codon deduplicated, closest distance, any exposure [works with CDF] #

  # RDT - closest residue, exposed w/i pdb context #
  # RDF* - closest residue, any exposure
  # RET - closest residue, merged exposure
  # REF* - closest residue, any exposure [works with RDF]

  # so six modes - and afterwards we can NA out buried seeds if merge_type E #

  # grab a smaller residue df -- only for those with patches #
  resolved_df = residue_df[!is.na(residue_df$codon_patch),]

  # CDF / CEF ~ exposure no influence #
  if(patch_mode == 'codon' && !only_exposed_in_patch){

    # any closest residue - only count unique codons #
    for(i in seq_len(nrow(codon_df))){
      if(is.na(codon_df$exposed_count[i]) || codon_df$exposed_count[i] == 0) next

      ro = which(resolved_df$codon_id == codon_df$codon_id[i])

      if(length(ro) == 1){
        codon_df$codon_patch[i] = resolved_df$codon_patch[i]
        next
      }

      if(length(ro)>1){

        dist_mats = lapply(ro, function(x){
          id = resolved_df$residue_id[x]
          pid = resolved_df$pdb[x]
          d = pdb_info_sets[[pid]]$residue_dist[id,]

          # drop if above dist_cutoff
          if(!is.na(dist_cutoff)){
            d = d[d <= dist_cutoff]
          }

          # replace ids with codon names #
          map = setNames(
            residue_df$codon_id[residue_df$pdb == pid],
            residue_df$residue_id[residue_df$pdb == pid]
          )

          # replace names
          names(d) = map[names(d)]
          d
        })

        # piece together, sort by dist #
        closest_codons = sort(unlist(dist_mats))

        # deduplicate codons #
        closest_codons = closest_codons[!duplicated(names(closest_codons))]

        # grab either max patch or as many as possible #
        closest_codons = closest_codons[seq_len(min(max_patch, length(closest_codons)))]

        codon_df$codon_patch[i] = paste0(names(closest_codons), collapse = '+')
    }

  }

  }

  # RDF / REF ~ exposure no influence. Use quota to add residues #
  if(patch_mode == 'residue' && !only_exposed_in_patch){

    # any closest residue - only count unique codons #
    for(i in seq_len(nrow(codon_df))){
      if(is.na(codon_df$exposed_count[i]) || codon_df$exposed_count[i] == 0) next

      ro = which(resolved_df$codon_id == codon_df$codon_id[i])

      if(length(ro) == 1){
        codon_df$codon_patch[i] = resolved_df$codon_patch[i]
        next
      }

      if(length(ro)>1){

        dist_mats = lapply(ro, function(x){
          id = resolved_df$residue_id[x]
          pid = resolved_df$pdb[x]
          d = pdb_info_sets[[pid]]$residue_dist[id,]

          # drop if above dist_cutoff
          if(!is.na(dist_cutoff)){
            d = d[d <= dist_cutoff]
          }

          # replace ids with codon names #
          map = setNames(
            residue_df$codon_id[residue_df$pdb == pid],
            residue_df$residue_id[residue_df$pdb == pid]
          )

          # replace names
          names(d) = map[names(d)]
          sort(d)
        })

        # contexts
        contexts = paste0('context', seq_along(dist_mats))

        # add context to dist_mats and unroll so we can sort by overall dist #
        dist_mats = lapply(seq_along(dist_mats), function(x){
          names(dist_mats[[x]]) = paste0('context', x, '_', names(dist_mats[[x]]))
          dist_mats[[x]]
        })



        # piece together, sort by dist #
        context_codons = sort(unlist(dist_mats))

        # and remove context to just have codon
        closest_codons = context_codons
        names(closest_codons) = gsub('context[0-9]+_', '', names(closest_codons))

        patch <- character()
        while(length(closest_codons)) {
          c1 <- closest_codons[1]
          patch <- c(patch, c1)

          for(c in contexts){
            id <- paste0(c, '_', names(c1))
            pos <- match(id, names(context_codons))
            if(!is.na(pos)) {
              closest_codons <- closest_codons[-pos[1]]
              context_codons <- context_codons[-pos[1]]
            }
          }

          if(length(patch) == max_patch) break
        }

        codon_df$codon_patch[i] = paste0(names(patch), collapse = '+')

    }

    }

  }

  # CDT ~ codon deduplicated - exposed with pdb context #
  if(patch_mode == 'codon' && only_exposed_in_patch && merge_type == 'distance'){

    # any closest residue - only count unique codons #
    for(i in seq_len(nrow(codon_df))){
      if(is.na(codon_df$exposed_count[i]) || codon_df$exposed_count[i] == 0) next

      ro = which(resolved_df$codon_id == codon_df$codon_id[i])

      if(length(ro) == 1){
        # good because it was exposed in pdb context #
        codon_df$codon_patch[i] = resolved_df$codon_patch[i]
        next
      }

      if(length(ro)>1){

        dist_mats = lapply(ro, function(x){
          id = resolved_df$residue_id[x]
          pid = resolved_df$pdb[x]
          d = pdb_info_sets[[pid]]$residue_dist[id,]

          # drop if above dist_cutoff
          if(!is.na(dist_cutoff)){
            d = d[d <= dist_cutoff]
          }

          # keep residues that are exposed in their pdb context #
          exp = residue_df$residue_id[residue_df$pdb == pid & residue_df$exposed]
          exp = exp[!is.na(exp)]

          d = d[names(d) %in% exp]

          # replace ids with codon names #
          map = setNames(
            residue_df$codon_id[residue_df$pdb == pid],
            residue_df$residue_id[residue_df$pdb == pid]
          )

          # replace names
          names(d) = map[names(d)]
          d
        })

        # piece together, sort by dist #
        closest_codons = sort(unlist(dist_mats))

        # deduplicate codons #
        closest_codons = closest_codons[!duplicated(names(closest_codons))]

        # grab either max patch or as many as possible #
        closest_codons = closest_codons[seq_len(min(max_patch, length(closest_codons)))]

        codon_df$codon_patch[i] = paste0(names(closest_codons), collapse = '+')
      }

    }

  }

  # RDT ~ exposed withing pdb context - use quota counting #
  if(patch_mode == 'residue' && only_exposed_in_patch && merge_type == 'distance'){

    # any closest residue - only count unique codons #
    for(i in seq_len(nrow(codon_df))){
      if(is.na(codon_df$exposed_count[i]) || codon_df$exposed_count[i] == 0) next

      ro = which(resolved_df$codon_id == codon_df$codon_id[i])

      if(length(ro) == 1){
        codon_df$codon_patch[i] = resolved_df$codon_patch[i]
        next
      }

      if(length(ro)>1){

        dist_mats = lapply(ro, function(x){
          id = resolved_df$residue_id[x]
          pid = resolved_df$pdb[x]
          d = pdb_info_sets[[pid]]$residue_dist[id,]

          # drop if above dist_cutoff
          if(!is.na(dist_cutoff)){
            d = d[d <= dist_cutoff]
          }

          # keep residues that are exposed in their pdb context #
          exp = residue_df$residue_id[residue_df$pdb == pid & residue_df$exposed]
          exp = exp[!is.na(exp)]

          d = d[names(d) %in% exp]

          # replace ids with codon names #
          map = setNames(
            residue_df$codon_id[residue_df$pdb == pid],
            residue_df$residue_id[residue_df$pdb == pid]
          )

          # replace names
          names(d) = map[names(d)]
          sort(d)
        })

        # contexts
        contexts = paste0('context', seq_along(dist_mats))

        # add context to dist_mats and unroll so we can sort by overall dist #
        dist_mats = lapply(seq_along(dist_mats), function(x){
          names(dist_mats[[x]]) = paste0('context', x, '_', names(dist_mats[[x]]))
          dist_mats[[x]]
        })



        # piece together, sort by dist #
        context_codons = sort(unlist(dist_mats))

        # and remove context to just have codon
        closest_codons = context_codons
        names(closest_codons) = gsub('context[0-9]+_', '', names(closest_codons))

        patch <- character()
        while(length(closest_codons)) {
          c1 <- closest_codons[1]
          patch <- c(patch, c1)

          for(c in contexts){
            id <- paste0(c, '_', names(c1))
            pos <- match(id, names(context_codons))
            if(!is.na(pos)) {
              closest_codons <- closest_codons[-pos[1]]
              context_codons <- context_codons[-pos[1]]
            }
          }

          if(length(patch) == max_patch) break
        }

        codon_df$codon_patch[i] = paste0(names(patch), collapse = '+')

      }

    }

  }

  # CET ~ codon deduplicated - exposed in merged table #
  if(patch_mode == 'codon' && only_exposed_in_patch && merge_type == 'exposure_distance'){

    # any closest residue - only count unique codons #
    for(i in seq_len(nrow(codon_df))){
      if(is.na(codon_df$exposed_count[i]) || codon_df$exposed_count[i] == 0) next

      ro = which(resolved_df$codon_id == codon_df$codon_id[i])

      if(length(ro) == 1){
        # good because it was exposed in pdb context #
        codon_df$codon_patch[i] = resolved_df$codon_patch[i]
        next
      }

      if(length(ro)>1){

        dist_mats = lapply(ro, function(x){
          id = resolved_df$residue_id[x]
          pid = resolved_df$pdb[x]
          d = pdb_info_sets[[pid]]$residue_dist[id,]

          # drop if above dist_cutoff
          if(!is.na(dist_cutoff)){
            d = d[d <= dist_cutoff]
          }

          # replace ids with codon names #
          map = setNames(
            residue_df$codon_id[residue_df$pdb == pid],
            residue_df$residue_id[residue_df$pdb == pid]
          )

          # replace names
          names(d) = map[names(d)]

          # filter exposed
          exp = codon_df$codon_id[codon_df$exposed]
          exp = exp[!is.na(exp)]

          d = d[names(d) %in% exp]
        })

        # piece together, sort by dist #
        closest_codons = sort(unlist(dist_mats))

        # deduplicate codons #
        closest_codons = closest_codons[!duplicated(names(closest_codons))]

        # grab either max patch or as many as possible #
        closest_codons = closest_codons[seq_len(min(max_patch, length(closest_codons)))]

        codon_df$codon_patch[i] = paste0(names(closest_codons), collapse = '+')
      }

    }

  }

  # RET ~ exposed in merged table
  if(patch_mode == 'residue' && only_exposed_in_patch && merge_type == 'exposure_distance'){

    # any closest residue - only count unique codons #
    for(i in seq_len(nrow(codon_df))){
      if(is.na(codon_df$exposed_count[i]) || codon_df$exposed_count[i] == 0) next

      ro = which(resolved_df$codon_id == codon_df$codon_id[i])

      if(length(ro) == 1){
        codon_df$codon_patch[i] = resolved_df$codon_patch[i]
        next
      }

      if(length(ro)>1){

        dist_mats = lapply(ro, function(x){
          id = resolved_df$residue_id[x]
          pid = resolved_df$pdb[x]
          d = pdb_info_sets[[pid]]$residue_dist[id,]

          # drop if above dist_cutoff
          if(!is.na(dist_cutoff)){
            d = d[d <= dist_cutoff]
          }

          # replace ids with codon names #
          map = setNames(
            residue_df$codon_id[residue_df$pdb == pid],
            residue_df$residue_id[residue_df$pdb == pid]
          )

          # replace names
          names(d) = map[names(d)]
          sort(d)

          # drop buried
          # keep residues that are exposed in their pdb context #
          exp = codon_df$codon_id[codon_df$exposed]
          exp = exp[!is.na(exp)]

          d = d[names(d) %in% exp]

        })

        # contexts
        contexts = paste0('context', seq_along(dist_mats))

        # add context to dist_mats and unroll so we can sort by overall dist #
        dist_mats = lapply(seq_along(dist_mats), function(x){
          names(dist_mats[[x]]) = paste0('context', x, '_', names(dist_mats[[x]]))
          dist_mats[[x]]
        })



        # piece together, sort by dist #
        context_codons = sort(unlist(dist_mats))

        # and remove context to just have codon
        closest_codons = context_codons
        names(closest_codons) = gsub('context[0-9]+_', '', names(closest_codons))

        patch <- character()
        while(length(closest_codons)) {
          c1 <- closest_codons[1]
          patch <- c(patch, c1)

          for(c in contexts){
            id <- paste0(c, '_', names(c1))
            pos <- match(id, names(context_codons))
            if(!is.na(pos)) {
              closest_codons <- closest_codons[-pos[1]]
              context_codons <- context_codons[-pos[1]]
            }
          }

          if(length(patch) == max_patch) break
        }

        codon_df$codon_patch[i] = paste0(names(patch), collapse = '+')

      }

    }

  }

  # do we need to NA out buried seeds #
  if(merge_type == 'exposure_distance'){
    codon_df$codon_patch[!codon_df$exposed] = NA
  }

  return(codon_df)

}

# collapse_to_codon() ----

#' collapse residue-level patches to codon windows
#'
#' merges residue-level patch assignments across one or more pdbs into a single
#' codon-level window per site. supports variable-length (union) or fixed-size
#' 3d windows, exposure-aware merging, multimers, and multi-pdb contexts.
#'
#' @param residue_df data frame of per-residue annotations (from \code{aln_msa_to_pdb()}),
#'   containing at least: \code{codon_id}, \code{codon}, \code{msa}, \code{pdb},
#'   \code{residue_id}, \code{exposed}, \code{codon_patch}, \code{codon_len},
#'   \code{unique_codon}, \code{max_dist}. rows whose \code{residue_id} starts with
#'   \code{"interface_"} are treated as interface pseudo-patches.
#' @param merge_type character. merging rule for variable/fixed windows:
#'   \code{"distance"} (geometry only) or \code{"exposure_distance"} (exposure + geometry).
#' @param merge_exposure numeric in [0,1]. fraction of pdb contexts in which a codon
#'   must be exposed to qualify as an exposed seed (default 0.5). computed as
#'   \code{exposed_count / resolved}.
#' @param max_patch integer or \code{NA}. if \code{NA}, produce union (variable length);
#'   otherwise build fixed-size windows capped at \code{max_patch}.
#' @param only_exposed_in_patch logical. if \code{TRUE}, all members in a patch must be
#'   exposed; if \code{FALSE}, only the seed must be exposed.
#' @param patch_mode character. \code{"codon"} (deduplicate by codon) or
#'   \code{"residue"} (quota by nearest residues) for fixed-size merging.
#' @param dist_cutoff numeric or \code{NA}. maximum Å distance when assembling fixed-size
#'   neighborhoods (ignored if variable-length union).
#' @param merge_interface_surface logical. when \code{TRUE} and using
#'   \code{merge_type = "exposure_distance"} with \code{only_exposed_in_patch = TRUE},
#'   interface pseudo-patches are filtered to exposed codons before inclusion.
#' @param pdb_info_sets named list of pdb-info objects (as from \code{pdb_to_patch()}),
#'   each containing \code{residue_dist}; required for fixed-size merging across pdbs.
#'
#' @details
#' gap-mapped residues (\code{codon == "-"}) are preserved by assigning a unique
#' placeholder \code{codon_id} and reinserting at their original relative position
#' using nearest upstream/downstream real codons.
#'
#' variable-length mode uses union builders:
#' \itemize{
#'   \item \code{merge_type = "distance"} → union of geometric neighbors
#'   \item \code{merge_type = "exposure_distance"} → union with exposure filtering
#'         (seed-only or all-members via \code{only_exposed_in_patch})
#' }
#' fixed-size mode uses nearest-neighbor ranking across pdb contexts with two schemes:
#' \itemize{
#'   \item \code{patch_mode = "codon"} – deduplicate by codon, keep closest unique codons
#'   \item \code{patch_mode = "residue"} – quota across contexts to balance sources
#' }
#' in exposure-aware merges, buried seeds are nulled post-merge when
#' \code{merge_type = "exposure_distance"}.
#'
#' @return data frame at codon resolution with updated windows, including (where present):
#'   \itemize{
#'     \item \code{codon_id}, \code{ref_aa}, \code{msa}, \code{codon}
#'     \item per-pdb columns \code{<pdb>_pdb_aa}, \code{<pdb>_residue_id}
#'     \item \code{resolved}, \code{exposed_count}, \code{exposed}
#'     \item \code{codon_patch}, \code{codon_len}, \code{unique_codon}, \code{max_dist}
#'     \item \code{msa_subset_id} (set to \code{NA} when \code{codon_patch} is \code{NA})
#'   }
#' interface pseudo-patches are appended as extra rows (with \code{msa_subset_id}
#' formed from interface id and pdb id) when applicable.
#' @export

collapse_to_codon = function(residue_df, merge_type = 'exposure_distance', merge_exposure = 0.5,
                              max_patch, only_exposed_in_patch = TRUE, patch_mode, dist_cutoff,
                              merge_interface_surface, pdb_info_sets){

  # options are distance and exposure_distance #
  # if max_patch is NA it is a union build #
  # if max_patch is set it is a fixed length build #
  # will we need patch_type = residue or codon?
  # if set codons are dedupped already so union will be fine #

  # STEP 0 - store interfaces and gap map residues add later # ----
  ro = grep('^interface', residue_df$residue_id)
  if(length(ro)>0){
    interf_df = residue_df[ro,]
    residue_df = residue_df[-ro,]
  } else {
    interf_df = residue_df[0,,drop=FALSE]
  } # lets try to keep (breaks .split_pdb_column - we will add later) #

  # keep gap mapped residues in their original positions #
  residue_df$upper_context = NA
  residue_df$lower_context = NA

  ro = which(residue_df$codon == '-')
  if(length(ro)>0){
    residue_df$codon_id[ro] = paste0(
      residue_df$msa[ro], '_', residue_df$pdb[ro], '_', residue_df$residue_id[ro]
    )

    for (i in ro) {
      # look upward until you hit a real codon
      up = NA
      j = i - 1
      while (j >= 1 && residue_df$codon[j] == "-") j = j - 1
      if (j >= 1) up = residue_df$codon_id[j]

      # look downward until you hit a real codon
      lo = NA
      j = i + 1
      while (j <= nrow(residue_df) && residue_df$codon[j] == "-") j = j + 1
      if (j <= nrow(residue_df)) lo = residue_df$codon_id[j]

      residue_df$upper_context[i] = up
      residue_df$lower_context[i] = lo
    }
  }

  # split pdb column into individual pdb columns ----
  codon_df = .split_pdb_column(residue_df)

  # add back any gap mapped positions ----
  if(length(ro)>0){
    for (i in rev(ro)) {
      id = residue_df$codon_id[i]
      up = residue_df$upper_context[i]
      lo = residue_df$lower_context[i]

      # extract the codon_df row for this id
      row_idx = which(codon_df$codon_id == id)
      row_data = codon_df[row_idx, , drop = FALSE]

      # drop it from current position
      codon_df = codon_df[-row_idx, ]

      # find new insertion point
      if (!is.na(up) && up %in% codon_df$codon_id) {
        insert_at = which(codon_df$codon_id == up) + 1
      } else if (!is.na(lo) && lo %in% codon_df$codon_id) {
        insert_at = which(codon_df$codon_id == lo)
      } else {
        insert_at = nrow(codon_df) + 1
      }

      # splice it back in
      codon_df = rbind(
        codon_df[seq_len(insert_at - 1), ],
        row_data,
        codon_df[seq(from = insert_at, to = nrow(codon_df)), ]
      )
    }
  }

  residue_df$upper_context = NULL
  residue_df$lower_context = NULL

  # STEP 1 -- update exposure based on merge_exposure threshold ----
  exp_count = aggregate(exposed ~ codon_id, data = residue_df, FUN = sum , na.rm = T)
  codon_df$exposed_count = exp_count$exposed[match(codon_df$codon_id, exp_count$codon_id)]

  # check against threshold (force NA to FALSE after checking against threshold)
  codon_df$exposed = (codon_df$exposed_count / codon_df$resolved) > merge_exposure

  # fix if exposed_count / resolved = 1 then include
  codon_df$exposed = ifelse((codon_df$exposed_count / codon_df$resolved) == 1, TRUE, codon_df$exposed)

  # STEP 3 -- uniting windows (distance, exposure_distance) ----
  # if variable length - do either union (or exposure filtered union) #
  # if fixed length and any position is more than once resolved - need to rebuild patch #
  max_res = max(codon_df$resolved, na.rm = TRUE)

  codon_df$codon_patch = NA
  codon_df$max_dist = NA

  if (max_res <= 1) {
    # no codon has multiple environments - nothing to merge - just copy over info
    codon_df$codon_patch[match(residue_df$codon_id, codon_df$codon_id)] = residue_df$codon_patch
    codon_df$codon_len[match(residue_df$codon_id, codon_df$codon_id)] = residue_df$codon_len
    codon_df$unique_codon[match(residue_df$codon_id, codon_df$codon_id)] = residue_df$unique_codon
    codon_df$max_dist[match(residue_df$codon_id, codon_df$codon_id)] = residue_df$max_dist

  } else {
    if (is.na(max_patch)) {
      codon_df = .variable_size_merge(codon_df = codon_df,
                                     residue_df = residue_df,
                                     merge_type = merge_type,
                                     only_exposed_in_patch = only_exposed_in_patch)
    } else {

      codon_df = .fixed_size_merge(codon_df = codon_df,
                                  residue_df = residue_df,
                                  patch_mode = patch_mode,
                                  dist_cutoff = dist_cutoff,
                                  merge_type = merge_type,
                                  max_patch = max_patch,
                                  pdb_info_sets = pdb_info_sets,
                                  only_exposed_in_patch = only_exposed_in_patch)
    }

  }

  # Step 4 - update msa_subset_id ----
  codon_df$msa_subset_id = codon_df$codon_id
  codon_df$msa_subset_id[is.na(codon_df$codon_patch)] = NA

  # STEP 5 handle interfaces ----
  # handle interfaces  -- skipped by above codon based merges #
  # only if exposure_distance and merge_surface_exposure and only_exposed_in_patch
  if(nrow(interf_df) && merge_type == 'exposure_distance' && merge_interface_surface && only_exposed_in_patch){

    exp = codon_df$codon_id[codon_df$exposed]
    exp = exp[!is.na(exp)]

    for(i in seq_len(nrow(interf_df))){
      if(is.na(interf_df$codon_patch[i])) next
      p = strsplit(interf_df$codon_patch[i], '\\+')[[1]]
      p = p[p %in% exp]
      interf_df$codon_patch[i] = paste0(p, collapse = '+')
    }
  }

  if(nrow(interf_df)){
    # reformat interf_df to match codon_df
    for(i in seq_len(nrow(interf_df))){
      # add new row to codon_df
      newro = nrow(codon_df)+1
      codon_df[newro,] = NA

      codon_df$codon_patch[newro] = interf_df$codon_patch[i]

      pdbid = interf_df$pdb[i]
      codon_df$msa_subset_id[newro] = paste0(interf_df$residue_id[i], '_', pdbid)
      pdbid = paste0(pdbid, '_residue_id')
      codon_df[newro, pdbid] = interf_df$residue_id[i]


    }


  }

  # update codon_len, unique_codon
  codon_df$codon_len[is.na(codon_df$codon_patch)] = NA
  codon_df$unique_codon[is.na(codon_df$codon_patch)] = NA

  for(i in seq_len(nrow(codon_df))){
    if(is.na(codon_df$codon_patch[i])) next

    # load patch split and count #
    p = strsplit(codon_df$codon_patch[i], '\\+')[[1]]
    codon_df$codon_len[i] = length(p)
    codon_df$unique_codon[i] = length(unique(p))

  }

  return(codon_df)
}





