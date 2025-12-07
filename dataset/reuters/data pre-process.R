# install.packages(c("tm", "XML"))
library(tm)
library(XML)

getwd()
# setwd("D:/000-202526/BS/project/reference/data reference/reuters+21578+text+categorization+collection")
files <- list.files("reuters21578", pattern="sgm$", full.names=TRUE)

extract_doc <- function(text){
  newid <- sub('.*NEWID="([0-9]+)".*', '\\1', text)
  title <- sub('.*<TITLE>(.*?)</TITLE>.*', '\\1', text)
  body  <- sub('.*<BODY>(.*?)</BODY>.*', '\\1', text)
  
  topics <- gregexpr("<D>(.*?)</D>", text)
  topics <- regmatches(text, topics)[[1]]
  topics <- gsub("</?D>", "", topics)
  topics <- paste(topics, collapse = ",")
  
  data.frame(
    id = newid,
    title = ifelse(grepl("<TITLE>", text), title, NA),
    body = ifelse(grepl("<BODY>", text), body, NA),
    topics = ifelse(length(topics) == 0, NA, topics),
    stringsAsFactors = FALSE
  )
}

all_docs <- list()

for (f in files) {
  raw <- paste(readLines(f, warn=FALSE), collapse="\n")
  docs <- unlist(strsplit(raw, "</REUTERS>"))
  
  for (d in docs){
    if (grepl("<REUTERS", d)){
      all_docs[[length(all_docs)+1]] <- extract_doc(d)
    }
  }
}

reuters_df <- do.call(rbind, all_docs)

dim(reuters_df)
head(reuters_df)


library(tm)

# --- View 1: body ---
corp_body <- VCorpus(VectorSource(reuters_df$body))

corp_body <- tm_map(corp_body, content_transformer(tolower))
corp_body <- tm_map(corp_body, removePunctuation)
corp_body <- tm_map(corp_body, removeNumbers)
corp_body <- tm_map(corp_body, removeWords, stopwords("en"))
corp_body <- tm_map(corp_body, stripWhitespace)

dtm_body <- DocumentTermMatrix(
  corp_body,
  control = list(
    wordLengths = c(3, Inf),
    bounds = list(global = c(5, Inf)) # delete really fewer words
  )
)

X_body <- as.matrix(dtm_body)  # N x V1

# --- View 2: title ---
corp_title <- VCorpus(VectorSource(reuters_df$title))

corp_title <- tm_map(corp_title, content_transformer(tolower))
corp_title <- tm_map(corp_title, removePunctuation)
corp_title <- tm_map(corp_title, removeNumbers)
corp_title <- tm_map(corp_title, removeWords, stopwords("en"))
corp_title <- tm_map(corp_title, stripWhitespace)

dtm_title <- DocumentTermMatrix(
  corp_title,
  control = list(
    wordLengths = c(3, Inf),
    bounds = list(global = c(3, Inf))
  )
)

X_title <- as.matrix(dtm_title)  # N x V2

# split topics string
topic_list <- strsplit(reuters_df$topics, ",")
topic_vocab <- sort(unique(unlist(topic_list)))
T <- length(topic_vocab)

N <- nrow(reuters_df)
X_topics <- matrix(0L, nrow = N, ncol = T)
colnames(X_topics) <- topic_vocab

for (i in seq_len(N)) {
  ts <- topic_list[[i]]
  ts <- ts[ts != "" & !is.na(ts)]  # delete the "space"
  if (length(ts) > 0) {
    X_topics[i, match(ts, topic_vocab)] <- 1L
  }
}

views <- list(
  body   = X_body,   # View 1
  title  = X_title,  # View 2
  topics = X_topics  # View 3
)

# use views[[v]][i, ] visit Pitman-Yor multiview inference


## install.packages("poweRlaw")

library(poweRlaw)

## Based on tail exponent alpha: PYP discount parameter d ---------------
## just experience, not seriously 
##   - smaller alpha  → much heavier tailed → bigger d
suggest_d_from_alpha <- function(alpha) {
  if (is.na(alpha)) return(c(NA_real_, NA_real_))
  
  if (alpha <= 1.2) {
    d_range <- c(0.7, 0.9)   # extremely heavy-tail
  } else if (alpha <= 1.6) {
    d_range <- c(0.4, 0.7)   # obviously heavy-tail
  } else if (alpha <= 2.0) {
    d_range <- c(0.2, 0.5)   # medium heavy-tail
  } else {
    d_range <- c(0.05, 0.2)  # relatively slight, use PYP or DP
  }
  d_range
}

## for these three view, analysize their heavy-tail ------------------------
## parameters:  
##   X_body, X_title, X_topics: files × feature matrix (count or 0/1);
##   view_names: for figures and tables,
##   n_boot: bootstrap times，for p-value test.
library(poweRlaw)

suggest_d_from_alpha <- function(alpha) {
  if (alpha <= 1.2) return(c(0.7, 0.9))
  if (alpha <= 1.6) return(c(0.4, 0.7))
  if (alpha <= 2.0) return(c(0.2, 0.5))
  return(c(0.05, 0.2))
}

analyze_powerlaw_views <- function(
    X_body, X_title, X_topics,
    view_names = c("body", "title", "topics"),
    xmin_fixed = 5)        # <-- key point：force xmin ≥ 5
{
  X_list <- list(X_body, X_title, X_topics)
  results <- list()
  
  par(mfrow = c(1, 3))
  
  for (i in seq_along(X_list)) {
    X <- X_list[[i]]
    name <- view_names[i]
    
    freq <- colSums(X)
    freq <- sort(freq[freq >= xmin_fixed], decreasing = TRUE)  # <--- big enough
    
    if (length(freq) < 20) {
      warning(paste("View", name, "features too few after filtering"))
      next
    }
    
    rank <- seq_along(freq)
    
    plot(rank, freq, log="xy",
         main = paste("View:", name),
         xlab="Rank (log)", ylab="Freq (log)")
    
    # Fit power-law with forced xmin
    pl <- displ$new(freq)
    pl$setXmin(xmin_fixed)
    est <- estimate_pars(pl)
    pl$setPars(est$pars)
    
    lines(pl, col=2)
    
    alpha_hat <- est$pars
    
    # Try lognormal comparison safely
    lr_stat <- NA; p_value <- NA; lr_sign <- NA
    
    try({
      ln <- dislnorm$new(freq)
      ln$setXmin(xmin_fixed)
      est_ln <- estimate_pars(ln)
      ln$setPars(est_ln$pars)
      
      cmp <- compare_distributions(pl, ln)
      lr_stat <- cmp$test_statistic
      p_value <- cmp$p_two_sided
      lr_sign <- ifelse(lr_stat > 0, 1, -1)
    }, silent = TRUE)
    
    d_range <- suggest_d_from_alpha(alpha_hat)
    
    results[[i]] <- data.frame(
      view = name,
      n_features = length(freq),
      xmin = xmin_fixed,
      alpha_hat = alpha_hat,
      LR = lr_stat,
      LR_sign = lr_sign,
      p_value = p_value,
      d_min = d_range[1],
      d_max = d_range[2]
    )
  }
  
  par(mfrow = c(1, 1))
  
  return(do.call(rbind, results))
}


res <- analyze_powerlaw_views(X_body, X_title, X_topics)
print(res)

