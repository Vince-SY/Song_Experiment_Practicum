
## Set seed
set.seed(12345)

## 1. Create experimental units
data <- data.frame(
  id = 1:30,
  supermarket = rep(c(1, 1, 2, 2, 2), each = 6),
  block = rep(1:5, each = 6)
)

## 2. Blocked treatment assignment
data$treatment <- unlist(
  lapply(1:5, function(b) {
    sample(rep(c(0, 1), each = 3))
  })
)

## 3. Randomize physical positions independently
positions <- expand.grid(
  row = 1:10,
  col = 1:3
)

positions <- positions[
  sample(nrow(positions)),
]

data$row <- positions$row
data$col <- positions$col

## 4. Display layout
layout <- matrix(
  NA_character_,
  nrow = 10,
  ncol = 3
)

for (i in 1:nrow(data)) {
  layout[data$row[i], data$col[i]] <- paste0(
    data$id[i],
    ifelse(data$treatment[i] == 1, "M", "D")
  )
}

print(layout)

## 5. Verify assignments
table(data$block, data$treatment)
table(data$row, data$treatment)
table(data$col, data$treatment)
