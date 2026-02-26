## code to prepare `demo_data` dataset goes here

# 生成模拟微生物组示例数据
set.seed(42)

# 分类信息表
tax_table <- data.frame(
  OTU     = paste0("OTU", 1:30),
  Kingdom = "Bacteria",
  Phylum  = rep(c("Firmicutes","Bacteroidetes","Proteobacteria",
                  "Actinobacteria","Verrucomicrobia"), each = 6),
  Class   = paste0("Class", rep(1:10, each = 3)),
  Order   = paste0("Order", rep(1:15, each = 2)),
  Family  = paste0("Family", 1:30),
  Genus   = paste0("Genus", 1:30),
  stringsAsFactors = FALSE
)

# 丰度表（样本 × OTU）
n_samples <- 10
otu_matrix <- matrix(
  rpois(n_samples * 30, lambda = rep(c(500,300,200,100,50), each = 6)),
  nrow  = n_samples,
  ncol  = 30,
  dimnames = list(
    paste0("Sample", 1:n_samples),
    paste0("OTU", 1:30)
  )
)

# 转换为相对丰度
otu_rel <- otu_matrix / rowSums(otu_matrix)

# 样本元数据
sample_data <- data.frame(
  SampleID = paste0("Sample", 1:n_samples),
  Group    = rep(c("Control", "Treatment"), each = 5),
  Subject  = paste0("Sub", 1:n_samples),
  stringsAsFactors = FALSE
)

# 打包成列表
demo_data <- list(
  otu_table   = otu_rel,
  tax_table   = tax_table,
  sample_data = sample_data
)

usethis::use_data(demo_data, overwrite = TRUE)
