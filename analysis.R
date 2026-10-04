r
# ============================================================
# CAP综合效果评价：熵平衡 + DID + 事件研究
# ============================================================

# 1. 加载包
library(readr)
library(dplyr)
library(WeightIt)
library(cobalt)
library(fixest)
library(tidyr)
library(ggplot2)

# 2. 读取期刊信息表
info <- read_csv("journal_info.csv", locale = locale(encoding = "UTF-8"))

# 3. 读取panel_data文件夹下所有面板数据
panel_files <- list.files("panel_data", pattern = "\\.csv$", full.names = TRUE)
panel_list <- lapply(panel_files, function(f) {
  read_csv(f, locale = locale(encoding = "UTF-8"))
})
panel <- bind_rows(panel_list)

# 4. 统一列名
panel <- panel %>%
  rename(year = Year) %>%
  mutate(journal_id = as.character(journal_id),
         year = as.numeric(year))

info <- info %>%
  mutate(journal_id = as.character(journal_id))

# 5. 合并
df <- panel %>%
  left_join(info, by = "journal_id")

cat("合并后观测数：", nrow(df), "\n")
cat("实验组期刊数：", n_distinct(df$journal_id[df$treat == 1]), "\n")
cat("对照组期刊数：", n_distinct(df$journal_id[df$treat == 0]), "\n")

# 6. 熵平衡（按学科分层）
info <- info %>%
  mutate(ln_tc_pre = log(tc_pre + 1),
         oa_status = factor(oa_status),
         subject = factor(subject))

subjects <- unique(info$subject)
info$weights <- NA

for (s in subjects) {
  sub <- info %>% filter(subject == s)
  if (sum(sub$treat == 1) < 1 | sum(sub$treat == 0) < 1) next
  
  w_out <- weightit(
    treat ~ jif_pre + ln_tc_pre + frequency + oa_status,
    data = sub,
    method = "ebal",
    estimand = "ATT"
  )
  info$weights[info$subject == s] <- w_out$weights
}

# 7. 权重合并到面板
df <- df %>%
  select(-any_of("weights")) %>%
  left_join(info %>% select(journal_id, weights), by = "journal_id")

cat("权重缺失数：", sum(is.na(df$weights)), "\n")

# 8. 基准DID
model_did <- feols(
  CEI ~ treat:post | journal_id + year,
  data = df,
  weights = ~weights,
  cluster = ~journal_id
)
summary(model_did)

# 9. 事件研究
df <- df %>%
  group_by(journal_id) %>%
  mutate(
    rel_time = ifelse(treat == 1, year - cap_year, NA),
    rel_time = ifelse(rel_time < -4, -4, rel_time),
    rel_time = ifelse(rel_time > 4, 4, rel_time)
  ) %>% ungroup()

df <- df %>%
  mutate(rel_factor = factor(rel_time, levels = c(-4:-2, 0:4)))

model_event <- feols(
  CEI ~ i(rel_factor, ref = "-1") | journal_id + year,
  data = df,
  weights = ~weights,
  cluster = ~journal_id
)
summary(model_event)

# 10. 动态效应图
iplot(model_event)

# 11. 缩短出版周期分析（若info中有delta_T列）
if ("delta_T" %in% names(info)) {
  df <- df %>% left_join(info %>% select(journal_id, delta_T), by = "journal_id")
  df <- df %>% mutate(delta_T = ifelse(is.na(delta_T), 0, delta_T))
  
  model_period <- feols(
    CEI ~ treat:post:delta_T | journal_id + year,
    data = df,
    weights = ~weights,
    cluster = ~journal_id
  )
  summary(model_period)
  
  model_immediacy <- feols(
    `Immediacy Index` ~ treat:post:delta_T | journal_id + year,
    data = df,
    weights = ~weights,
    cluster = ~journal_id
  )
  summary(model_immediacy)
}

# 12. 保存合并后数据
write_csv(df, "merged_data_with_weights.csv")